functions {
  // The propensity function for the bursting gene model
  vector bursting_gene_propensity(real[] x, real[] theta) {
    real kon = theta[1];
    real koff = theta[2];
    real kP = theta[3];
    real gam = theta[4];
    
    vector[4] propensities;
    propensities[1] = kon * x[1];  // R1: OFF -> ON
    propensities[2] = koff * x[2]; // R2: ON -> OFF
    propensities[3] = kP * x[2];   // R3: ON -> ON + Protein
    propensities[4] = gam * x[3];  // R4: Protein -> null
    
    return propensities;
  }

  // The ODE function for the bursting gene model
  real[] bursting_gene_odes(real t, real[] y, real[] theta, real[] x_r, int[] x_i) {
    real dydt[3];
    
    // Stoichiometry matrix
    real S[3, 4] = {
      {-1,  1,  0,  0},
      { 1, -1,  0,  0},
      { 0,  0,  1, -1}
    };
    
    // Calculate propensities
    vector[4] w = bursting_gene_propensity(y, theta);
    
    // Initialize derivatives
    dydt[1] = 0;
    dydt[2] = 0;
    dydt[3] = 0;
    
    // Compute derivatives using stoichiometry matrix and propensities
    for (i in 1:4) {
      dydt[1] += S[1, i] * w[i];
      dydt[2] += S[2, i] * w[i];
      dydt[3] += S[3, i] * w[i];
    }
    
    return dydt;
  }
}
data {
  int<lower=1> T; // Number of time points
  real ts[T]; // Observation times
  real y0[3]; // Initial conditions for OFF, ON, and Protein
  int<lower=0> y[T, 3]; // Observed data (OFF, ON, and Protein)
}
parameters {
  real<lower=0> kon; // Rate of OFF -> ON
  real<lower=0> koff; // Rate of ON -> OFF
  real<lower=0> kP; // Rate of ON -> ON + Protein
  real<lower=0> gam; // Rate of Protein -> null
  real<lower=0> sigma; // Measurement noise
}
transformed parameters {
  real theta[4];
  theta[1] = kon;
  theta[2] = koff;
  theta[3] = kP;
  theta[4] = gam;

  // Solve the ODE
  real y_hat[T, 3];
  y_hat = integrate_ode_rk45(bursting_gene_odes, y0, 0, ts, theta, rep_array(0.0, 0), rep_array(0, 0));
}
model {
  // Priors using lognormal distributions
  kon ~ lognormal(0, 1); // Prior for OFF -> ON rate
  koff ~ lognormal(0, 1); // Prior for ON -> OFF rate
  kP ~ lognormal(0, 1); // Prior for ON -> ON + Protein rate
  gam ~ lognormal(0, 1); // Prior for Protein -> null rate
  sigma ~ lognormal(0, 1); // Prior for measurement noise

  // Likelihood
  for (t in 1:T) {
    y[t] ~ normal(y_hat[t], sigma);
  }
}