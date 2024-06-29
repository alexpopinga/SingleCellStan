functions {
  real[] bursting_gene_odes(real t,
                            real[] y,
                            real[] theta,
                            real[] x_r,
                            int[] x_i) {
    real dydt[3];
    real k_on = theta[1];
    real k_off = theta[2];
    real k_m = theta[3];
    real gamma_m = theta[4];
    real k_p = theta[5];
    real gamma_p = theta[6];

    real G_on = y[1];
    real m = y[2];
    real p = y[3];
    real G_off = 1.0 - G_on;

    dydt[1] = k_on * G_off - k_off * G_on; // d[G_on]/dt
    dydt[2] = k_m * G_on - gamma_m * m;    // dm/dt
    dydt[3] = k_p * m - gamma_p * p;       // dp/dt

    return dydt;
  }
}

data {
  int<lower=1> N;           // Number of time points
  real t0;                  // Initial time
  real ts[N];               // Time points
  real y0[3];               // Initial conditions for G_on, m, and p
  int y_obs[N, 2];          // Observed data for m and p
}

parameters {
  real<lower=0> k_on;
  real<lower=0> k_off;
  real<lower=0> k_m;
  real<lower=0> gamma_m;
  real<lower=0> k_p;
  real<lower=0> gamma_p;
}

transformed parameters {
  real y_hat[N, 3];          // Predicted concentrations
  real theta[6];
  theta[1] = k_on;
  theta[2] = k_off;
  theta[3] = k_m;
  theta[4] = gamma_m;
  theta[5] = k_p;
  theta[6] = gamma_p;
  y_hat = integrate_ode_rk45(bursting_gene_odes, y0, t0, ts, theta, rep_array(0.0, 0), rep_array(0, 0));
}

model {
  // Priors
  k_on ~ lognormal(0, 1);
  k_off ~ lognormal(0, 1);
  k_m ~ lognormal(0, 1);
  gamma_m ~ lognormal(0, 1);
  k_p ~ lognormal(0, 1);
  gamma_p ~ lognormal(0, 1);

  // Likelihood
  for (i in 1:N) {
    y_obs[i, 1] ~ poisson(y_hat[i, 2]); // mRNA counts
    y_obs[i, 2] ~ poisson(y_hat[i, 3]); // Protein counts
  }
}
