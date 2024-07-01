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

  // Build the infinitesimal generator matrix for a continuous-time Markov chain
  matrix build_inf_gen(matrix S, real[] theta, matrix states) {
    int nSpecies = rows(states);
    int nStates = cols(states);
    int nReactions = cols(S);
    
    matrix[nStates + 1, nStates + 1] infGen;
    matrix[nStates, nStates] sink;
    
    infGen = rep_matrix(0, nStates + 1, nStates + 1);
    sink = rep_matrix(0, nStates, nStates);
    
    // Compute the propensities for all states
    for (mu in 1:nReactions) {
      for (i in 1:nStates) {
        real state[nSpecies] = to_array_1d(states[, i]);
        vector[4] propensities = bursting_gene_propensity(state, theta);
        
        // Flow of probability out of all states due to reaction mu
        infGen[i, i] -= propensities[mu];
        
        // Compute the states after reaction mu
        real newState[nSpecies];
        for (j in 1:nSpecies) {
          newState[j] = state[j] + S[j, mu];
        }
        
        // Check if the new state is non-negative
        if (min(newState) >= 0) {
          // Find the index of the new state
          int idx = 0;
          for (k in 1:nStates) {
            int match = 1;
            for (j in 1:nSpecies) {
              if (newState[j] != states[j, k]) {
                match = 0;
                break;
              }
            }
            if (match == 1) {
              idx = k;
              break;
            }
          }
          
          if (idx > 0) {
            infGen[idx, i] += propensities[mu];
          } else {
            sink[i, i] += propensities[mu];
          }
        } else {
          sink[i, i] += propensities[mu];
        }
      }
    }
    
    // Add the sink as the final row of the infinitesimal generator
    for (i in 1:nStates) {
      infGen[nStates + 1, i] = sink[i, i];
    }
    
    return infGen;
  }
  
  // Solve the CME using the FSP method
  vector solve_fsp(matrix infGen, real t) {
    int nStates = rows(infGen) - 1;
    vector[nStates + 1] p;
    
    // Initialize the state probabilities (initial state is known)
    p = rep_vector(0, nStates + 1);
    p[1] = 1.0; // Assuming initial state is the first state
    
    // Matrix exponential for the solution (simplified approach)
    matrix[nStates + 1, nStates + 1] expGen = matrix_exp(infGen * t);
    p = expGen * p;
    
    return p;
  }
}
data {
  int<lower=1> T; // Number of time points
  real t0; // Initial time
  real ts[T]; // Observation times
  real y0[3]; // Initial conditions for OFF, ON, and Protein
  int<lower=0> y[T, 3]; // Observed data (OFF, ON, and Protein)
  int nStates; // Number of states in FSP approximation
  matrix[3, 4] S; // Stoichiometry matrix
  matrix[3, nStates] states; // States in FSP approximation
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

  // Build the infinitesimal generator matrix
  matrix[nStates + 1, nStates + 1] infGen = build_inf_gen(S, theta, states);
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
    vector[nStates + 1] p = solve_fsp(infGen, ts[t]);
    
    // Find the marginal probability for the observed number of protein molecules
    int protein_count = y[t, 3];
    real p_obs = 0.0;
    
    for (i in 1:nStates) {
      if (states[3, i] == protein_count) {
        p_obs += p[i];
      }
    }
    
    // Use the marginal probability as the likelihood
    target += log(p_obs + 1e-10); // Add a small number to avoid log(0)
  }
}
