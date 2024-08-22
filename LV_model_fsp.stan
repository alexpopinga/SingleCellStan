functions {
  // The propensity function for a birth-death model
  vector bursting_gene_propensity(real[] x, real[] theta) {
    real birth = theta[1];
    real chomp = theta[2];
    real death = theta[3];
    
    vector[3] propensities;
    propensities[1] = birth * x[1];  // R1: 0 -> x
    propensities[1] = chomp * x[1] * x[2];  // R2: x -> y
    propensities[2] = death * x[2]; // R3: y -> 0
    
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
        vector[3] propensities = bursting_gene_propensity(state, theta);
        
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
  real y0[2]; // Initial conditions for prey and predators
  int<lower=0> y[T, 2]; // Observed data (prey and predators)
  int nStates; // Number of states in FSP approximation
  matrix[2, 3] S; // Stoichiometry matrix
  matrix[2, nStates] states; // States in FSP approximation
}
parameters {
  real<lower=0> birth; // Rate of 0 -> prey
  real<lower=0> chomp; // Rate of prey -> predator
  real<lower=0> death; // Rate of predator -> 0
  real<lower=0> sigma; // Measurement noise
}
transformed parameters {
  real theta[3];
  theta[1] = birth;
  theta[2] = chomp;
  theta[3] = death;

  // Build the infinitesimal generator matrix
  matrix[nStates + 1, nStates + 1] infGen = build_inf_gen(S, theta, states);
}
model {
  // Priors using lognormal distributions
  birth ~ lognormal(0, 1); // Prior for 0 -> prey
  chomp ~ lognormal(0, 1); // Prior for prey -> predator
  death ~ lognormal(0, 1); // Prior for predator -> 0
  sigma ~ lognormal(0, 1); // Prior for measurement noise

  // Likelihood
  for (t in 1:T) {
    vector[nStates + 1] p = solve_fsp(infGen, ts[t]);
    
    // Find the marginal probability for the observed number of protein molecules
    int predator_count = y[t, 2];
    real p_obs = 0.0;
    
    for (i in 1:nStates) {
      if (states[2, i] == predator_count) {
        p_obs += p[i];
      }
    }
    
    // Use the marginal probability as the likelihood
    target += log(p_obs + 1e-10); // Add a small number to avoid log(0)
  }
}
