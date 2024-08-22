library(rstan)

# Define the stoichiometry matrix
S <- matrix(c(
  1, -1, 0, 
  0, 1, -1
), nrow = 2, byrow = TRUE)

# Define the states in the FSP approximation
# This is a simplification for illustrative purposes; you may need to generate this based on your specific problem.
states <- expand.grid(prey = 0:100, predator = 0:100)
states <- as.matrix(states[states$prey + states$predator == 1, ])

# Example data
time_points <- c(0.1, 10, 20, 30, 40, 50) # Use a small positive value for the initial time
observed_data <- matrix(c(
  100, 10,
  80, 20,
  70, 30,
  80, 20,
  75, 25,
  85, 15
), ncol = 2, byrow = TRUE)
initial_conditions <- c(100, 10)  # Initial OFF, ON, Protein

data <- list(
  T = length(time_points),
  t0 = 0.0,
  ts = time_points,
  y0 = initial_conditions,
  y = observed_data,
  nStates = nrow(states),
  S = S,
  states = t(states) # Transpose to match the Stan format
)

# Custom initialization function
init_fun <- function() {
  list(
    birth = exp(rnorm(1, 0, 1)),  # Lognormal(0, 1)
    chomp = exp(rnorm(1, 0, 1)),   # Lognormal(0, 1)
    death = exp(rnorm(1, 0, 1)),  # Lognormal(0, 1)
    sigma = exp(rnorm(1, 0, 1))   # Lognormal(0, 1)
  )
}

# Compile and fit the model
stan_model <- stan_model(file = 'LV_model_fsp.stan')
fit <- sampling(stan_model, data = data, init = init_fun, iter = 400, warmup = 200, chains = 4, verbose = TRUE)

# Print the results
print(fit)

# Plot trace plots to check convergence
traceplot(fit)

# Display summary statistics
summary_stats <- summary(fit)$summary
print(summary_stats)

# Save the plot of the chains to a file
library(ggplot2)
ggsave("traceplot.png", traceplot(fit))

# Additional diagnostics
check_hmc_diagnostics(fit)

# Print the results
print(fit)

# Check the sampling algorithm used
stan_args <- fit@stan_args
print(stan_args)