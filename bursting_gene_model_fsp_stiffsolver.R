library(rstan)

# Define the stoichiometry matrix
S <- matrix(c(
  -1, 1, 0, 0,
   1, -1, 0, 0,
   0, 0, 1, -1
), nrow = 3, byrow = TRUE)

# Define the states in the FSP approximation
states <- expand.grid(OFF = 0:1, ON = 0:1, Protein = 0:180)
states <- as.matrix(states[states$OFF + states$ON <= 1, ])

# Example data
time_points <- c(0.1, 10, 20, 30, 40, 50) # Use a small positive value for the initial time
observed_data <- matrix(c(
  1, 0, 0,
  0, 1, 20,
  1, 0, 40,
  0, 1, 50,
  1, 0, 20,
  1, 0, 15
), ncol = 3, byrow = TRUE)
initial_conditions <- c(1, 0, 0)  # Initial OFF, ON, Protein

theta <- c(0.02, 0.05, 10, 0.1)

data <- list(
  T = length(time_points),
  t0 = 0.0,
  ts = time_points,
  y0 = initial_conditions,
  y = observed_data,
  nStates = nrow(states), 
  S = S,
  states = t(states), # Transpose to match the Stan format
  theta = theta
)

# Custom initialization function
init_fun <- function() {
  list(
    kon = exp(rnorm(1, 0, 1)),  # Lognormal(0, 1)
    koff = exp(rnorm(1, 0, 1)),   # Lognormal(0, 1)
    kP = exp(rnorm(1, 0, 1)),  # Lognormal(0, 1)
    gam = exp(rnorm(1, 0, 1)),   # Lognormal(0, 1)
    sigma = exp(rnorm(1, 0, 1))   # Lognormal(0, 1)
  )
}

# Compile and fit the model
stan_model <- stan_model(file = 'bursting_gene_model_fsp_stiffsolver.stan')
fit <- sampling(stan_model, data = data, init = init_fun, iter = 2000, warmup = 1000, chains = 4, verbose = TRUE)

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

# Check the sampling algorithm used
stan_args <- fit@stan_args
print(stan_args)
