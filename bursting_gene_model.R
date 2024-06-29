library(rstan)

# Example data
time_points <- c(0, 1, 2, 3, 4, 5)
observed_data <- matrix(c(
  0, 1, 4, 10, 20, 40,
  0, 2, 6, 12, 24, 50
), ncol = 2, byrow = TRUE)
initial_conditions <- c(0.1, 0, 0)  # Initial G_on, m, p

data <- list(
  N = length(time_points),
  t0 = 0.0,
  ts = time_points,
  y0 = initial_conditions,
  y_obs = observed_data
)

# Compile and fit the model
stan_model <- stan_model(file = 'bursting_gene_model.stan')
fit <- sampling(stan_model, data = data, iter = 2000, warmup = 1000, chains = 4)

# Print the results
print(fit)