library(rstan)

# Example data
time_points <- c(0.1, 10, 20, 30, 40, 50)
observed_data <- matrix(c(
  1, 0, 0,
  0, 1, 20,
  1, 0, 40,
  0, 1, 50,
  1, 0, 20,
  1, 0, 15
), ncol = 3, byrow = TRUE)
initial_conditions <- c(1, 0, 0)  # Initial OFF, ON, Protein

data <- list(
  T = length(time_points),
  ts = time_points,
  y0 = initial_conditions,
  y = observed_data
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
stan_model <- stan_model(file = 'bursting_gene_model.stan')
fit <- sampling(stan_model, data = data, init = init_fun, iter = 2000, warmup = 1000, chains = 4, verbose = TRUE)

# Print the results
print(fit)