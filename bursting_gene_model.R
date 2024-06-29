library(rstan)
library(ggplot2)

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
fit <- sampling(stan_model, data = data, init = init_fun, iter = 400, warmup = 200, chains = 4, verbose = TRUE)

# Print the results
print(fit)

# Plot traces to check convergence and save the plot of the chains to a file
traceplot(fit)
ggsave("traceplot.png", traceplot(fit))

# Display summary statistics
summary(fit)$summary
