# Simulation configuration -----------------------------------------------

mu <- 0
sigma <- 1

sample_sizes <- c(30L, 100L, 500L)

loq_probabilities <- c(
  LOQ_10 = 0.10,
  LOQ_40 = 0.40
)

loq_values <- stats::qlnorm(
  loq_probabilities,
  meanlog = mu,
  sdlog = sigma
)

n_monte_carlo <- 2000L

population_parameters <- list(
  mean = exp(mu + sigma^2 / 2),
  median = exp(mu),
  sd = sqrt((exp(sigma^2) - 1) * exp(2 * mu + sigma^2)),
  cv = sqrt(exp(sigma^2) - 1)
)

simulation_config <- list(
  distribution = "lognormal",
  meanlog = mu,
  sdlog = sigma,
  sample_sizes = sample_sizes,
  loq_probabilities = loq_probabilities,
  loq_values = loq_values,
  n_monte_carlo = n_monte_carlo,
  population_parameters = population_parameters
)

rm(mu, sigma)
