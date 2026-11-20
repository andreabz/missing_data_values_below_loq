# Basic checks for the sample generator -----------------------------------

source("R/01_generate_sample.R")

set.seed(20261120)

x <- generate_left_censored_sample(
  n = 100L,
  loq = stats::qlnorm(0.40, meanlog = 0, sdlog = 1),
  meanlog = 0,
  sdlog = 1
)

stopifnot(
  nrow(x$observed) == 100L,
  length(x$truth) == 100L,
  all(x$observed$censored == (x$truth < x$metadata$loq)),
  all(x$observed$value[!x$observed$censored] == x$truth[!x$observed$censored]),
  all(x$observed$value[x$observed$censored] == x$metadata$loq),
  x$metadata$n_censored == sum(x$observed$censored),
  isTRUE(all.equal(
    x$metadata$theoretical_censoring_probability,
    0.40,
    tolerance = 1e-12
  ))
)

cat("Generator checks passed.\n")
print(x$metadata)
