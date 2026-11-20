# Independent, small-scale statistical validation for censored estimators.
# This module is deliberately separate from the production Monte Carlo path.

censored_lognormal_loglik <- function(meanlog, log_sdlog, observed) {
  validate_observed_data(observed)
  if (length(meanlog) != 1L || !is.finite(meanlog) ||
      length(log_sdlog) != 1L || !is.finite(log_sdlog))
    stop("meanlog and log_sdlog must be finite scalars")
  if (any(observed$loq <= 0) || any(observed$value[!observed$censored] <= 0))
    stop("The censored lognormal likelihood requires positive LOQs and detections")
  sdlog <- exp(log_sdlog)
  if (!is.finite(sdlog) || sdlog <= 0) return(-Inf)
  detected <- !observed$censored
  sum(stats::dlnorm(observed$value[detected], meanlog, sdlog, log = TRUE)) +
    sum(stats::plnorm(observed$loq[observed$censored], meanlog, sdlog, log.p = TRUE))
}

fit_independent_censored_lognormal <- function(observed) {
  validate_observed_data(observed)
  detected <- observed$value[!observed$censored]
  if (any(observed$loq <= 0) || any(detected <= 0))
    return(list(status = "not_applicable", reason = "Positive LOQ and detections are required."))
  if (length(detected) < 2L || length(unique(detected)) < 2L)
    return(list(status = "non_identifiable", reason =
      "At least two distinct detections are required; all-censored samples are not identifiable."))

  start <- c(mean(log(detected)), log(max(stats::sd(log(detected)), 0.15)))
  objective <- function(par) {
    if (par[2L] < -8 || par[2L] > 5) return(1e100)
    value <- censored_lognormal_loglik(par[1L], par[2L], observed)
    if (!is.finite(value)) 1e100 else -value
  }
  fit <- stats::optim(start, objective, method = "BFGS",
    control = list(maxit = 2000, reltol = 1e-11), hessian = TRUE)
  boundary <- fit$par[2L] <= -7.999 || fit$par[2L] >= 4.999
  ok <- fit$convergence == 0L && is.finite(fit$value) && !boundary
  params <- c(meanlog = fit$par[1L], sdlog = exp(fit$par[2L]))
  list(status = if (ok) "success" else "non_converged", parameters = params,
    loglik = -fit$value, convergence = fit$convergence, message = fit$message,
    iterations = fit$counts[["function"]], boundary = boundary,
    moments = if (ok) lognormal_moments(params[["meanlog"]], params[["sdlog"]]) else NULL)
}

validation_summary <- function(estimates, truth) {
  estimates <- estimates[is.finite(estimates)]
  if (!length(estimates)) return(c(bias = NA_real_, rmse = NA_real_, empirical_sd = NA_real_))
  error <- estimates - truth
  c(bias = mean(error), rmse = sqrt(mean(error^2)),
    empirical_sd = if (length(estimates) > 1L) stats::sd(estimates) else NA_real_)
}
