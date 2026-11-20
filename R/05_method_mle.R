# Left-censored lognormal maximum likelihood ---------------------------------

lognormal_moments <- function(meanlog, sdlog) {
  if (length(meanlog) != 1L || !is.finite(meanlog) || length(sdlog) != 1L ||
      !is.finite(sdlog) || sdlog <= 0) stop("meanlog must be finite and sdlog must be finite and positive")
  mean <- exp(meanlog + sdlog^2 / 2)
  median <- exp(meanlog)
  sd <- sqrt((exp(sdlog^2) - 1) * exp(2 * meanlog + sdlog^2))
  cv <- sqrt(exp(sdlog^2) - 1)
  list(mean = mean, median = median, sd = sd, cv = cv)
}

apply_mle_method <- function(observed, replicate_id = NA_integer_) {
  ensure_censoring_packages_attached()
  validate_observed_data(observed)
  if (any(observed$loq <= 0) || any(observed$value[!observed$censored] <= 0))
    return(new_method_result(replicate_id, "mle_lognormal", observed, "not_applicable",
      error = "Lognormal MLE requires a positive LOQ and positive detected observations."))
  detected <- observed$value[!observed$censored]
  if (length(detected) < 2L || length(unique(detected)) < 2L) {
    msg <- if (!length(detected)) "All observations are censored; lognormal parameters are not identifiable." else
      "At least two distinct detected values are required for a stable lognormal fit."
    return(new_method_result(replicate_id, "mle_lognormal", observed,
      estimate_status = if (!length(detected)) "non_identifiable" else "not_applicable",
      inference_status = "not_implemented", error = msg))
  }
  if (any(detected <= 0))
    return(new_method_result(replicate_id, "mle_lognormal", observed, "not_applicable",
      error = "Lognormal MLE requires positive detected observations."))

  captured <- capture_conditions(NADA::cenmle(
    obs = observed$value, censored = observed$censored,
    dist = "lognormal", conf.int = 0.95
  ))
  fit <- captured$value
  if (inherits(fit, "method_error"))
    return(new_method_result(replicate_id, "mle_lognormal", observed, "failed",
      warnings = captured$warnings, error = fit$message))

  sr <- fit@survreg
  mu <- unname(stats::coef(sr)[1L])
  sigma <- unname(sr$scale)
  fail <- sr$fail
  converged <- is.null(fail) || (length(fail) == 1L && is.numeric(fail) && is.finite(fail) && fail == 0)
  if (any(grepl("did not converge|ran out of iterations", captured$warnings, ignore.case = TRUE)))
    converged <- FALSE
  finite_params <- length(mu) == 1L && is.finite(mu) && length(sigma) == 1L && is.finite(sigma) && sigma > 0
  estimates <- if (finite_params) tryCatch(lognormal_moments(mu, sigma), error = function(e) NULL) else NULL
  finite_estimates <- !is.null(estimates) && all(is.finite(unlist(estimates)))
  status <- if (!converged) "non_converged" else if (!finite_params || !finite_estimates) "failed" else "success"
  if (!finite_estimates) estimates <- list(mean = NA_real_, median = NA_real_, sd = NA_real_, cv = NA_real_)
  vc <- tryCatch(stats::vcov(sr), error = function(e) NULL)
  se <- if (!is.null(vc)) sqrt(diag(vc)) else numeric()
  parameters <- list(meanlog = mu, sdlog = sigma,
    se_meanlog = if (length(se) >= 1L) unname(se[1L]) else NA_real_,
    se_log_sdlog = if (length(se) >= 2L) unname(se[2L]) else NA_real_)
  new_method_result(replicate_id, "mle_lognormal", observed, status,
    inference_status = "not_implemented", warnings = captured$warnings,
    error = if (status == "non_converged") "survreg did not converge" else NA_character_,
    estimates = estimates, parameters = parameters,
    details = list(fit = fit, converged = converged,
      inference_note = "Mean interval not implemented pending methodological validation"))
}
