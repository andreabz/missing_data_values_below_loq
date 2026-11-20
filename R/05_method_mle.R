# Left-censored lognormal maximum likelihood ---------------------------------

apply_mle_method <- function(observed, replicate_id = NA_integer_) {
  ensure_censoring_packages_attached()
  validate_observed_data(observed)
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
  converged <- is.null(sr$fail) || identical(sr$fail, 0L)
  finite <- length(mu) == 1L && is.finite(mu) && length(sigma) == 1L && is.finite(sigma) && sigma > 0
  status <- if (!converged) "non_converged" else if (!finite) "partial" else "success"
  estimates <- if (finite) list(mean = exp(mu + sigma^2 / 2), median = exp(mu),
      sd = sqrt((exp(sigma^2) - 1) * exp(2 * mu + sigma^2)),
      cv = sqrt(exp(sigma^2) - 1)) else
    list(mean = NA_real_, median = NA_real_, sd = NA_real_, cv = NA_real_)
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
