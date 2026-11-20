# Reverse Kaplan-Meier distribution and identifiable quantiles ----------------

apply_reverse_km_method <- function(observed, replicate_id = NA_integer_,
                                    probs = c(0.025, 0.10, 0.50, 0.90, 0.975)) {
  ensure_censoring_packages_attached()
  validate_observed_data(observed)
  if (any(!is.finite(probs)) || any(probs < 0 | probs > 1)) stop("probs must be in [0, 1]")
  captured <- capture_conditions(NADA::cenfit(observed$value, observed$censored))
  fit <- captured$value
  if (inherits(fit, "method_error"))
    return(new_method_result(replicate_id, "reverse_km", observed, "failed",
      warnings = captured$warnings, error = fit$message))

  # A common left-censoring threshold leaves quantiles below the censored
  # proportion unidentified. Keep them NA instead of using package moments.
  p_cen <- mean(observed$censored)
  identifiable <- probs > p_cen
  q <- rep(NA_real_, length(probs)); names(q) <- paste0("q", probs)
  quantile_warnings <- character()
  if (any(identifiable)) {
    qresult <- capture_conditions(quantile(fit, probs = probs[identifiable]))
    quantile_warnings <- qresult$warnings
    qfit <- qresult$value
    if (!inherits(qfit, "method_error")) q[identifiable] <- as.numeric(qfit)
  }
  estimates <- list(quantiles = q, mean = NA_real_, median = if (0.5 > p_cen) unname(q[which.min(abs(probs - 0.5))]) else NA_real_,
                    sd = NA_real_, cv = NA_real_)
  status <- if (!any(identifiable)) "non_identifiable" else if (all(identifiable)) "success" else "partial"
  new_method_result(replicate_id, "reverse_km", observed, status,
    inference_status = "not_implemented", warnings = unique(c(captured$warnings, quantile_warnings)),
    estimates = estimates,
    details = list(fit = fit, identifiable_quantiles = setNames(identifiable, names(q)),
      non_identifiable_reason = "Quantiles at or below the common censoring proportion are not identified."))
}
