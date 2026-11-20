# Reverse Kaplan-Meier distribution and identifiable quantiles ----------------

apply_reverse_km_method <- function(observed, replicate_id = NA_integer_,
                                    probs = c(0.025, 0.10, 0.50, 0.90, 0.975)) {
  ensure_censoring_packages_attached()
  validate_observed_data(observed)
  if (!is.numeric(probs) || !length(probs) || any(!is.finite(probs)) ||
      any(probs < 0 | probs > 1) || anyDuplicated(probs))
    stop("probs must be a non-empty numeric vector of unique finite values in [0, 1]")
  captured <- capture_conditions(NADA::cenfit(observed$value, observed$censored))
  fit <- captured$value
  if (inherits(fit, "method_error"))
    return(new_method_result(replicate_id, "reverse_km", observed, "failed",
      warnings = captured$warnings, error = fit$message))

  # NADA's reverse-KM quantile method can return NA even above the censoring
  # fraction (for example, beyond the estimable tail). Query each requested
  # probability separately and retain only values actually returned finite.
  p_cen <- mean(observed$censored)
  q <- rep(NA_real_, length(probs)); names(q) <- paste0("q", probs)
  quantile_warnings <- character()
  quantile_errors <- character()
  for (i in seq_along(probs)) {
    if (probs[i] <= p_cen) next
    qresult <- capture_conditions(quantile(fit, probs = probs[i]))
    quantile_warnings <- c(quantile_warnings, qresult$warnings)
    if (inherits(qresult$value, "method_error")) {
      quantile_errors <- c(quantile_errors, qresult$value$message)
      next
    }
    value <- suppressWarnings(as.numeric(qresult$value)[1L])
    if (length(value) == 1L && is.finite(value)) q[i] <- value
  }
  identifiable <- is.finite(q) & probs > p_cen
  median_index <- which(probs == 0.5)
  median_value <- if (length(median_index) == 1L && identifiable[median_index]) q[median_index] else NA_real_
  estimates <- list(quantiles = q, mean = NA_real_, median = median_value,
                    sd = NA_real_, cv = NA_real_)
  status <- if (length(quantile_errors) && !any(identifiable)) "failed" else
    if (all(identifiable)) "success" else if (any(identifiable)) "partial" else "non_identifiable"
  new_method_result(replicate_id, "reverse_km", observed, status,
    inference_status = "not_implemented", warnings = unique(c(captured$warnings, quantile_warnings)),
    error = if (length(quantile_errors)) paste(unique(quantile_errors), collapse = " | ") else NA_character_,
    estimates = estimates,
    details = list(fit = fit, identifiable_quantiles = setNames(identifiable, names(q)),
      non_identifiable_reason = "Quantiles at or below the common censoring proportion, and quantiles NADA does not return finitely, are unavailable.",
      quantile_errors = unique(quantile_errors)))
}
