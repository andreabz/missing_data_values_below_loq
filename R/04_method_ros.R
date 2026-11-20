# Log-ROS point estimation. Inference is intentionally deferred. --------------

apply_ros_method <- function(observed, replicate_id = NA_integer_) {
  ensure_censoring_packages_attached()
  validate_observed_data(observed)
  if (any(observed$loq <= 0) || any(observed$value[!observed$censored] <= 0))
    return(new_method_result(replicate_id, "ros_log", observed, "not_applicable",
      error = "Log-ROS requires a positive LOQ and positive detected observations."))
  n_detected <- sum(!observed$censored)
  if (n_detected < 2L || length(unique(observed$value[!observed$censored])) < 2L) {
    return(new_method_result(replicate_id, "ros_log", observed,
      estimate_status = if (n_detected == 0L) "non_identifiable" else "not_applicable",
      inference_status = "not_implemented",
      error = if (n_detected == 0L) "No uncensored observations identify the ROS regression." else
        "ROS requires at least two distinct uncensored values."))
  }
  if (any(observed$value[!observed$censored] <= 0)) {
    return(new_method_result(replicate_id, "ros_log", observed, "not_applicable",
      error = "Log-ROS requires positive uncensored observations."))
  }
  captured <- capture_conditions(NADA2::ros(
    obs = observed$value, censored = observed$censored,
    forwardT = "log", reverseT = "exp"
  ))
  fit <- captured$value
  if (inherits(fit, "method_error"))
    return(new_method_result(replicate_id, "ros_log", observed, "failed",
      warnings = captured$warnings, error = fit$message))

  x <- fit$modeled
  coefficients <- tryCatch(stats::coef(fit$model), error = function(e) numeric())
  valid_fit <- is.numeric(x) && length(x) == nrow(observed) &&
    all(is.finite(x)) && all(x > 0) && length(coefficients) >= 2L &&
    all(is.finite(coefficients))
  if (!valid_fit) {
    return(new_method_result(replicate_id, "ros_log", observed, "failed",
      warnings = captured$warnings,
      error = "ROS returned non-finite, non-positive, or incomplete fitted values/parameters.",
      details = list(fit = fit)))
  }
  estimates <- list(mean = mean(x), median = stats::median(x),
    sd = if (length(x) > 1L) stats::sd(x) else NA_real_,
    cv = if (mean(x) > 0) stats::sd(x) / mean(x) else NA_real_)
  new_method_result(replicate_id, "ros_log", observed,
    estimate_status = if (all(is.finite(unlist(estimates)))) "success" else "failed",
    inference_status = "not_implemented", warnings = captured$warnings,
    estimates = estimates, parameters = as.list(coefficients),
    details = list(fit = fit, bootstrap = "not_implemented; future bootstrap must refit ROS each draw"))
}

bootstrap_ros <- function(...) {
  stop("ROS bootstrap is not implemented; when added it must refit ROS in every resample.")
}
