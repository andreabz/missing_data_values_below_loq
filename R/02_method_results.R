# Shared result contract and condition capture --------------------------------

new_method_result <- function(
  replicate_id, method, observed, estimate_status = "not_applicable",
  inference_status = "not_implemented", warnings = character(), error = NA_character_,
  estimates = list(), parameters = list(), interval = c(lower = NA_real_, upper = NA_real_),
  ucl95 = NA_real_, details = list()
) {
  allowed <- c("success", "partial", "non_identifiable", "non_converged",
               "not_applicable", "failed", "not_implemented")
  stopifnot(estimate_status %in% allowed, inference_status %in% allowed)
  n <- nrow(observed)
  n_censored <- sum(observed$censored)
  interval <- setNames(as.numeric(interval), c("lower", "upper"))
  width <- if (all(is.finite(interval))) diff(interval) else NA_real_
  list(
    replicate_id = replicate_id, n = n, loq = unique(observed$loq), method = method,
    n_censored = n_censored, censoring_proportion = n_censored / n,
    estimates = estimates, parameters = parameters,
    inference = list(status = inference_status, ci_mean_95 = interval,
                     ci_width = width, ucl95 = ucl95),
    estimate_status = estimate_status, warnings = warnings,
    error = error, details = details
  )
}
capture_conditions <- function(expr) {
  warnings <- character()
  value <- tryCatch(
    withCallingHandlers(expr, warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }),
    error = function(e) structure(list(message = conditionMessage(e)), class = "method_error")
  )
  list(value = value, warnings = unique(warnings))
}

ensure_censoring_packages_attached <- function() {
  # These package versions resolve a few legacy helper names from the search
  # path even when called with ::; attach only the required method packages.
  if (!"package:NADA" %in% search()) suppressPackageStartupMessages(library(package = "NADA"))
  if (!"package:NADA2" %in% search()) suppressPackageStartupMessages(library(package = "NADA2"))
  invisible(TRUE)
}

validate_observed_data <- function(observed) {
  required <- c("value", "censored", "loq")
  if (!is.data.frame(observed) || !all(required %in% names(observed)))
    stop("observed must be a data.frame with value, censored, and loq columns")
  if (!is.numeric(observed$value) || !is.logical(observed$censored) ||
      !is.numeric(observed$loq) || anyNA(observed$censored) ||
      any(!is.finite(observed$value)) || any(!is.finite(observed$loq)) ||
      length(unique(observed$loq)) != 1L || nrow(observed) < 1L ||
      length(observed$value) != nrow(observed) || length(observed$censored) != nrow(observed) ||
      length(observed$loq) != nrow(observed))
    stop("observed data have invalid values, indicators, or non-common LOQ")
  if (any(observed$value[observed$censored] != observed$loq[observed$censored]))
    stop("censored observations must be recorded at their LOQ")
  if (any(observed$value[!observed$censored] <= observed$loq[!observed$censored]))
    stop("uncensored observations must be strictly greater than their LOQ")
  invisible(TRUE)
}
