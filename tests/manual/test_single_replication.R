# Functional checks for one-replicate architecture; no Monte Carlo loop.
for (file in c("R/02_method_results.R", "R/03_methods_simple.R", "R/04_method_ros.R",
               "R/05_method_mle.R", "R/06_method_reverse_km.R", "R/07_inference.R",
               "R/08_single_replication.R", "R/09_evaluation.R")) source(file)

make_sample <- function(censored, detected, loq = 0.5) {
  n <- length(censored)
  values <- rep(loq, n)
  values[!censored] <- detected
  list(observed = data.frame(id = seq_len(n), value = values,
    censored = censored, loq = rep(loq, n)), truth = seq_len(n) / n,
    metadata = list(n = n, loq = loq, n_censored = sum(censored)))
}

scenarios <- list(
  none = make_sample(rep(FALSE, 8), seq(0.7, 2.1, length.out = 8)),
  partial = make_sample(c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
                        seq(0.7, 2.1, length.out = 6)),
  elevated = make_sample(c(rep(TRUE, 6), FALSE, FALSE), c(1.1, 1.8)),
  all = make_sample(rep(TRUE, 8), numeric()),
  no_positive_detects = make_sample(c(TRUE, TRUE, FALSE, FALSE), c(0, 0)),
  nonfinite = make_sample(c(FALSE, FALSE, FALSE), c(1, Inf, 2))
)

checks <- list()
for (nm in names(scenarios)) {
  s <- scenarios[[nm]]
  if (nm == "nonfinite") {
    stopifnot(inherits(try(validate_observed_data(s$observed), silent = TRUE), "try-error"))
    checks[[nm]] <- "input rejected as expected"
    next
  }
  r <- run_single_replication(s, replicate_id = match(nm, names(scenarios)))
  stopifnot(identical(names(r$methods), c("complete_case", "substitute_zero",
    "substitute_loq_half", "substitute_loq", "ros_log", "mle_lognormal", "reverse_km")))
  stopifnot(all(vapply(r$methods, function(z) z$n == nrow(s$observed), logical(1))))
  stopifnot(all(vapply(r$methods, function(z) z$n_censored == sum(s$observed$censored), logical(1))))
  stopifnot(!any(vapply(r$methods, function(z) "truth" %in% names(z), logical(1))))
  if (nm == "all") {
    stopifnot(r$methods$mle_lognormal$estimate_status == "non_identifiable")
    stopifnot(r$methods$ros_log$estimate_status == "non_identifiable")
  }
  if (nm == "no_positive_detects") {
    stopifnot(r$methods$mle_lognormal$estimate_status == "not_applicable")
    stopifnot(r$methods$ros_log$estimate_status == "not_applicable")
  }
  ev <- evaluate_replication(r, s$truth)
  stopifnot(is.finite(ev$sample_truth$mean))
  checks[[nm]] <- vapply(r$methods, `[[`, character(1), "estimate_status")
}

# The public wrapper identifies the all-censored case before fitting. Confirm
# separately that the installed backend itself warns/fails to converge there.
ensure_censoring_packages_attached()
raw_failure <- capture_conditions(NADA::cenmle(rep(9.5, 8), rep(TRUE, 8), dist = "lognormal"))
stopifnot(!inherits(raw_failure$value, "method_error"))
raw_fit <- raw_failure$value@survreg
stopifnot(length(raw_failure$warnings) > 0L || !is.null(raw_fit$fail))

sample <- scenarios$partial
single <- run_single_replication(sample, replicate_id = 99L)
stopifnot(all(vapply(single$methods, function(x) is.na(x$inference$ci_width), logical(1))))
stopifnot(is.na(single$methods$reverse_km$estimates$mean))
stopifnot(is.na(single$methods$mle_lognormal$inference$ucl95))
stopifnot(!"truth" %in% names(single))

cat("Single replication functional checks passed.\n")
print(checks)
