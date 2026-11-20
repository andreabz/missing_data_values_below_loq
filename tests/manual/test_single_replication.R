# Functional checks for one-replicate architecture; no Monte Carlo loop.
for (file in c("R/02_method_results.R", "R/03_methods_simple.R", "R/04_method_ros.R",
               "R/05_method_mle.R", "R/06_method_reverse_km.R", "R/07_inference.R",
               "R/08_single_replication.R", "R/09_evaluation.R")) source(file)

make_sample <- function(censored, detected, loq = 0.5) {
  n <- length(censored)
  values <- rep(loq, n)
  values[!censored] <- detected
  truth <- values
  n_censored <- sum(censored)
  if (n_censored) truth[censored] <- loq * seq(0.1, 0.9, length.out = n_censored)
  list(observed = data.frame(id = seq_len(n), value = values,
    censored = censored, loq = rep(loq, n)), truth = truth,
    metadata = list(n = n, loq = loq, n_censored = sum(censored)))
}

scenarios <- list(
  none = make_sample(rep(FALSE, 8), seq(0.7, 2.1, length.out = 8)),
  partial = make_sample(c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
                        seq(0.7, 2.1, length.out = 6)),
  elevated = make_sample(c(rep(TRUE, 6), FALSE, FALSE), c(1.1, 1.8)),
  all = make_sample(rep(TRUE, 8), numeric()),
  nonfinite = make_sample(c(FALSE, FALSE, FALSE), c(1, Inf, 2))
)

robustness_cases <- list(
  nonpositive_loq = make_sample(rep(FALSE, 3), c(0.5, 1, 1.5), loq = 0),
  malformed_censored_value = make_sample(c(TRUE, FALSE), 1),
  malformed_detected_value = make_sample(c(FALSE, FALSE), c(0.5, 1))
)
robustness_cases$malformed_censored_value$observed$value[1] <- 0.4

expected_observed <- data.frame(value = c(0.5, 0.5, 0.8, 1.2, 1.6),
  censored = c(TRUE, TRUE, FALSE, FALSE, FALSE), loq = rep(0.5, 5))
expected_stats <- list(
  complete_case = c(mean = 1.2, median = 1.2, sd = 0.4),
  substitute_zero = c(mean = 0.72, median = 0.8, sd = sqrt(0.512)),
  substitute_loq_half = c(mean = 0.82, median = 0.8, sd = sqrt(0.35075)),
  substitute_loq = c(mean = 0.92, median = 0.8, sd = sqrt(0.227))
)
for (method in names(expected_stats)) {
  result <- apply_simple_method(expected_observed, method, replicate_id = 7L)
  stopifnot(isTRUE(all.equal(unlist(result$estimates[c("mean", "median", "sd")]),
    expected_stats[[method]], tolerance = 1e-12)))
}
known_lognormal <- lognormal_moments(0, 1)
stopifnot(isTRUE(all.equal(unlist(known_lognormal),
  c(mean = exp(0.5), median = 1, sd = sqrt((exp(1) - 1) * exp(1)),
    cv = sqrt(exp(1) - 1)), tolerance = 1e-12)))
validate_observed_data(expected_observed)
validate_observed_data(robustness_cases$nonpositive_loq$observed)

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
  stopifnot(all(s$observed$value[s$observed$censored] == s$observed$loq[s$observed$censored]))
  stopifnot(all(s$truth[s$observed$censored] < s$observed$loq[s$observed$censored]))
  stopifnot(all(s$truth[!s$observed$censored] == s$observed$value[!s$observed$censored]))
  stopifnot(!any(vapply(r$methods, function(z) "truth" %in% names(z), logical(1))))
  if (nm == "all") {
    stopifnot(r$methods$mle_lognormal$estimate_status == "non_identifiable")
    stopifnot(r$methods$ros_log$estimate_status == "non_identifiable")
  }
  ev <- evaluate_replication(r, s$truth)
  stopifnot(is.finite(ev$sample_truth$mean))
  checks[[nm]] <- vapply(r$methods, `[[`, character(1), "estimate_status")
}

for (nm in names(robustness_cases)) {
  sample <- robustness_cases[[nm]]
  if (nm == "nonpositive_loq") {
    z <- run_single_replication(sample)
    stopifnot(z$methods$mle_lognormal$estimate_status == "not_applicable")
    stopifnot(z$methods$ros_log$estimate_status == "not_applicable")
  } else {
    stopifnot(inherits(try(validate_observed_data(sample$observed), silent = TRUE), "try-error"))
  }
}

# NADA returns only actually estimable reverse-KM quantiles. For this sample
# q(.5) is finite while q(.025) is below the censoring fraction and q(.975) is
# beyond the finite tail of the fitted reverse-KM curve.
reverse_sample <- make_sample(c(rep(TRUE, 16), rep(FALSE, 24)),
  exp(seq(log(0.8), log(9), length.out = 24)))
reverse_result <- apply_reverse_km_method(as.data.frame(reverse_sample$observed))
stopifnot(is.finite(reverse_result$estimates$quantiles[["q0.5"]]))
stopifnot(is.na(reverse_result$estimates$quantiles[["q0.025"]]))
stopifnot(is.na(reverse_result$estimates$quantiles[["q0.975"]]))
stopifnot(is.finite(reverse_result$estimates$median))
median_not_requested <- apply_reverse_km_method(as.data.frame(reverse_sample$observed), probs = 0.6)
stopifnot(is.na(median_not_requested$estimates$median))

base_sample <- scenarios$partial
rep_a <- run_single_replication(base_sample, replicate_id = 41L, methods = "substitute_loq_half")
rep_b <- run_single_replication(base_sample, replicate_id = 42L, methods = "substitute_loq")
collected <- collect_replication_results(list(rep_a, rep_b))
stopifnot(inherits(collected, "replication_result_collection"))
stopifnot(identical(names(collected), c("41", "42")))
stopifnot(identical(collected[[1]]$substitute_loq_half$replicate_id, 41L))
stopifnot(identical(collected[[1]]$substitute_loq_half$method, "substitute_loq_half"))
stopifnot(identical(collected[[2]]$substitute_loq$replicate_id, 42L))

eval_sample <- make_sample(c(TRUE, TRUE, FALSE, FALSE), c(1, 1.5))
eval_rep <- run_single_replication(eval_sample, replicate_id = 55L, methods = "substitute_loq")
evaluation <- evaluate_replication(eval_rep, eval_sample$truth,
  population_parameters = list(mean = 1.2, median = 1, sd = 0.5, cv = 0.4))
row <- evaluation$methods[[1]]
stopifnot(is.finite(row$errors_vs_sample_truth$mean))
stopifnot(isTRUE(all.equal(row$errors_vs_population_parameters$mean,
  row$estimates$mean - 1.2)))
missing_reference <- evaluate_replication(eval_rep, eval_sample$truth)
stopifnot(is.na(missing_reference$methods[[1]]$errors_vs_population_parameters$mean))

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


# Input-contract edge cases introduced in the review.
stopifnot(inherits(try(apply_reverse_km_method(reverse_sample$observed, probs = numeric()),
  silent = TRUE), "try-error"))
stopifnot(inherits(try(apply_reverse_km_method(reverse_sample$observed, probs = "0.5"),
  silent = TRUE), "try-error"))
stopifnot(inherits(try(apply_reverse_km_method(reverse_sample$observed, probs = c(0.5, 0.5)),
  silent = TRUE), "try-error"))
stopifnot(inherits(try(evaluate_replication(eval_rep, replace(eval_sample$truth, 1, Inf)),
  silent = TRUE), "try-error"))
duplicate_id_error <- try(collect_replication_results(list(rep_a, rep_a)), silent = TRUE)
stopifnot(inherits(duplicate_id_error, "try-error"))

# Deterministic API cross-checks: wrappers must reproduce the corresponding
# installed NADA/NADA2 point estimates for the same fixed observed sample.
n_check <- 40L
truth_check <- stats::qlnorm((seq_len(n_check) - 0.5) / n_check, meanlog = 0, sdlog = 1)
loq_check <- stats::qlnorm(0.25, meanlog = 0, sdlog = 1)
censored_check <- truth_check < loq_check
observed_check <- data.frame(
  value = ifelse(censored_check, loq_check, truth_check),
  censored = censored_check, loq = rep(loq_check, n_check)
)
ros_direct <- NADA2::ros(observed_check$value, observed_check$censored,
  forwardT = "log", reverseT = "exp")
ros_wrapped <- apply_ros_method(observed_check)
stopifnot(ros_wrapped$estimate_status == "success")
stopifnot(isTRUE(all.equal(ros_wrapped$estimates$mean, mean(ros_direct$modeled),
  tolerance = 1e-10)))
stopifnot(isTRUE(all.equal(ros_wrapped$estimates$median, stats::median(ros_direct$modeled),
  tolerance = 1e-10)))

mle_direct <- NADA::cenmle(observed_check$value, observed_check$censored,
  dist = "lognormal", conf.int = 0.95)
mle_direct_moments <- lognormal_moments(unname(stats::coef(mle_direct@survreg)[1L]),
  unname(mle_direct@survreg$scale))
mle_wrapped <- apply_mle_method(observed_check)
stopifnot(mle_wrapped$estimate_status == "success")
stopifnot(isTRUE(all.equal(mle_wrapped$estimates$mean, mle_direct_moments$mean,
  tolerance = 1e-10)))
stopifnot(isTRUE(all.equal(mle_wrapped$estimates$sd, mle_direct_moments$sd,
  tolerance = 1e-10)))

# Reverse-KM behavior across several fixed censoring thresholds.
for (censor_probability in c(0, 0.10, 0.25, 0.40, 0.75)) {
  threshold <- stats::qlnorm(censor_probability, meanlog = 0, sdlog = 1)
  censored_level <- truth_check < threshold
  if (censor_probability == 0) censored_level[] <- FALSE
  level_data <- data.frame(
    value = ifelse(censored_level, threshold, truth_check),
    censored = censored_level, loq = rep(threshold, n_check)
  )
  level_result <- apply_reverse_km_method(level_data)
  stopifnot(length(level_result$estimates$quantiles) == 5L)
  stopifnot(all(is.na(level_result$estimates$quantiles[
    c(0.025, 0.10, 0.50, 0.90, 0.975) <= mean(censored_level)])))
}

stopifnot(mle_wrapped$parameters$vcov_status %in% c("valid", "unavailable_or_invalid"))
bad_numeric_id <- run_single_replication(base_sample, replicate_id = Inf,
  methods = "substitute_loq")
stopifnot(inherits(try(collect_replication_results(list(bad_numeric_id)), silent = TRUE),
  "try-error"))
bad_character_id <- run_single_replication(base_sample, replicate_id = "",
  methods = "substitute_loq")
stopifnot(inherits(try(collect_replication_results(list(bad_character_id)), silent = TRUE),
  "try-error"))
