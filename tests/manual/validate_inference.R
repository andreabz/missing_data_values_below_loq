# Deterministic checks for method-specific inferential procedures.
for (file in c("R/02_method_results.R", "R/04_method_ros.R", "R/05_method_mle.R",
               "R/06_method_reverse_km.R", "R/07_inference.R", "R/08_single_replication.R",
               "R/11_inference_validation_decisions.R"))
  source(file)
ensure_censoring_packages_attached()

stopifnot(inherits(try(validate_confidence_level(1), silent = TRUE), "try-error"))
stopifnot(inherits(try(validate_bootstrap_replicates(1), silent = TRUE), "try-error"))
stopifnot(identical(validate_bootstrap_replicates(20), 20L))

# Exact coverage decisions separate significant undercoverage, equivalence,
# inconclusive precision, and interval availability.
undercovered <- classify_inference_coverage(430, 500, availability = 1)
stopifnot(undercovered$status == "not validated", undercovered$undercoverage_p < 0.05)
well_calibrated <- classify_inference_coverage(475, 500, availability = 1)
stopifnot(well_calibrated$status == "validated")
well_calibrated_ci <- stats::binom.test(475, 500, conf.level = 0.90)$conf.int
stopifnot(isTRUE(all.equal(well_calibrated$lower, unname(well_calibrated_ci[[1L]]))),
  isTRUE(all.equal(well_calibrated$upper, unname(well_calibrated_ci[[2L]]))),
  well_calibrated$lower >= 0.93, well_calibrated$upper <= 0.97)
expected_undercoverage_p <- stats::binom.test(430, 500, p = 0.95,
  alternative = "less")$p.value
stopifnot(isTRUE(all.equal(undercovered$undercoverage_p, expected_undercoverage_p)))
weak_undercoverage <- classify_inference_coverage(468, 500, availability = 1)
stopifnot(weak_undercoverage$undercoverage_p >= 0.05,
  weak_undercoverage$status == "inconclusive", weak_undercoverage$lower < 0.93)
imprecise <- classify_inference_coverage(9, 10, availability = 1)
stopifnot(imprecise$status == "inconclusive")
unavailable <- classify_inference_coverage(475, 500, availability = 0.90)
stopifnot(unavailable$status == "inconclusive")
unconfirmed_stability <- classify_inference_coverage(475, 500, availability = 1,
  stability_confirmed = FALSE)
stopifnot(unconfirmed_stability$status == "inconclusive")
structural <- try(classify_inference_coverage(NA_real_, 0, availability = NA_real_), silent = TRUE)
stopifnot(inherits(structural, "try-error"))
structural_rkm <- classify_reverse_km_coverage(NA_real_, 0, NA_real_, applicable = FALSE)
stopifnot(structural_rkm$status == "not applicable", is.na(structural_rkm$coverage))
no_finite_rkm <- classify_reverse_km_coverage(0, 0, 0, applicable = TRUE)
stopifnot(no_finite_rkm$status == "inconclusive", is.na(no_finite_rkm$coverage),
  no_finite_rkm$availability == 0)

invalid_decisions <- list(
  function() classify_inference_coverage(1.5, 10),
  function() classify_inference_coverage(-1, 10),
  function() classify_inference_coverage(11, 10),
  function() classify_inference_coverage(1, 0),
  function() classify_inference_coverage(1, 10, availability = NA_real_),
  function() classify_inference_coverage(1, 10, availability = 1.1),
  function() classify_inference_coverage(1, 10, nominal = 1),
  function() classify_inference_coverage(1, 10, nominal = NA_real_),
  function() classify_inference_coverage(1, 10, tolerance = 0),
  function() classify_inference_coverage(1, 10, alpha = 0.5),
  function() classify_inference_coverage(1, 10, min_availability = -0.1),
  function() classify_inference_coverage(1, 10, stability_confirmed = NA))
stopifnot(all(vapply(invalid_decisions, function(f) inherits(try(f(), silent = TRUE),
  "try-error"), logical(1))))

ci_counts <- coverage_metrics(c(0, 2, NA_real_, 4), c(2, 4, 4, 6), c(1, 3, 3, 5))
stopifnot(ci_counts[["coverage_successes"]] == 3, ci_counts[["n_total"]] == 4,
  ci_counts[["available_intervals"]] == 3,
  ci_counts[["conditional_successes"]] == 3,
  ci_counts[["conditional_trials"]] == 3)
ucl_counts <- ucl_metrics(c(2, 4, NA_real_, 5), c(1, 5, 2, 7))
stopifnot(ucl_counts[["ucl_successes"]] == 1, ucl_counts[["n_total"]] == 4,
  ucl_counts[["available_ucls"]] == 3)
count_based_decision <- classify_inference_coverage(ci_counts[["coverage_successes"]],
  ci_counts[["n_total"]], ci_counts[["finite_interval_rate"]])
stopifnot(count_based_decision$successes == 3L, count_based_decision$trials == 4L)

set.seed(731)
n <- 80L
loq <- stats::qlnorm(0.40)
latent <- stats::rlnorm(n)
censored <- latent <= loq
observed <- data.frame(value = ifelse(censored, loq, latent), censored = censored,
  loq = rep(loq, n))

# Verify the survreg covariance order/transformation, the log-delta interval,
# and a distinct one-sided UCL.
mle <- apply_mle_method(observed)
mle_inf <- infer_mean(mle, observed, conf.level = 0.90, bootstrap_replicates = 29L)
stopifnot(mle_inf$inference$status %in% c("success", "partial"))
sr <- mle$details$fit@survreg
cov_names <- c("(Intercept)", "Log(scale)")
V <- stats::vcov(sr)[cov_names, cov_names, drop = FALSE]
sigma <- sr$scale
eta <- stats::coef(sr)[1L] + sigma^2 / 2
grad <- c(1, sigma^2)
se <- sqrt(drop(t(grad) %*% V %*% grad))
expected_ci <- exp(eta + stats::qnorm(c(0.05, 0.95)) * se)
expected_ucl <- exp(eta + stats::qnorm(0.90) * se)
stopifnot(isTRUE(all.equal(unname(mle_inf$inference$delta_log_wald),
  unname(expected_ci), tolerance = 1e-12)))
stopifnot(isTRUE(all.equal(unname(mle_inf$inference$ucl),
  unname(expected_ucl), tolerance = 1e-12)))
stopifnot(is.na(mle_inf$inference$ucl95))
stopifnot(mle_inf$inference$ucl != mle_inf$inference$delta_log_wald[["upper"]])
stopifnot(mle_inf$inference$bootstrap$requested == 29L)
stopifnot(mle_inf$inference$bootstrap$valid + mle_inf$inference$bootstrap$failed +
  mle_inf$inference$bootstrap$non_identifiable == 29L)
bad_covariance <- mle
bad_fit <- bad_covariance$details$fit
bad_fit@survreg$var[,] <- NA_real_
bad_covariance$details$fit <- bad_fit
bad_inf <- infer_mean(bad_covariance, observed, bootstrap_replicates = 5L)
stopifnot(bad_inf$inference$status == "failed")

# Each ROS resample passes through apply_ros_method; failures and non-identifiability
# are retained, while reconstructed values are never resampled as raw observations.
ros <- apply_ros_method(observed)
ros_boot <- bootstrap_ros(observed, B = 31L, conf.level = 0.90)
stopifnot(ros_boot$requested == 31L)
stopifnot(ros_boot$valid + ros_boot$failed + ros_boot$non_identifiable == 31L)
stopifnot(all(c("lower", "upper") %in% names(ros_boot$percentile)))
stopifnot(all(c("lower", "upper") %in% names(ros_boot$basic)))
stopifnot(is.finite(ros_boot$ucl_percentile) || ros_boot$valid < 2L)
stopifnot(is.finite(ros_boot$ucl_basic) || ros_boot$valid < 2L)
if (ros_boot$valid >= 2L) {
  alpha <- 0.10
  stopifnot(isTRUE(all.equal(unname(ros_boot$ucl_percentile),
    unname(stats::quantile(ros_boot$distribution, 0.90, names = FALSE, type = 6)), tolerance = 1e-12)))
  stopifnot(isTRUE(all.equal(unname(ros_boot$ucl_basic),
    unname(2 * ros$estimates$mean - stats::quantile(ros_boot$distribution, alpha,
      names = FALSE, type = 6)), tolerance = 1e-12)))
}

all_censored <- data.frame(value = rep(100, 20), censored = rep(TRUE, 20), loq = rep(100, 20))
ros_nonid <- bootstrap_ros_mean(all_censored, 1, B = 11L)
stopifnot(ros_nonid$valid == 0L, ros_nonid$non_identifiable == 11L)
mle_nonid <- apply_mle_method(all_censored)
stopifnot(mle_nonid$estimate_status == "non_identifiable")
mle_boot_nonid <- bootstrap_mle_mean(all_censored, meanlog = 0, sdlog = 1, B = 11L)
stopifnot(mle_boot_nonid$valid == 0L, mle_boot_nonid$non_identifiable == 11L)

# Reverse-KM quantiles below and at the observed censoring fraction are
# non-identifiable; finite point estimates above it receive inverted-band CIs.
km <- apply_reverse_km_method(observed)
p_censored <- mean(observed$censored)
km_probs <- c(p_censored - 0.01, p_censored, 0.50, 0.75, 0.999)
km_inf <- infer_reverse_km_quantiles(km, observed, probs = km_probs,
  theoretical_censoring_probability = 0.40)
keys <- names(km_inf$inference$quantile_intervals)
stopifnot(km_inf$inference$quantile_intervals[[keys[1]]]$status == "theoretically_non_identifiable")
stopifnot(km_inf$inference$quantile_intervals[[keys[2]]]$status == "theoretically_non_identifiable")
stopifnot(km_inf$inference$mean_status == "non_identifiable")
stopifnot(km_inf$inference$quantile_intervals[[keys[3]]]$status %in%
  c("success", "non_finite_interval", "not_returned"))
stopifnot(km_inf$inference$quantile_intervals[[keys[5]]]$status %in%
  c("not_returned", "non_finite_interval"))
high_censor_obs <- data.frame(value = c(rep(loq, 44), exp(seq(log(loq * 1.01), log(9), length.out = 36))),
  censored = c(rep(TRUE, 44), rep(FALSE, 36)), loq = rep(loq, 80))
high_censor_km <- apply_reverse_km_method(high_censor_obs)
high_censor_inf <- infer_reverse_km_quantiles(high_censor_km, high_censor_obs,
  probs = 0.50, theoretical_censoring_probability = 0.40)
stopifnot(high_censor_inf$inference$quantile_intervals[["q0.5"]]$status ==
  "not_estimable_in_sample")

# The main replication path does not trigger bootstrap unless explicitly asked.
sample <- list(observed = observed, metadata = list(n = n, loq = loq))
default_run <- run_single_replication(sample, methods = "ros_log")
stopifnot(default_run$methods$ros_log$inference$status == "not_implemented")

cat("Inference deterministic checks passed.\n")
