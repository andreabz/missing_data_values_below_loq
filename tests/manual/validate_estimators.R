# Independent estimator validation; intentionally small and separate from the
# primary Monte Carlo design in R/00_config.R.
source("R/02_method_results.R")
source("R/04_method_ros.R")
source("R/05_method_mle.R")
source("R/06_method_reverse_km.R")
source("R/10_independent_validation.R")

set.seed(20261120)
validation_settings <- list(seed = 20261120L, n = 60L, replicates = 250L,
  censoring_probabilities = c(0.10, 0.40), reverse_km_replicates = 250L)
ensure_censoring_packages_attached()

observed_from <- function(x, loq) {
  censored <- x <= loq
  data.frame(value = ifelse(censored, loq, x), censored = censored,
    loq = rep(loq, length(x)))
}
latent_grid <- stats::qlnorm((seq_len(80L) - 0.5) / 80L)
mle_cases <- list(
  none = observed_from(latent_grid, stats::qlnorm(0.001)),
  partial = observed_from(latent_grid, stats::qlnorm(0.40)),
  elevated = observed_from(latent_grid, stats::qlnorm(0.75)),
  total = observed_from(latent_grid, max(latent_grid) * 1.01)
)
mle_comparisons <- lapply(mle_cases[c("none", "partial", "elevated")], function(obs) {
  independent <- fit_independent_censored_lognormal(obs)
  package <- capture_conditions(NADA::cenmle(obs$value, obs$censored,
    dist = "lognormal", conf.int = 0.95))
  if (independent$status != "success" || inherits(package$value, "method_error"))
    return(list(independent = independent, package_error =
      if (inherits(package$value, "method_error")) package$value$message else NA_character_))
  sr <- package$value@survreg
  p <- c(meanlog = unname(stats::coef(sr)[1L]), sdlog = unname(sr$scale))
  list(independent = independent, package = p,
    package_loglik = censored_lognormal_loglik(p[[1L]], log(p[[2L]]), obs),
    parameter_difference = independent$parameters - p,
    loglik_difference = independent$loglik -
      censored_lognormal_loglik(p[[1L]], log(p[[2L]]), obs),
    package_warnings = package$warnings)
})
stopifnot(all(vapply(mle_comparisons, function(z)
  !is.null(z$parameter_difference) && max(abs(z$parameter_difference)) < 1e-4 &&
    abs(z$loglik_difference) < 1e-4, logical(1))))
stopifnot(fit_independent_censored_lognormal(mle_cases$total)$status == "non_identifiable")

ros_cases <- lapply(mle_cases, function(obs) {
  z <- apply_ros_method(obs)
  list(result = z, observed = obs)
})
stopifnot(ros_cases$none$result$estimate_status == "success")
stopifnot(max(abs(ros_cases$none$result$details$fit$modeled - mle_cases$none$value)) < 1e-8)
stopifnot(ros_cases$total$result$estimate_status == "non_identifiable")
for (case in ros_cases[c("none", "partial", "elevated")]) {
  modeled <- case$result$details$fit$modeled
  stopifnot(length(modeled) == nrow(case$observed), all(is.finite(modeled)), all(modeled > 0))
  stopifnot(identical(order(modeled), order(case$observed$value)))
  stopifnot(all(modeled[!case$observed$censored] == case$observed$value[!case$observed$censored]))
}

# Simulation measures both reconstruction of the same latent sample and
# estimation of population moments; failures/non-identifiability stay in rates.
ros_sim <- list()
for (cp in validation_settings$censoring_probabilities) {
  loq <- stats::qlnorm(cp)
  draws <- lapply(seq_len(validation_settings$replicates), function(i) {
    truth <- stats::rlnorm(validation_settings$n)
    obs <- observed_from(truth, loq)
    result <- apply_ros_method(obs, i)
    if (result$estimate_status != "success") return(list(ok = FALSE, reason = result$estimate_status))
    modeled <- result$details$fit$modeled
    stats_for <- function(x) c(mean = mean(x), median = stats::median(x),
      sd = stats::sd(x), cv = stats::sd(x) / mean(x))
    list(ok = TRUE, modeled = stats_for(modeled), sample = stats_for(truth))
  })
  valid <- vapply(draws, function(x) isTRUE(x$ok), logical(1))
  pop <- unlist(lognormal_moments(0, 1))
  metrics <- list()
  for (stat in names(pop)) {
    est <- vapply(draws[valid], function(x) x$modeled[[stat]], numeric(1))
    ref <- vapply(draws[valid], function(x) x$sample[[stat]], numeric(1))
    metrics[[stat]] <- list(vs_sample = validation_summary(est - ref, 0),
      vs_population = validation_summary(est, pop[[stat]]))
  }
  ros_sim[[paste0("censor_", cp)]] <- list(probability = cp,
    success = sum(valid), failure = sum(!valid), non_identifiable = sum(vapply(draws,
      function(x) identical(x$reason, "non_identifiable"), logical(1))), metrics = metrics)
}

# Reverse-KM empirical quantile validation. Probabilities at/below F(LOQ) are
# structurally unidentified; above it, only finite returned values are scored.
reverse_probs <- c(0.10, 0.25, 0.40, 0.50, 0.75, 0.90)
reverse_sim <- list()
for (cp in validation_settings$censoring_probabilities) {
  loq <- stats::qlnorm(cp)
  draws <- lapply(seq_len(validation_settings$reverse_km_replicates), function(i) {
    truth <- stats::rlnorm(validation_settings$n)
    obs <- observed_from(truth, loq)
    fit <- capture_conditions(NADA::cenfit(obs$value, obs$censored))
    if (inherits(fit$value, "method_error")) return(list(failed = TRUE, values = rep(NA_real_, length(reverse_probs))))
    qs <- rep(NA_real_, length(reverse_probs))
    for (j in seq_along(reverse_probs)) if (reverse_probs[j] > cp) {
      q <- capture_conditions(quantile(fit$value, probs = reverse_probs[j]))
      if (!inherits(q$value, "method_error")) {
        v <- suppressWarnings(as.numeric(q$value)[1L]); if (is.finite(v)) qs[j] <- v
      }
    }
    list(failed = FALSE, values = qs)
  })
  vals <- do.call(rbind, lapply(draws, `[[`, "values"))
  identified <- reverse_probs > cp
  metrics <- lapply(which(identified), function(j) validation_summary(vals[, j],
    stats::qlnorm(reverse_probs[j])))
  names(metrics) <- paste0("q", reverse_probs[identified])
  reverse_sim[[paste0("censor_", cp)]] <- list(probability = cp,
    theoretical_non_identifiable = if (any(!identified)) paste0("q", reverse_probs[!identified]) else character(),
    failed_fits = sum(vapply(draws, `[[`, logical(1), "failed")),
    finite_return_rate = stats::setNames(vapply(which(identified),
      function(j) mean(is.finite(vals[, j])), numeric(1)),
      paste0("q", reverse_probs[identified])),
    metrics = metrics)
}

format_table <- function(x, digits = 5L) paste(capture.output(print(round(x, digits))), collapse = "\n")
lines <- c("# Independent validation of censored estimators", "",
  "Generated by `tests/manual/validate_estimators.R`.", "",
  "## Reproducibility", "",
  paste0("- R: ", R.version.string),
  paste0("- Packages: NADA ", packageVersion("NADA"), "; NADA2 ", packageVersion("NADA2"),
    "; survival ", packageVersion("survival")),
  paste0("- Seed: ", validation_settings$seed, "; simulation n = ", validation_settings$n,
    "; replicates per condition = ", validation_settings$replicates,
    "; reverse KM replicates = ", validation_settings$reverse_km_replicates),
  "- Simulation settings are local to this script and do not read or modify `n_monte_carlo`.",
  "- Data-generating distribution: Lognormal(meanlog = 0, sdlog = 1); common LOQs are its 10th and 40th percentiles.",
  "", "## Criteria", "",
  "The independent likelihood sums log densities for detections and log CDF probabilities at each censored LOQ. It optimizes meanlog and log(sdlog) by BFGS; log(sdlog) is constrained to [-8, 5], maxit=2000, reltol=1e-11. Boundary solutions and nonzero optimizer convergence codes are reported as non-converged. Agreement tolerance for parameter and log-likelihood comparisons is 1e-4 (optimizer scale and package parameterization); uncensored ROS reproduction tolerance is 1e-8 (floating-point reconstruction).",
  "", "## MLE independent likelihood comparison", "",
  "Mean, median, SD and CV are derived from the fitted lognormal parameters. Likelihood differences evaluate the same independent observed-data likelihood at each fit.", "")
for (nm in names(mle_comparisons)) {
  z <- mle_comparisons[[nm]]
  lines <- c(lines, paste0("### ", nm), "")
  if (!is.null(z$package_error)) lines <- c(lines, paste0("Fit issue: ", z$package_error), "") else {
    lines <- c(lines, paste0("Censored: ", sum(mle_cases[[nm]]$censored), "/", nrow(mle_cases[[nm]])),
      paste0("Independent parameters: ", paste(names(z$independent$parameters), signif(z$independent$parameters, 7), collapse = "; ")),
      paste0("NADA::cenmle parameters: ", paste(names(z$package), signif(z$package, 7), collapse = "; ")),
      paste0("Parameter differences (independent - NADA): ", paste(signif(z$parameter_difference, 4), collapse = ", ")),
      paste0("Independent log-likelihood: ", signif(z$independent$loglik, 8), "; at NADA fit: ",
        signif(z$package_loglik, 8), "; difference: ", signif(z$loglik_difference, 4)),
      paste0("Derived moments: ", paste(names(z$independent$moments), signif(unlist(z$independent$moments), 6), collapse = "; ")), "")
  }
}
lines <- c(lines, "All-censored sample: independently classified as non-identifiable because the likelihood increases as the location moves below the LOQ while scale is unconstrained by detections.",
  "", "## ROS deterministic checks and simulation", "",
  "Deterministic cases include no, partial, elevated (75%), and total censoring. In uncensored data, modeled values reproduce observations; in all identifiable cases, reconstructed values were checked for positivity, finiteness, length, ordering, and exact preservation of detections.", "")
for (nm in names(ros_sim)) {
  z <- ros_sim[[nm]]
  lines <- c(lines, paste0("### ", nm, " (F(LOQ)=", z$probability, ")"), "",
    paste0("Successful: ", z$success, "/", validation_settings$replicates,
      "; failures: ", z$failure, "; non-identifiable: ", z$non_identifiable), "",
    "| Statistic | Reference | Bias | RMSE | Empirical SD |", "|---|---|---:|---:|---:|")
  for (stat in names(z$metrics)) {
    a <- z$metrics[[stat]]$vs_sample; b <- z$metrics[[stat]]$vs_population
    lines <- c(lines, paste0("| ", stat, " vs same-replicate latent sample | 0 | ",
      paste(signif(a[c("bias", "rmse", "empirical_sd")], 5), collapse = " | "), " |"),
      paste0("| ", stat, " vs population parameter | ", signif(unlist(lognormal_moments(0, 1))[[stat]], 5), " | ",
      paste(signif(b[c("bias", "rmse", "empirical_sd")], 5), collapse = " | "), " |"))
  }
  lines <- c(lines, "")
}
lines <- c(lines, "## Reverse KM quantiles", "",
  "For each common LOQ, probabilities p <= F(LOQ) are listed as structurally non-identifiable and are not scored. For p > F(LOQ), quantiles are scored only when NADA returns finite values; the finite-return rate captures numerical/tail limitations separately.", "")
for (nm in names(reverse_sim)) {
  z <- reverse_sim[[nm]]
  lines <- c(lines, paste0("### ", nm, " (F(LOQ)=", z$probability, ")"), "",
    paste0("Theoretically non-identifiable: ", if (length(z$theoretical_non_identifiable))
      paste(z$theoretical_non_identifiable, collapse = ", ") else "none among requested probabilities",
      "; failed fits: ", z$failed_fits), "")
  for (q in names(z$metrics)) lines <- c(lines, paste0("- ", q, " finite return rate: ",
    signif(z$finite_return_rate[[q]], 4), "; bias/RMSE/empirical SD: ",
    paste(signif(z$metrics[[q]], 5), collapse = "/")))
  lines <- c(lines, "")
}
lines <- c(lines, "## Decisions and limitations", "",
  "- **MLE: validated with reservations.** The independently coded observed-data likelihood and package fit are compared on deterministic uncensored, partial, and elevated-censoring samples; this validates objective/point-fit agreement for these cases, not interval coverage. The optimizer uses explicit bounds to avoid degenerate scales.",
  "- **ROS: validated with reservations.** Deterministic structural checks and a small targeted simulation evaluate reconstruction and estimation separately. This is not confidence-interval validation; results depend on the ROS regression and its package implementation.",
  "- **Reverse KM: validated with reservations.** Theoretical non-identifiability below the common LOQ is separated from finite-return rates above it. No distribution below LOQ is attributed to reverse KM.",
  "", "The simulations are diagnostic (250 samples per condition), not the primary 2,000-replicate design. Monte Carlo error is material, especially at elevated censoring and for extreme statistics. No confidence intervals, UCLs, bootstrap, or inferential coverage are implemented or assessed.", "")
writeLines(lines, "reports/estimator_validation.md")
cat("Wrote reports/estimator_validation.md\n")
