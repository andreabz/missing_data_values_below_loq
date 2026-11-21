# Contained simulation of candidate inference procedures. It uses the
# protocol's lognormal model, n, and absolute LOQs, but never reads n_monte_carlo.
for (file in c("R/02_method_results.R", "R/04_method_ros.R", "R/05_method_mle.R",
               "R/06_method_reverse_km.R", "R/07_inference.R",
               "R/11_inference_validation_decisions.R")) source(file)
ensure_censoring_packages_attached()

env_integer <- function(name, default) {
  x <- Sys.getenv(name, unset = as.character(default)); as.integer(x)
}
detected_cores <- parallel::detectCores(logical = TRUE)
if (!is.finite(detected_cores) || detected_cores < 1L) detected_cores <- 1L
settings <- list(seed = 20261121L,
  n_simulations = env_integer("INFERENCE_VALIDATION_N", 500L),
  bootstrap_replicates = env_integer("INFERENCE_BOOTSTRAP_B", 99L),
  stability_replicates = env_integer("INFERENCE_STABILITY_B", 199L),
  stability_cases_per_condition = env_integer("INFERENCE_STABILITY_N", 10L),
  cores = env_integer("INFERENCE_VALIDATION_CORES", min(8L,
    max(1L, detected_cores - 1L))), conf.level = 0.95,
  sample_sizes = c(30L, 100L, 500L), censoring_probabilities = c(0.10, 0.40),
  reverse_probs = c(0.10, 0.25, 0.40, 0.50, 0.75, 0.90))
stopifnot(settings$n_simulations >= 1L, settings$bootstrap_replicates >= 2L,
  settings$stability_replicates > settings$bootstrap_replicates,
  settings$stability_cases_per_condition >= 1L, settings$cores >= 1L)
mu <- 0; sigma <- 1
loqs <- stats::qlnorm(settings$censoring_probabilities, meanlog = mu, sdlog = sigma)
true_mean <- exp(mu + sigma^2 / 2)
true_quantiles <- stats::qlnorm(settings$reverse_probs, meanlog = mu, sdlog = sigma)

one_replication <- function(id, n, loq, cp, condition_id) {
  set.seed(settings$seed + condition_id * 100000L + id)
  truth <- stats::rlnorm(n, mu, sigma)
  cens <- truth <= loq
  observed <- data.frame(value = ifelse(cens, loq, truth), censored = cens,
    loq = rep(loq, n))

  m <- apply_mle_method(observed, id)
  mi <- infer_mean(m, observed, conf.level = settings$conf.level,
    bootstrap_replicates = settings$bootstrap_replicates)
  ros <- apply_ros_method(observed, id)
  ri <- if (ros$estimate_status == "success") infer_mean(ros, observed,
    conf.level = settings$conf.level, bootstrap_replicates = settings$bootstrap_replicates) else ros
  km <- apply_reverse_km_method(observed, id, probs = settings$reverse_probs)
  ki <- infer_reverse_km_quantiles(km, observed, probs = settings$reverse_probs,
    conf.level = settings$conf.level,
    theoretical_censoring_probability = cp)
  list(n_censored = sum(cens),
    mle = list(status = m$estimate_status, point = m$estimates$mean,
      delta_ci = mi$inference$delta_log_wald, delta_ucl = mi$inference$ucl95,
      bootstrap_ci = mi$inference$parametric_bootstrap_percentile,
      bootstrap_ucl = mi$inference$ucl95_parametric_bootstrap,
      bootstrap = mi$inference$bootstrap %||% list(requested = 0L, valid = 0L,
        failed = 0L, non_identifiable = 0L, valid_rate = NA_real_)),
    ros = list(status = ros$estimate_status, point = ros$estimates$mean,
      percentile = ri$inference$bootstrap_percentile %||% c(lower = NA_real_, upper = NA_real_),
      basic = ri$inference$bootstrap_basic %||% c(lower = NA_real_, upper = NA_real_),
      ucl_percentile = ri$inference$ucl95_percentile %||% NA_real_,
      ucl_basic = ri$inference$ucl95_basic %||% NA_real_,
      bootstrap = ri$inference$bootstrap %||% list(requested = 0L,
        valid = 0L, failed = 0L, non_identifiable = 0L, valid_rate = NA_real_)),
    reverse_km = list(status = ki$inference$quantile_status,
      intervals = lapply(ki$inference$quantile_intervals, `[[`, "ci"),
      estimates = lapply(ki$inference$quantile_intervals, `[[`, "estimate")))
}

`%||%` <- function(x, y) if (is.null(x)) y else x

point_metrics <- function(x, target) {
  finite <- is.finite(x)
  e <- x[finite] - target
  c(bias = if (length(e)) mean(e) else NA_real_,
    rmse = if (length(e)) sqrt(mean(e^2)) else NA_real_,
    empirical_sd = if (length(e) > 1L) stats::sd(x[finite]) else NA_real_,
    finite_rate = mean(finite))
}

condition_results <- list()
condition_id <- 0L
for (n in settings$sample_sizes) for (j in seq_along(settings$censoring_probabilities)) {
  condition_id <- condition_id + 1L
  cp <- settings$censoring_probabilities[j]; loq <- unname(loqs[j])
  message("Validating n=", n, ", F(LOQ)=", cp, " (", settings$n_simulations, " replicates)")
  ids <- seq_len(settings$n_simulations)
  draws <- if (.Platform$OS.type == "windows" || settings$cores == 1L)
    lapply(ids, one_replication, n = n, loq = loq, cp = cp, condition_id = condition_id) else
      parallel::mclapply(ids, one_replication, n = n, loq = loq, cp = cp,
        condition_id = condition_id, mc.cores = settings$cores, mc.preschedule = TRUE,
        mc.set.seed = FALSE)
  getv <- function(path) vapply(draws, function(z) {
    tryCatch({
      x <- z
      for (p in path) x <- x[[p]]
      if (length(x) == 1L) as.numeric(x) else NA_real_
    }, error = function(e) NA_real_)
  }, numeric(1))
  all_draws <- length(draws)
  mle_point <- getv(c("mle", "point")); ros_point <- getv(c("ros", "point"))
  mle_delta_lo <- getv(c("mle", "delta_ci", "lower")); mle_delta_hi <- getv(c("mle", "delta_ci", "upper"))
  mle_boot_lo <- getv(c("mle", "bootstrap_ci", "lower")); mle_boot_hi <- getv(c("mle", "bootstrap_ci", "upper"))
  ros_per_lo <- getv(c("ros", "percentile", "lower")); ros_per_hi <- getv(c("ros", "percentile", "upper"))
  ros_basic_lo <- getv(c("ros", "basic", "lower")); ros_basic_hi <- getv(c("ros", "basic", "upper"))
  mle_delta_ucl <- getv(c("mle", "delta_ucl")); mle_boot_ucl <- getv(c("mle", "bootstrap_ucl"))
  ros_per_ucl <- getv(c("ros", "ucl_percentile")); ros_basic_ucl <- getv(c("ros", "ucl_basic"))
  true_target <- rep(true_mean, all_draws)
  reverse_metrics <- lapply(seq_along(settings$reverse_probs), function(k) {
    qname <- paste0("q", settings$reverse_probs[k])
    estimate <- vapply(draws, function(z) z$reverse_km$estimates[[qname]], numeric(1))
    interval <- lapply(draws, function(z) z$reverse_km$intervals[[qname]])
    lo <- vapply(interval, function(x) x[["lower"]], numeric(1))
    hi <- vapply(interval, function(x) x[["upper"]], numeric(1))
    applicable <- settings$reverse_probs[k] > cp
    na_metrics <- c(coverage_all = NA_real_, mcse = NA_real_, conditional_coverage = NA_real_,
      finite_interval_rate = NA_real_, mean_width = NA_real_, coverage_successes = NA_real_,
      n_total = NA_real_, available_intervals = NA_real_, conditional_successes = NA_real_,
      conditional_trials = NA_real_)
    list(prob = settings$reverse_probs[k], applicable = applicable,
      point = if (applicable) point_metrics(estimate, true_quantiles[k]) else
        c(bias = NA_real_, rmse = NA_real_, empirical_sd = NA_real_, finite_rate = 0),
      interval = if (applicable) coverage_metrics(lo, hi, rep(true_quantiles[k], all_draws)) else na_metrics,
      statuses = sort(table(vapply(draws, function(z) z$reverse_km$status[[qname]], character(1)))))
  })
  names(reverse_metrics) <- paste0("q", settings$reverse_probs)
  boot_status <- function(method) {
    x <- lapply(draws, `[[`, method)
    counts <- lapply(x, function(y) y$bootstrap %||% list(requested = 0L,
      valid = 0L, failed = 0L, non_identifiable = 0L))
    c(requested = sum(vapply(counts, `[[`, numeric(1), "requested")),
      valid = sum(vapply(counts, `[[`, numeric(1), "valid")),
      failed = sum(vapply(counts, `[[`, numeric(1), "failed")),
      non_identifiable = sum(vapply(counts, `[[`, numeric(1), "non_identifiable")))
  }
  condition_results[[paste0("n", n, "_cp", cp)]] <- list(n = n, cp = cp,
    observed_censoring_mean = mean(vapply(draws, `[[`, numeric(1), "n_censored")) / n,
    mle = list(point = point_metrics(mle_point, true_mean),
      fit_statuses = sort(table(vapply(draws, function(d) d$mle$status, character(1)))),
      delta = coverage_metrics(mle_delta_lo, mle_delta_hi, true_target),
      parametric_bootstrap = coverage_metrics(mle_boot_lo, mle_boot_hi, true_target),
      ucl_delta = ucl_metrics(mle_delta_ucl, true_target),
      ucl_parametric_bootstrap = ucl_metrics(mle_boot_ucl, true_target),
      bootstrap_counts = boot_status("mle")),
    ros = list(point = point_metrics(ros_point, true_mean),
      fit_statuses = sort(table(vapply(draws, function(d) d$ros$status, character(1)))),
      percentile = coverage_metrics(ros_per_lo, ros_per_hi, true_target),
      basic = coverage_metrics(ros_basic_lo, ros_basic_hi, true_target),
      ucl_percentile = ucl_metrics(ros_per_ucl, true_target),
      ucl_basic = ucl_metrics(ros_basic_ucl, true_target),
      bootstrap_counts = boot_status("ros")),
    reverse_km = reverse_metrics, draws = draws)
}

# Stability check on the selected bootstrap procedures with twice the primary
# B, on a reproducible subset from each condition.
stability <- list()
for (id in seq_len(condition_id)) {
  key <- names(condition_results)[id]; z <- condition_results[[key]]
  cases <- z$draws[seq_len(min(settings$stability_cases_per_condition, length(z$draws)))]
  differences <- list(mle = c(), ros = c())
  for (i in seq_along(cases)) {
    # Reconstruct the original observed sample deterministically from its seed.
    set.seed(settings$seed + id * 100000L + i)
    truth <- stats::rlnorm(z$n, mu, sigma)
    cens <- truth <= unname(loqs[match(z$cp, settings$censoring_probabilities)])
    obs <- data.frame(value = ifelse(cens, unname(loqs[match(z$cp, settings$censoring_probabilities)]), truth),
      censored = cens, loq = rep(unname(loqs[match(z$cp, settings$censoring_probabilities)]), z$n))
    set.seed(settings$seed + 9000000L + id * 10000L + i)
    mf <- apply_mle_method(obs)
    mb1 <- bootstrap_mle_mean(obs, mf$parameters$meanlog, mf$parameters$sdlog,
      settings$bootstrap_replicates, settings$conf.level)
    set.seed(settings$seed + 9100000L + id * 10000L + i)
    mb2 <- bootstrap_mle_mean(obs, mf$parameters$meanlog, mf$parameters$sdlog,
      settings$stability_replicates, settings$conf.level)
    differences$mle <- c(differences$mle, abs(mb1$ci - mb2$ci))
    set.seed(settings$seed + 9200000L + id * 10000L + i)
    rf <- apply_ros_method(obs)
    rb1 <- bootstrap_ros_mean(obs, rf$estimates$mean, settings$bootstrap_replicates, settings$conf.level)
    set.seed(settings$seed + 9300000L + id * 10000L + i)
    rb2 <- bootstrap_ros_mean(obs, rf$estimates$mean, settings$stability_replicates, settings$conf.level)
    differences$ros <- c(differences$ros, abs(rb1$percentile - rb2$percentile))
  }
  stability[[key]] <- lapply(differences, function(x) c(median_abs_endpoint_change = stats::median(x),
    p90_abs_endpoint_change = unname(stats::quantile(x, 0.90)), max_abs_endpoint_change = max(x)))
}

fmt <- function(x) if (length(x) == 1L && is.finite(x)) formatC(x, digits = 4, format = "f") else "NA"
lines <- c("# Validation of inference procedures for censored estimators", "",
  "Generated by `tests/manual/validate_inference_simulation.R`.", "",
  "## Reproducibility and settings", "",
  paste0("- R: ", R.version.string),
  paste0("- Packages: NADA ", packageVersion("NADA"), "; NADA2 ", packageVersion("NADA2"),
    "; survival ", packageVersion("survival")),
  paste0("- Seed: ", settings$seed, "; n per condition: ", settings$n_simulations,
    "; bootstrap B: ", settings$bootstrap_replicates, "; stability B: ", settings$stability_replicates,
    " on ", settings$stability_cases_per_condition, " samples per condition."),
  paste0("- Parallel workers: ", settings$cores, "; RNG is reset deterministically per condition and replicate."),
  "- The settings are local to this script. `R/00_config.R` and `n_monte_carlo` are not sourced or read.",
  "- Data-generating model: Lognormal(meanlog=0, sdlog=1); n=30, 100, 500; absolute LOQs are qlnorm(0.10) and qlnorm(0.40).",
  paste0("- Nominal level: ", settings$conf.level, ". Target Monte Carlo SE for 95% coverage is <=0.01; at n=", settings$n_simulations,
    " it is about ", fmt(sqrt(0.95 * 0.05 / settings$n_simulations)), "."),
  "", "## Procedures and formulas", "",
  "**MLE.** The arithmetic mean is exp(mu + sigma^2/2). `survreg` estimates (mu, log(sigma)); the installed covariance order is `(Intercept)`, `Log(scale)`. For eta=log(mean), gradient in that order is (1, sigma^2), so Var(eta) is g' V g. The primary candidate is the exponentiated log-Wald interval eta +/- z(0.975) SE(eta), with a separately calibrated one-sided UCL exp(eta + z(0.95) SE). The alternative is a parametric percentile bootstrap that generates new latent lognormal values, censors at the same absolute LOQ, and refits MLE each draw. Bootstrap failures and all-censored non-identifiable fits remain counted.",
  "", "**ROS.** Estimand: arithmetic mean of the distribution reconstructed by ROS, used as an estimate of the population mean. Primary bootstrap resamples whole observed (value, censoring, LOQ) rows and refits ROS in every draw. Percentile and basic intervals are evaluated; the percentile UCL is the bootstrap 0.95 quantile, while the basic UCL reflects the bootstrap 0.05 quantile around the original estimate. No modeled ROS value is treated as an original observation. BCa was excluded because ties and non-identifiable resamples make the jackknife acceleration unstable.",
  "", "**Reverse KM.** The estimand is each requested quantile with p > F(LOQ); no nonparametric mean interval is produced. Quantile bands are formed by inverting the plain Greenwood lower/upper CDF bands from `NADA::cenfit`. The installed `quantile(cenfit, conf.int=TRUE)` method was inspected and not used: its source passes the probability ordinate in place of the time coordinate while calculating its endpoints. The upper CDF band gives the lower quantile limit; the lower CDF band gives the upper limit. The simulation labels p <= known theoretical F(LOQ) as structurally non-identifiable; p > F(LOQ) but at/below the realized sample censoring fraction as not estimable in that sample; and missing point/band crossings as numerical/tail return limitations.",
  "", "## Results by sample size and LOQ", "",
  "For MLE and ROS, coverage is unconditional over all simulation replicates: missing/non-finite intervals count as uncovered. Conditional coverage among finite intervals and finite-return rates are also shown. For Reverse KM, decision coverage is conditional on a finite returned band for an identifiable quantile, while availability is the proportion of all simulation replicates returning such a band; conditional coverage is not the unconditional probability of obtaining a usable covering interval. MCSE is the binomial standard error of unconditional coverage. Width is averaged over finite intervals.", "")
for (key in names(condition_results)) {
  z <- condition_results[[key]]
  lines <- c(lines, paste0("### n=", z$n, ", F(LOQ)=", z$cp), "",
    paste0("Mean realized censoring proportion: ", fmt(z$observed_censoring_mean), "."), "",
    "| Procedure | Point bias | Point RMSE | 95% coverage (MCSE) | Conditional coverage | Finite interval rate | Mean width | IC status | UCL coverage (MCSE) | Finite UCL rate | UCL status |",
    "|---|---:|---:|---:|---:|---:|---:|---|---:|---:|---|")
  rows <- list(
    list("MLE log-Wald", z$mle$point, z$mle$delta, z$mle$ucl_delta),
    list("MLE parametric percentile", z$mle$point, z$mle$parametric_bootstrap, z$mle$ucl_parametric_bootstrap),
    list("ROS percentile", z$ros$point, z$ros$percentile, z$ros$ucl_percentile),
    list("ROS basic", z$ros$point, z$ros$basic, z$ros$ucl_basic))
  for (row in rows) {
    p <- row[[2L]]; ci <- row[[3L]]; u <- row[[4L]]
    stability_confirmed <- row[[1L]] == "MLE log-Wald"
    interval_decision <- classify_inference_coverage(
      ci[["coverage_successes"]], ci[["n_total"]],
      ci[["finite_interval_rate"]], stability_confirmed = stability_confirmed)
    ucl_decision <- classify_inference_coverage(
      u[["ucl_successes"]], u[["n_total"]],
      u[["finite_ucl_rate"]], stability_confirmed = stability_confirmed)
    lines <- c(lines, paste0("| ", row[[1L]], " | ", fmt(p[["bias"]]), " | ", fmt(p[["rmse"]]),
      " | ", fmt(ci[["coverage_all"]]), " (", fmt(ci[["mcse"]]), ") | ",
      fmt(ci[["conditional_coverage"]]), " | ", fmt(ci[["finite_interval_rate"]]), " | ",
      fmt(ci[["mean_width"]]), " | ", interval_decision$status, " | ",
      fmt(u[["coverage_all"]]), " (", fmt(u[["mcse"]]), ") | ",
      fmt(u[["finite_ucl_rate"]]), " | ", ucl_decision$status, " |"))
  }
  lines <- c(lines, "", paste0("MLE bootstrap: requested ", z$mle$bootstrap_counts[["requested"]],
    "; valid ", z$mle$bootstrap_counts[["valid"]], "; failed ", z$mle$bootstrap_counts[["failed"]],
    "; non-identifiable ", z$mle$bootstrap_counts[["non_identifiable"]], "."),
    paste0("ROS bootstrap: requested ", z$ros$bootstrap_counts[["requested"]],
      "; valid ", z$ros$bootstrap_counts[["valid"]], "; failed ", z$ros$bootstrap_counts[["failed"]],
      "; non-identifiable ", z$ros$bootstrap_counts[["non_identifiable"]], "."), "",
    paste0("Original MLE fit statuses: ", paste(names(z$mle$fit_statuses),
      as.integer(z$mle$fit_statuses), sep = ": ", collapse = ", "), "."),
    paste0("Original ROS fit statuses: ", paste(names(z$ros$fit_statuses),
      as.integer(z$ros$fit_statuses), sep = ": ", collapse = ", "), "."), "",
    paste0("Theoretical reverse-KM quantiles with p <= F(LOQ): ",
      if (any(settings$reverse_probs <= z$cp)) paste0("not identifiable: ",
        paste(paste0("q", settings$reverse_probs[settings$reverse_probs <= z$cp]), collapse = ", "))
      else "none among requested probabilities."), "",
    "Reverse-KM quantiles (error and coverage among finite estimates/intervals; finite rates use all replicates):", "",
    "| Quantile | Point bias | Point RMSE | Finite point rate | 95% coverage (MCSE) | Conditional coverage | Finite interval rate | Mean width | Conditional coverage status | Status counts |",
    "|---|---:|---:|---:|---:|---:|---:|---:|---|---|")
  for (qname in names(z$reverse_km)) {
    q <- z$reverse_km[[qname]]
    status_text <- paste(paste(names(q$statuses), as.integer(q$statuses), sep = ": "), collapse = ", ")
    if (!q$applicable) coverage_status <- "not applicable: structurally non-identifiable"
    else {
      coverage_status <- classify_reverse_km_coverage(
        q$interval[["conditional_successes"]], q$interval[["conditional_trials"]],
        q$interval[["finite_interval_rate"]], applicable = TRUE)$status
    }
    lines <- c(lines, paste0("| ", qname, " | ", fmt(q$point[["bias"]]), " | ", fmt(q$point[["rmse"]]),
      " | ", fmt(q$point[["finite_rate"]]), " | ", fmt(q$interval[["coverage_all"]]), " (",
      fmt(q$interval[["mcse"]]), ") | ", fmt(q$interval[["conditional_coverage"]]), " | ",
      fmt(q$interval[["finite_interval_rate"]]), " | ", fmt(q$interval[["mean_width"]]), " | ",
      coverage_status, " | ", status_text, " |"))
  }
  lines <- c(lines, "")
}
lines <- c(lines, "## Bootstrap B stability (verification already run)", "",
  paste0("On ", settings$stability_cases_per_condition, " samples per condition, compare endpoint changes for B=",
    settings$bootstrap_replicates, " vs B=", settings$stability_replicates, "."), "",
  "| Condition | Procedure | Median absolute endpoint change | 90th percentile | Maximum |",
  "|---|---|---:|---:|---:|")
for (key in names(stability)) for (method in names(stability[[key]])) {
  x <- stability[[key]][[method]]
  lines <- c(lines, paste0("| ", key, " | ", method, " | ", fmt(x[[1L]]), " | ", fmt(x[[2L]]), " | ", fmt(x[[3L]]), " |"))
}
lines <- c(lines, "", "## Decision rules", "",
  "Each experimental condition, method, estimand and interval type is classified separately. The coverage target is the nominal 0.95. A one-sided exact binomial test (`binom.test`, alternative = `less`, alpha = 0.05) identifies statistically significant undercoverage; that result is **not validated**. Validation also requires a two-sided exact Clopper–Pearson interval at confidence 1 - 2*alpha = 0.90 wholly within 0.93–0.97. This is the equal-tailed confidence-interval form of two one-sided equivalence tests at alpha=0.05 for the +/-0.02 margin. Availability is reported separately and must be at least 0.95 for validation. Bootstrap procedures also require demonstrated endpoint stability; otherwise they remain **inconclusive** unless undercoverage is significant. All other cases without significant undercoverage are **inconclusive**; failure to reject undercoverage alone is not evidence of validity.",
  "", "The alpha=0.05 classifications are scenario-specific and descriptive. No Holm correction is applied because these rows make no simultaneous or family-wide validation claim. The per-stratum statuses do not guarantee family-wise error control and must not be combined into a global claim that a procedure family is validated.",
  "", "Two-sided 95% mean intervals and one-sided 95% UCLs have separate decisions. Reverse-KM coverage is conditional on returning a finite interval for a theoretically identifiable quantile; its unconditional interval availability remains a separate metric. Structural non-identifiability is not assessed as coverage.",
  "", "B=99 results remain preliminary pending targeted endpoint stability verification. The targeted B=99/B=499 supplement is reported separately in `reports/inference_bootstrap_stability.md`; it does not replace or recompute the six-condition coverage simulation.",
  "", "No Gamma misspecification scenario or primary Monte Carlo was run.", "")
writeLines(lines, "reports/inference_validation.md")
cat("Wrote reports/inference_validation.md\n")
