# Method-specific inference. This layer is opt-in and never runs a bootstrap
# implicitly in the primary Monte Carlo path.

infer_mean <- function(method_result, observed = NULL, conf.level = 0.95,
                       bootstrap_replicates = 999L, enabled = TRUE) {
  stopifnot(is.list(method_result), !is.null(method_result$method))
  validate_confidence_level(conf.level)
  if (!isTRUE(enabled)) {
    method_result$inference$status <- if (method_result$method == "reverse_km")
      "non_identifiable" else "not_implemented"
    method_result$inference$note <- if (method_result$method == "reverse_km")
      "A nonparametric mean is not identified under a common LOQ; request quantile inference explicitly."
      else "Inference is disabled for this run; request it explicitly for validation or analysis."
    return(method_result)
  }
  if (method_result$method == "reverse_km") {
    if (is.null(observed)) stop("observed data are required for reverse-KM quantile inference")
    return(infer_reverse_km_quantiles(method_result, observed, conf.level = conf.level))
  }
  if (method_result$estimate_status != "success") {
    method_result$inference$status <- if (method_result$estimate_status == "non_identifiable")
      "non_identifiable" else "not_applicable"
    method_result$inference$note <- "Point-estimation fit is not successful; no interval was constructed."
    return(method_result)
  }
  if (is.null(observed)) stop("observed data are required for method-specific inference")
  validate_observed_data(observed)
  switch(method_result$method,
    mle_lognormal = infer_mle_mean(method_result, observed, conf.level, bootstrap_replicates),
    ros_log = infer_ros_mean(method_result, observed, conf.level, bootstrap_replicates),
    {
      method_result$inference$status <- "not_applicable"
      method_result$inference$note <- "No inferential procedure is defined for this method."
      method_result
    })
}

validate_confidence_level <- function(conf.level) {
  if (length(conf.level) != 1L || !is.numeric(conf.level) ||
      !is.finite(conf.level) || conf.level <= 0 || conf.level >= 1)
    stop("conf.level must be a finite scalar strictly between 0 and 1")
  invisible(TRUE)
}

validate_bootstrap_replicates <- function(B) {
  if (length(B) != 1L || !is.numeric(B) || !is.finite(B) || B < 2 || B != as.integer(B))
    stop("bootstrap_replicates must be an integer of at least 2")
  as.integer(B)
}

inference_interval <- function(lower, upper, conf.level, method, status = "success") {
  c(lower = lower, upper = upper)
}

infer_mle_mean <- function(method_result, observed, conf.level, bootstrap_replicates) {
  B <- validate_bootstrap_replicates(bootstrap_replicates)
  fit <- method_result$details$fit
  sr <- if (!is.null(fit) && methods::is(fit, "cenmle")) fit@survreg else NULL
  estimate <- method_result$estimates$mean
  result <- method_result
  result$inference$estimand <- "Arithmetic population mean of the fitted lognormal distribution"
  result$inference$delta_log_wald <- c(lower = NA_real_, upper = NA_real_)
  result$inference$ci_mean_level <- list(level = conf.level,
    interval = c(lower = NA_real_, upper = NA_real_))
  result$inference$parametric_bootstrap_percentile <- c(lower = NA_real_, upper = NA_real_)
  result$inference$ucl <- NA_real_
  result$inference$ucl_level <- conf.level
  result$inference$ucl95 <- NA_real_
  result$inference$ucl95_method <- "one-sided log-Wald at requested confidence level"
  if (is.null(sr) || !is.finite(estimate)) {
    result$inference$status <- "failed"
    result$inference$note <- "The NADA cenmle/survreg fit or arithmetic mean is unavailable."
    return(result)
  }
  theta <- c(meanlog = unname(stats::coef(sr)[1L]), log_sdlog = log(unname(sr$scale)))
  if (!is.finite(theta[["log_sdlog"]]) || exp(theta[["log_sdlog"]]) <= sqrt(.Machine$double.eps)) {
    result$inference$status <- "failed"
    result$inference$note <- "The survreg scale is at the numerical boundary; no Wald inference was reported."
    result$inference$diagnostics <- list(meanlog = theta[["meanlog"]],
      log_sdlog = theta[["log_sdlog"]], scale_boundary = TRUE)
    return(result)
  }
  V <- tryCatch(stats::vcov(sr), error = function(e) NULL)
  nm <- c("(Intercept)", "Log(scale)")
  valid <- is.matrix(V) && all(nm %in% rownames(V)) && all(nm %in% colnames(V))
  if (valid) {
    V <- V[nm, nm, drop = FALSE]
    eig <- tryCatch(eigen(V, symmetric = TRUE, only.values = TRUE)$values,
      error = function(e) NA_real_)
    valid <- all(is.finite(V)) && all(is.finite(eig)) && min(eig) > 0
  }
  if (!valid) {
    result$inference$status <- "failed"
    result$inference$note <- "The (meanlog, log_sdlog) covariance matrix is missing, non-finite, or not positive definite."
    result$inference$diagnostics <- list(covariance = V, parameter_order = nm,
      meanlog = theta[1L], log_sdlog = theta[2L])
    return(result)
  }
  sigma <- exp(theta[["log_sdlog"]])
  log_mean <- theta[["meanlog"]] + sigma^2 / 2
  gradient <- c(1, sigma^2)
  se_log_mean <- sqrt(drop(t(gradient) %*% V %*% gradient))
  if (!is.finite(se_log_mean) || se_log_mean <= 0 || !is.finite(log_mean)) {
    result$inference$status <- "failed"
    result$inference$note <- "Delta-method standard error for log arithmetic mean is invalid."
    return(result)
  }
  z2 <- stats::qnorm(1 - (1 - conf.level) / 2)
  z1 <- stats::qnorm(conf.level)
  result$inference$delta_log_wald <- c(
    lower = exp(log_mean - z2 * se_log_mean),
    upper = exp(log_mean + z2 * se_log_mean))
  result$inference$ci_mean_level <- list(level = conf.level,
    interval = result$inference$delta_log_wald)
  result$inference$ci_mean_95 <- if (conf.level == 0.95) result$inference$delta_log_wald else
    c(lower = NA_real_, upper = NA_real_)
  result$inference$ci_width <- diff(result$inference$delta_log_wald)
  result$inference$ucl <- exp(log_mean + z1 * se_log_mean)
  result$inference$ucl95 <- if (conf.level == 0.95) result$inference$ucl else NA_real_
  result$inference$delta_diagnostics <- list(parameterization = c("meanlog", "log(sdlog)"),
    covariance = V, gradient_log_mean = gradient, se_log_mean = se_log_mean,
    conf_level = conf.level, two_sided_critical = z2, one_sided_critical = z1)

  boot <- bootstrap_mle_mean(observed, theta[["meanlog"]], sigma, B, conf.level)
  result$inference$parametric_bootstrap_percentile <- boot$ci
  result$inference$parametric_bootstrap_ucl <- boot$ucl
  result$inference$ucl95_parametric_bootstrap <- if (conf.level == 0.95) boot$ucl else NA_real_
  result$inference$bootstrap <- boot[c("requested", "valid", "failed", "non_identifiable", "valid_rate")]
  result$inference$bootstrap_diagnostics <- boot[c("failure_reasons", "distribution")]
  result$inference$status <- if (all(is.finite(result$inference$delta_log_wald)) &&
      boot$valid_rate >= 0.8) "success" else "partial"
  result$inference$note <- paste0("Delta log-Wald interval and separately calibrated one-sided UCL; ",
    "parametric percentile bootstrap is a candidate, with full latent regeneration, censoring, and refitting.")
  result
}

bootstrap_mle_mean <- function(observed, meanlog, sdlog, B, conf.level = 0.95) {
  B <- validate_bootstrap_replicates(B)
  validate_observed_data(observed)
  loq <- observed$loq[[1L]]
  estimates <- rep(NA_real_, B)
  reasons <- character(B)
  for (b in seq_len(B)) {
    latent <- stats::rlnorm(nrow(observed), meanlog, sdlog)
    censored <- latent <= loq
    boot_obs <- data.frame(value = ifelse(censored, loq, latent), censored = censored,
      loq = rep(loq, length(latent)))
    z <- apply_mle_method(boot_obs, b)
    if (z$estimate_status == "success" && is.finite(z$estimates$mean)) estimates[b] <- z$estimates$mean
    else reasons[b] <- if (z$estimate_status == "non_identifiable") "non_identifiable" else
      if (!is.na(z$error)) z$error else z$estimate_status
  }
  valid <- is.finite(estimates)
  alpha <- (1 - conf.level) / 2
  ci <- if (sum(valid) >= 2L) stats::quantile(estimates[valid], c(alpha, 1 - alpha),
    names = FALSE, type = 6) else c(NA_real_, NA_real_)
  ucl <- if (sum(valid) >= 2L) unname(stats::quantile(estimates[valid], conf.level,
    names = FALSE, type = 6)) else NA_real_
  list(ci = setNames(ci, c("lower", "upper")), ucl = ucl, requested = B,
    valid = sum(valid), failed = sum(!valid & reasons != "non_identifiable"),
    non_identifiable = sum(reasons == "non_identifiable"), valid_rate = mean(valid),
    failure_reasons = sort(table(reasons[nzchar(reasons)]), decreasing = TRUE),
    distribution = estimates)
}

infer_ros_mean <- function(method_result, observed, conf.level, bootstrap_replicates) {
  boot <- bootstrap_ros_mean(observed, method_result$estimates$mean,
    validate_bootstrap_replicates(bootstrap_replicates), conf.level)
  result <- method_result
  result$inference$estimand <- "Arithmetic mean of the ROS-reconstructed distribution, used as an estimate of the population mean"
  result$inference$bootstrap_percentile <- boot$percentile
  result$inference$bootstrap_basic <- boot$basic
  result$inference$ucl_percentile <- boot$ucl_percentile
  result$inference$ucl_basic <- boot$ucl_basic
  result$inference$ucl95_percentile <- if (conf.level == 0.95) boot$ucl_percentile else NA_real_
  result$inference$ucl95_basic <- if (conf.level == 0.95) boot$ucl_basic else NA_real_
  result$inference$ci_mean_level <- list(level = conf.level, interval = boot$percentile)
  result$inference$ci_mean_95 <- if (conf.level == 0.95) boot$percentile else
    c(lower = NA_real_, upper = NA_real_)
  result$inference$ci_width <- if (all(is.finite(boot$percentile))) diff(boot$percentile) else NA_real_
  result$inference$bootstrap <- boot[c("requested", "valid", "failed", "non_identifiable", "valid_rate")]
  result$inference$bootstrap_diagnostics <- boot[c("failure_reasons", "distribution")]
  result$inference$status <- if (boot$valid_rate >= 0.8 && all(is.finite(boot$percentile)))
    "success" else "partial"
  result$inference$note <- paste0("Nonparametric unit bootstrap: each draw resamples observed (value, censoring, LOQ) rows and refits ROS; ",
    "reconstructed values are never resampled as raw data. BCa was not used because failed/non-identifiable resamples and ties make jackknife acceleration unstable.")
  result
}

bootstrap_ros_mean <- function(observed, estimate, B = 999L, conf.level = 0.95) {
  B <- validate_bootstrap_replicates(B)
  validate_confidence_level(conf.level)
  validate_observed_data(observed)
  if (!is.finite(estimate)) stop("A finite original ROS mean estimate is required")
  n <- nrow(observed)
  values <- rep(NA_real_, B)
  reasons <- character(B)
  for (b in seq_len(B)) {
    index <- sample.int(n, n, replace = TRUE)
    boot_obs <- observed[index, c("value", "censored", "loq"), drop = FALSE]
    z <- apply_ros_method(boot_obs, b)
    if (z$estimate_status == "success" && is.finite(z$estimates$mean)) values[b] <- z$estimates$mean
    else reasons[b] <- if (z$estimate_status == "non_identifiable") "non_identifiable" else
      if (!is.na(z$error)) z$error else z$estimate_status
  }
  valid_values <- values[is.finite(values)]
  valid <- length(valid_values)
  alpha <- (1 - conf.level) / 2
  q <- if (valid >= 2L) stats::quantile(valid_values,
    c(alpha, 1 - alpha, 1 - conf.level, conf.level),
    names = FALSE, type = 6) else rep(NA_real_, 4L)
  percentile <- c(lower = q[1L], upper = q[2L])
  basic <- c(lower = 2 * estimate - q[2L], upper = 2 * estimate - q[1L])
  list(percentile = percentile, basic = basic, ucl_percentile = q[4L],
    ucl_basic = 2 * estimate - q[3L], requested = B, valid = valid,
    failed = sum(!is.finite(values) & reasons != "non_identifiable"),
    non_identifiable = sum(reasons == "non_identifiable"), valid_rate = valid / B,
    failure_reasons = sort(table(reasons[nzchar(reasons)]), decreasing = TRUE),
    distribution = values)
}

infer_reverse_km_quantiles <- function(method_result, observed,
                                       probs = c(0.025, 0.10, 0.50, 0.90, 0.975),
                                       conf.level = 0.95,
                                       theoretical_censoring_probability = NULL) {
  validate_observed_data(observed)
  validate_confidence_level(conf.level)
  if (!is.numeric(probs) || !length(probs) || any(!is.finite(probs)) ||
      any(probs <= 0 | probs >= 1) || anyDuplicated(probs))
    stop("probs must be unique probabilities strictly between 0 and 1")
  if (!is.null(theoretical_censoring_probability) &&
      (length(theoretical_censoring_probability) != 1L ||
       !is.finite(theoretical_censoring_probability) ||
       theoretical_censoring_probability < 0 || theoretical_censoring_probability > 1))
    stop("theoretical_censoring_probability must be NULL or a scalar in [0, 1]")
  result <- method_result
  fit <- method_result$details$fit
  result$inference$estimand <- "Identifiable quantiles of the reverse-KM distribution above the common LOQ"
  result$inference$mean_status <- "non_identifiable"
  result$inference$ci_mean_95 <- c(lower = NA_real_, upper = NA_real_)
  result$inference$quantile_intervals <- vector("list", length(probs))
  result$inference$quantile_confidence_level <- conf.level
  names(result$inference$quantile_intervals) <- paste0("q", probs)
  p_cen <- mean(observed$censored)
  fit_for_band <- fit
  fitted_level <- tryCatch(fit@survfit$conf.int, error = function(e) NA_real_)
  if (methods::is(fit, "cenfit") && (!is.finite(fitted_level) || fitted_level != conf.level)) {
    refit <- capture_conditions(NADA::cenfit(observed$value, observed$censored,
      conf.int = conf.level, conf.type = "plain"))
    if (!inherits(refit$value, "method_error")) fit_for_band <- refit$value
  }
  statuses <- character(length(probs)); names(statuses) <- names(result$inference$quantile_intervals)
  for (i in seq_along(probs)) {
    p <- probs[i]; key <- names(statuses)[i]
    if (!is.null(theoretical_censoring_probability) &&
        p <= theoretical_censoring_probability) {
      statuses[i] <- "theoretically_non_identifiable"
      result$inference$quantile_intervals[[i]] <- list(status = statuses[i],
        estimate = NA_real_, ci = c(lower = NA_real_, upper = NA_real_),
        reason = "Probability is at or below the theoretical F(LOQ); the quantile is structurally unidentified.")
      next
    }
    if (p <= p_cen) {
      statuses[i] <- "not_estimable_in_sample"
      result$inference$quantile_intervals[[i]] <- list(status = statuses[i],
        estimate = NA_real_, ci = c(lower = NA_real_, upper = NA_real_),
        reason = "Probability is above theoretical F(LOQ) or unknown, but at/below the observed censoring fraction in this sample.")
      next
    }
    q <- tryCatch(suppressWarnings(as.numeric(quantile(fit, probs = p))[1L]),
      error = function(e) NA_real_)
    if (!is.finite(q)) {
      statuses[i] <- "not_returned"
      result$inference$quantile_intervals[[i]] <- list(status = statuses[i],
        estimate = NA_real_, ci = c(lower = NA_real_, upper = NA_real_),
        reason = "The reverse-KM point quantile was not returned finitely.")
      next
    }
    ci <- reverse_km_inverted_band_ci(fit_for_band, p, conf.level)
    finite_ci <- all(is.finite(ci)) && ci[1L] <= ci[2L]
    statuses[i] <- if (finite_ci) "success" else "non_finite_interval"
    result$inference$quantile_intervals[[i]] <- list(status = statuses[i],
      estimate = q, ci = ci,
      reason = if (finite_ci) NA_character_ else "The inverted reverse-KM confidence band does not yield two finite ordered endpoints.")
  }
  result$inference$quantile_status <- statuses
  result$inference$status <- if (all(statuses == "success")) "success" else
    if (any(statuses == "success")) "partial" else
      if (all(statuses %in% c("theoretically_non_identifiable", "non_identifiable")))
        "non_identifiable" else "failed"
  result$inference$note <- "The mean is not identified nonparametrically under a common LOQ. Theoretical non-identifiability (p <= F(LOQ)) is distinct from sample non-estimability (p <= observed censoring fraction) and a quantile not returned by the fitted curve. Quantile limits invert the plain Greenwood confidence band from NADA::cenfit; the package quantile(conf.int=TRUE) method is not used because it applies the probability vector as the ordinate when calculating its endpoints."
  result
}

reverse_km_inverted_band_ci <- function(fit, prob, conf.level = 0.95) {
  validate_confidence_level(conf.level)
  if (!methods::is(fit, "cenfit")) return(c(lower = NA_real_, upper = NA_real_))
  sf <- fit@survfit
  if (!is.null(sf$strata) || !all(c("time", "surv", "lower", "upper") %in% names(sf)))
    return(c(lower = NA_real_, upper = NA_real_))
  # Reverse KM's transformed survival is the original CDF. Invert the upper
  # band for the lower quantile endpoint and the lower band for the upper.
  invert <- function(curve, p) {
    keep <- is.finite(curve) & is.finite(sf$time)
    curve <- curve[keep]; time <- sf$time[keep]
    if (!length(curve) || any(diff(curve) > 1e-10)) return(NA_real_)
    idx <- which(curve >= p)
    if (!length(idx)) return(NA_real_)
    j <- max(idx)
    if (j >= length(time)) time[j] else time[j + 1L]
  }
  lower <- invert(sf$upper, prob)
  upper <- invert(sf$lower, prob)
  c(lower = lower, upper = upper)
}
