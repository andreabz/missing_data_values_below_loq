# Conservative classification of empirical interval coverage.
# A 90% exact binomial interval implements two one-sided 5% equivalence tests
# against the prespecified nominal +/- tolerance band.
classify_inference_coverage <- function(successes, trials, availability = 1,
                                       nominal = 0.95, tolerance = 0.02,
                                       min_availability = 0.95, alpha = 0.05,
                                       stability_confirmed = TRUE) {
  scalar_finite_number <- function(x) is.numeric(x) && !is.logical(x) && !is.complex(x) &&
    length(x) == 1L && is.finite(x)
  if (!scalar_finite_number(successes) || !scalar_finite_number(trials) ||
      successes != floor(successes) || trials != floor(trials) ||
      successes < 0 || trials < 0 || successes > trials ||
      successes > .Machine$integer.max || trials > .Machine$integer.max)
    stop("successes and trials must be non-negative integer counts with successes <= trials")
  if (trials == 0 && successes != 0)
    stop("successes must be zero when trials is zero")
  if (!scalar_finite_number(availability) || availability < 0 || availability > 1)
    stop("availability must be a finite scalar in [0, 1]")
  if (!scalar_finite_number(nominal) || nominal <= 0 || nominal >= 1)
    stop("nominal must be a finite scalar strictly between 0 and 1")
  if (!scalar_finite_number(tolerance) || tolerance <= 0 ||
      nominal - tolerance < 0 || nominal + tolerance > 1)
    stop("tolerance must be positive and keep nominal +/- tolerance within [0, 1]")
  if (!scalar_finite_number(alpha) || alpha <= 0 || alpha >= 0.5)
    stop("alpha must be a finite scalar strictly between 0 and 0.5")
  if (!scalar_finite_number(min_availability) || min_availability < 0 || min_availability > 1)
    stop("min_availability must be a finite scalar in [0, 1]")
  if (!is.logical(stability_confirmed) || length(stability_confirmed) != 1L ||
      is.na(stability_confirmed))
    stop("stability_confirmed must be TRUE or FALSE")
  if (trials == 0) return(list(status = "inconclusive", coverage = NA_real_,
    lower = NA_real_, upper = NA_real_, undercoverage_p = NA_real_,
    successes = 0L, trials = 0L, availability = availability))

  exact <- stats::binom.test(as.integer(successes), as.integer(trials),
    p = nominal, alternative = "less")
  equivalence <- stats::binom.test(as.integer(successes), as.integer(trials),
    conf.level = 1 - 2 * alpha)$conf.int
  status <- if (exact$p.value < alpha) "not validated" else
    if (!isTRUE(stability_confirmed)) "inconclusive" else
    if (availability < min_availability) "inconclusive" else
      if (equivalence[[1L]] >= nominal - tolerance &&
          equivalence[[2L]] <= nominal + tolerance) "validated" else "inconclusive"
  list(status = status, coverage = successes / trials,
    lower = unname(equivalence[[1L]]), upper = unname(equivalence[[2L]]),
    undercoverage_p = unname(exact$p.value), successes = as.integer(successes),
    trials = as.integer(trials), availability = availability)
}

classify_reverse_km_coverage <- function(successes, trials, availability,
                                         applicable, ...) {
  if (!is.logical(applicable) || length(applicable) != 1L || is.na(applicable))
    stop("applicable must be TRUE or FALSE")
  if (!applicable) return(list(status = "not applicable", coverage = NA_real_,
    lower = NA_real_, upper = NA_real_, undercoverage_p = NA_real_,
    successes = NA_integer_, trials = NA_integer_, availability = NA_real_))
  classify_inference_coverage(successes, trials, availability, ...)
}

# Keep integer counts alongside rates so decisions never reconstruct counts
# from rounded report values.
coverage_metrics <- function(lower, upper, target) {
  if (!is.numeric(lower) || !is.numeric(upper) || !is.numeric(target) ||
      length(lower) != length(upper) || length(upper) != length(target) ||
      !length(target) || any(!is.finite(target)))
    stop("lower, upper, and target must be equal-length numeric vectors with finite targets")
  n <- length(target)
  finite <- is.finite(lower) & is.finite(upper) & lower <= upper
  covered <- finite & lower <= target & upper >= target
  rate <- mean(covered)
  available <- sum(finite)
  successes <- sum(covered)
  c(coverage_all = rate, mcse = sqrt(rate * (1 - rate) / n),
    conditional_coverage = if (available) successes / available else NA_real_,
    finite_interval_rate = available / n,
    mean_width = if (available) mean(upper[finite] - lower[finite]) else NA_real_,
    coverage_successes = successes, n_total = n, available_intervals = available,
    conditional_successes = successes, conditional_trials = available)
}

ucl_metrics <- function(ucl, target) {
  if (!is.numeric(ucl) || !is.numeric(target) || length(ucl) != length(target) ||
      !length(target) || any(!is.finite(target)))
    stop("ucl and target must be equal-length numeric vectors with finite targets")
  n <- length(target)
  finite <- is.finite(ucl)
  covered <- finite & ucl >= target
  successes <- sum(covered)
  available <- sum(finite)
  rate <- successes / n
  c(coverage_all = rate, mcse = sqrt(rate * (1 - rate) / n),
    finite_ucl_rate = available / n, mean_ucl = if (available) mean(ucl[finite]) else NA_real_,
    ucl_successes = successes, n_total = n, available_ucls = available)
}
