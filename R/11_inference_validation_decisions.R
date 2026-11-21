# Conservative classification of empirical interval coverage.
# A 90% exact binomial interval implements two one-sided 5% equivalence tests
# against the prespecified nominal +/- tolerance band.
classify_inference_coverage <- function(successes, trials, availability = 1,
                                       nominal = 0.95, tolerance = 0.02,
                                       min_availability = 0.95, alpha = 0.05,
                                       stability_confirmed = TRUE) {
  valid_counts <- length(successes) == 1L && length(trials) == 1L &&
    is.finite(successes) && is.finite(trials) && trials >= 1 &&
    successes >= 0 && successes <= trials && successes == as.integer(successes) &&
    trials == as.integer(trials)
  if (!valid_counts || !is.finite(availability) || availability < 0 || availability > 1)
    return(list(status = "inconclusive", coverage = NA_real_, lower = NA_real_,
      upper = NA_real_, undercoverage_p = NA_real_, availability = availability))

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
    undercoverage_p = unname(exact$p.value), availability = availability)
}
