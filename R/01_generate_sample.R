# Generate one left-censored sample --------------------------------------
#
# Latent values are returned separately for simulation evaluation.
# The observed table contains the censoring threshold and an explicit
# censoring indicator. A threshold value is not a detected concentration.

generate_left_censored_sample <- function(
  n,
  loq,
  meanlog = 0,
  sdlog = 1
) {
  stopifnot(
    length(n) == 1L,
    is.finite(n),
    n >= 1,
    n == as.integer(n),
    length(loq) == 1L,
    is.finite(loq),
    loq > 0,
    length(meanlog) == 1L,
    is.finite(meanlog),
    length(sdlog) == 1L,
    is.finite(sdlog),
    sdlog > 0
  )

  n <- as.integer(n)

  x_true <- stats::rlnorm(
    n = n,
    meanlog = meanlog,
    sdlog = sdlog
  )

  censored <- x_true < loq

  observed <- data.table::data.table(
    id = seq_len(n),
    value = ifelse(censored, loq, x_true),
    censored = censored,
    loq = loq
  )

  list(
    observed = observed,
    truth = x_true,
    metadata = list(
      n = n,
      loq = loq,
      theoretical_censoring_probability = stats::plnorm(
        loq,
        meanlog = meanlog,
        sdlog = sdlog
      ),
      n_censored = sum(censored),
      observed_censoring_proportion = mean(censored)
    )
  )
}
