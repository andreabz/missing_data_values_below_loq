# Complete-case and substitution methods -------------------------------------

apply_simple_method <- function(observed, method, replicate_id = NA_integer_) {
  validate_observed_data(observed)
  choices <- c("complete_case", "substitute_zero", "substitute_loq_half", "substitute_loq")
  if (length(method) != 1L || !method %in% choices) stop("Unknown simple method")

  x <- observed$value
  censored <- observed$censored
  if (method == "complete_case") x <- x[!censored]
  if (method == "substitute_zero") x[censored] <- 0
  if (method == "substitute_loq_half") x[censored] <- observed$loq[censored] / 2
  if (method == "substitute_loq") x[censored] <- observed$loq[censored]

  estimates <- if (!length(x)) {
    list(mean = NA_real_, median = NA_real_, sd = NA_real_, cv = NA_real_)
  } else {
    m <- mean(x)
    s <- if (length(x) >= 2L) stats::sd(x) else NA_real_
    list(mean = m, median = stats::median(x), sd = s,
         cv = if (is.finite(m) && m != 0 && is.finite(s)) s / m else NA_real_)
  }
  status <- if (!length(x)) "not_applicable" else if (length(x) < 2L) "partial" else "success"
  new_method_result(replicate_id, method, observed, status,
                    estimates = estimates,
                    details = list(n_analysis = length(x), sample_type =
                      if (method == "complete_case") "uncensored_only" else "substituted"))
}
