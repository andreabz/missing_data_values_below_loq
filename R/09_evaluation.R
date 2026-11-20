# Evaluation is the only layer that receives latent values --------------------

evaluate_replication <- function(replication, truth, population_parameters = NULL) {
  if (!inherits(replication, "single_replication_result")) stop("replication has invalid class")
  if (!is.numeric(truth) || length(truth) != nrow(replication$observed))
    stop("truth must be numeric and match the observed sample length")
  truth_stats <- list(mean = mean(truth), median = stats::median(truth),
    sd = stats::sd(truth), cv = stats::sd(truth) / mean(truth))
  method_rows <- lapply(replication$methods, function(x) {
    list(replicate_id = x$replicate_id, method = x$method,
      estimate_status = x$estimate_status, inference_status = x$inference$status,
      estimates = x$estimates, warnings = x$warnings, error = x$error,
      mean_error_vs_sample_truth = if (length(x$estimates$mean) == 1L &&
        is.finite(x$estimates$mean)) x$estimates$mean - truth_stats$mean else NA_real_)
  })
  list(replicate_id = replication$replicate_id,
    sample_truth = truth_stats, population_parameters = population_parameters,
    methods = method_rows)
}

collect_replication_results <- function(replications) {
  if (!is.list(replications) || any(!vapply(replications, inherits, logical(1),
      what = "single_replication_result"))) stop("All inputs must be single_replication_result objects")
  unlist(lapply(replications, `[[`, "methods"), recursive = FALSE)
}
