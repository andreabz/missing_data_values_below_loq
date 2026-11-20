# Evaluation is the only layer that receives latent values --------------------

evaluate_replication <- function(replication, truth, population_parameters = NULL) {
  if (!inherits(replication, "single_replication_result")) stop("replication has invalid class")
  if (!is.numeric(truth) || length(truth) != nrow(replication$observed))
    stop("truth must be numeric and match the observed sample length")
  truth_stats <- list(mean = safe_stat(mean(truth)), median = safe_stat(stats::median(truth)),
    sd = safe_stat(stats::sd(truth)), cv = safe_stat(stats::sd(truth) / mean(truth)))
  if (is.null(population_parameters)) population_parameters <- list()
  safe_error <- function(estimate, reference) {
    if (length(estimate) != 1L || length(reference) != 1L ||
        !is.finite(estimate) || !is.finite(reference)) return(NA_real_)
    estimate - reference
  }
  method_rows <- lapply(replication$methods, function(x) {
    statistics <- c("mean", "median", "sd", "cv")
    sample_errors <- population_errors <- stats::setNames(vector("list", length(statistics)), statistics)
    for (stat in statistics) {
      estimate <- x$estimates[[stat]]
      sample_errors[[stat]] <- safe_error(estimate, truth_stats[[stat]])
      population_errors[[stat]] <- safe_error(estimate, population_parameters[[stat]])
    }
    list(replicate_id = x$replicate_id, method = x$method,
      estimate_status = x$estimate_status, inference_status = x$inference$status,
      estimates = x$estimates, warnings = x$warnings, error = x$error,
      errors_vs_sample_truth = sample_errors,
      errors_vs_population_parameters = population_errors)
  })
  list(replicate_id = replication$replicate_id,
    sample_truth = truth_stats, population_parameters = population_parameters,
    methods = method_rows)
}

safe_stat <- function(x) if (length(x) == 1L && is.finite(x)) x else NA_real_

collect_replication_results <- function(replications) {
  if (!is.list(replications) || !length(replications) || any(!vapply(replications, inherits, logical(1),
      what = "single_replication_result"))) stop("All inputs must be single_replication_result objects")
  rows <- lapply(replications, function(replication) {
    lapply(replication$methods, function(result) {
      result$replicate_id <- replication$replicate_id
      result
    })
  })
  names(rows) <- vapply(replications, function(x) as.character(x$replicate_id), character(1))
  structure(rows, class = c("replication_result_collection", "list"))
}
