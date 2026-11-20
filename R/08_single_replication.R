# One replicate: methods see observed data only --------------------------------

run_single_replication <- function(sample, replicate_id = 1L,
                                   methods = c("complete_case", "substitute_zero",
                                     "substitute_loq_half", "substitute_loq",
                                     "ros_log", "mle_lognormal", "reverse_km")) {
  if (!is.list(sample) || is.null(sample$observed)) stop("sample must contain observed data")
  observed <- as.data.frame(sample$observed)
  validate_observed_data(observed)
  allowed <- c("complete_case", "substitute_zero", "substitute_loq_half",
               "substitute_loq", "ros_log", "mle_lognormal", "reverse_km")
  if (any(!methods %in% allowed)) stop("Unknown method requested")

  runners <- list(
    complete_case = function() apply_simple_method(observed, "complete_case", replicate_id),
    substitute_zero = function() apply_simple_method(observed, "substitute_zero", replicate_id),
    substitute_loq_half = function() apply_simple_method(observed, "substitute_loq_half", replicate_id),
    substitute_loq = function() apply_simple_method(observed, "substitute_loq", replicate_id),
    ros_log = function() apply_ros_method(observed, replicate_id),
    mle_lognormal = function() apply_mle_method(observed, replicate_id),
    reverse_km = function() apply_reverse_km_method(observed, replicate_id)
  )
  results <- lapply(methods, function(method) {
    result <- tryCatch(runners[[method]](), error = function(e)
      new_method_result(replicate_id, method, observed, "failed", error = conditionMessage(e)))
    infer_mean(result)
  })
  names(results) <- methods
  structure(list(replicate_id = replicate_id,
    metadata = sample$metadata,
    observed = observed,
    methods = results), class = "single_replication_result")
}
