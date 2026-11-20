# One replicate: methods see observed data only --------------------------------

run_single_replication <- function(sample, replicate_id = 1L,
                                   methods = c("complete_case", "substitute_zero",
                                     "substitute_loq_half", "substitute_loq",
                                     "ros_log", "mle_lognormal", "reverse_km"),
                                   inference = FALSE, inference_control = list()) {
  if (!is.list(sample) || is.null(sample$observed)) stop("sample must contain observed data")
  observed <- as.data.frame(sample$observed)
  validate_observed_data(observed)
  allowed <- c("complete_case", "substitute_zero", "substitute_loq_half",
               "substitute_loq", "ros_log", "mle_lognormal", "reverse_km")
  if (any(!methods %in% allowed)) stop("Unknown method requested")
  if (!is.logical(inference) || length(inference) != 1L || is.na(inference))
    stop("inference must be TRUE or FALSE")
  if (!is.list(inference_control)) stop("inference_control must be a list")
  if (length(inference_control) && (is.null(names(inference_control)) ||
      anyDuplicated(names(inference_control)) ||
      any(!names(inference_control) %in% c("conf.level", "bootstrap_replicates"))))
    stop("inference_control may contain only conf.level and bootstrap_replicates")

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
    do.call(infer_mean, c(list(method_result = result, observed = observed,
      enabled = inference), inference_control))
  })
  names(results) <- methods
  structure(list(replicate_id = replicate_id,
    metadata = sample$metadata,
    observed = observed,
    methods = results), class = "single_replication_result")
}
