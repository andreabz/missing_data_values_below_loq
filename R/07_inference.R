# Method-specific inference interfaces ---------------------------------------

infer_mean <- function(method_result, ...) {
  stopifnot(is.list(method_result), !is.null(method_result$method))
  # No generic CI is imposed: ROS needs a validated full-refit bootstrap,
  # MLE needs a validated interval for the arithmetic mean, and reverse KM
  # does not identify the mean below a common LOQ.
  method_result$inference$status <- switch(method_result$method,
    complete_case = "not_implemented",
    substitute_zero = "not_implemented",
    substitute_loq_half = "not_implemented",
    substitute_loq = "not_implemented",
    ros_log = "not_implemented",
    mle_lognormal = "not_implemented",
    reverse_km = "non_identifiable",
    "not_applicable")
  method_result$inference$note <- switch(method_result$method,
    ros_log = "Pending bootstrap that refits the full ROS procedure.",
    mle_lognormal = "Pending methodological validation of the mean interval.",
    reverse_km = "Mean is not identified nonparametrically below a common LOQ.",
    "Method-specific interval procedure not yet implemented.")
  method_result
}
