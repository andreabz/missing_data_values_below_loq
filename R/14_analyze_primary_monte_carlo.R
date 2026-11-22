# Deterministic analysis of the completed primary Monte Carlo.
# Run from the repository root: Rscript R/14_analyze_primary_monte_carlo.R

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Package 'data.table' is required. No input files were changed.", call. = FALSE)
}

results_dir <- "results/primary_monte_carlo_lognormal_seed_20261122"
audit_report <- "reports/primary_monte_carlo_data_quality.md"
analysis_report <- "reports/primary_monte_carlo_analysis.md"
figures_dir <- "reports/figures/primary_monte_carlo"
expected_commit <- "c29e92a3ea5b0116e46dbbdfb7465e60d96c88cf"
expected_methods <- c("complete_case", "substitute_zero", "substitute_loq_half",
  "substitute_loq", "ros_log", "mle_lognormal", "reverse_km")
moment_estimands <- c("mean", "median", "sd", "cv")
reverse_estimands <- c("mean", "q0.025", "q0.1", "q0.5", "q0.9", "q0.975")
key_cols <- c("scenario_id", "replicate_id", "method", "estimand")
group_cols <- c("scenario_id", "sample_size", "loq_name", "method", "estimand")
relative_reference_tolerance <- sqrt(.Machine$double.eps)

fail <- function(message) stop(message, call. = FALSE)
assert <- function(ok, message) if (!isTRUE(ok)) fail(message)
same_numeric <- function(x, y, tolerance = 1e-12) {
  length(x) == length(y) && all(is.na(x) == is.na(y)) &&
    all(is.finite(x) == is.finite(y)) &&
    (!any(is.finite(x)) || max(abs(x[is.finite(x)] - y[is.finite(y)])) <= tolerance)
}
sha256_files <- function(paths) {
  output <- tryCatch(system2("sha256sum", paths, stdout = TRUE, stderr = FALSE),
    error = function(e) character())
  if (length(output) != length(paths)) fail("Could not calculate SHA-256 hashes for the primary inputs.")
  values <- sub("  .*", "", output)
  stats::setNames(values, normalizePath(paths, winslash = "/", mustWork = TRUE))
}
fmt <- function(x, digits = 4L) {
  ifelse(is.finite(x), formatC(x, digits = digits, format = "f"), "—")
}
md_table <- function(x, columns, labels = columns, digits = 4L) {
  stopifnot(length(columns) == length(labels))
  header <- paste0("| ", paste(labels, collapse = " | "), " |")
  rule <- paste0("|", paste(rep("---", length(labels)), collapse = "|"), "|")
  rows <- vapply(seq_len(nrow(x)), function(i) {
    values <- lapply(columns, function(nm) {
      value <- x[[nm]][i]
      if (is.numeric(value) && grepl("^(n_|rel_n|^n$|distinct_proportions)", nm) &&
          is.finite(value) && value == floor(value)) format(value, scientific = FALSE, trim = TRUE) else
        if (is.numeric(value)) fmt(value, digits) else if (is.na(value)) "—" else as.character(value)
    })
    paste0("| ", paste(unlist(values), collapse = " | "), " |")
  }, character(1))
  c(header, rule, rows)
}
row_quantiles <- function(x) {
  if (!length(x)) return(rep(NA_real_, 4L))
  unname(stats::quantile(x, c(0.025, 0.25, 0.75, 0.975), names = FALSE, type = 7))
}
summarize_errors <- function(g, error_col, ref_col, benchmark) {
  eligible_status <- !is.na(g$estimand_status) & g$estimand_status == "success" &
    g$estimate_status %in% c("success", "partial") &
    g$method_eval_status %in% c("success", "partial")
  usable_estimate <- eligible_status & is.finite(g$estimate)
  err <- g[[error_col]]
  valid_error <- usable_estimate & is.finite(err)
  e <- err[valid_error]
  estimate_values <- g$estimate[usable_estimate]
  rel_reason <- rep("used", nrow(g))
  rel_reason[!eligible_status] <- "excluded_status"
  rel_reason[eligible_status & !is.finite(g$estimate)] <- "excluded_nonfinite_estimate"
  finite_estimate <- eligible_status & is.finite(g$estimate)
  rel_reason[finite_estimate & !is.finite(g[[ref_col]])] <- "excluded_nonfinite_reference"
  finite_reference <- finite_estimate & is.finite(g[[ref_col]])
  rel_reason[finite_reference & !is.finite(err)] <- "excluded_nonfinite_error"
  rel_reason[finite_reference & abs(g[[ref_col]]) <= relative_reference_tolerance] <-
    "excluded_zero_or_near_zero_reference"
  rel_ok <- rel_reason == "used" & is.finite(err)
  relative_error <- err[rel_ok] / g[[ref_col]][rel_ok]
  q <- row_quantiles(e)
  rel_q <- row_quantiles(relative_error)
  reference_label <- benchmark
  data.table::data.table(
    benchmark = reference_label,
    scenario_id = g$scenario_id[1L], sample_size = g$sample_size[1L],
    loq_name = g$loq_name[1L], method = g$method[1L], estimand = g$estimand[1L],
    n_expected = nrow(g), n_finite = sum(is.finite(g$estimate)),
    n_usable = sum(usable_estimate), n_missing_nonfinite = sum(!is.finite(g$estimate)),
    n_missing_na = sum(is.na(g$estimate) & !is.nan(g$estimate)),
    n_nan = sum(is.nan(g$estimate)), n_infinite = sum(is.infinite(g$estimate)),
    n_finite_excluded_by_status = sum(is.finite(g$estimate) & !eligible_status),
    n_metric = length(e), n_relative = length(relative_error),
    n_relative_excluded_status = sum(rel_reason == "excluded_status"),
    n_relative_excluded_estimate = sum(rel_reason == "excluded_nonfinite_estimate"),
    n_relative_excluded_reference = sum(rel_reason == "excluded_nonfinite_reference"),
    n_relative_excluded_error = sum(rel_reason == "excluded_nonfinite_error"),
    n_relative_excluded_small_reference = sum(rel_reason == "excluded_zero_or_near_zero_reference"),
    n_theoretically_nonidentifiable = sum(g$estimand_status == "theoretically_non_identifiable", na.rm = TRUE),
    n_structurally_nonidentifiable = sum(g$estimand_status == "non_identifiable", na.rm = TRUE),
    n_not_estimable_in_sample = sum(g$estimand_status == "not_estimable_in_sample", na.rm = TRUE),
    n_not_returned = sum(g$estimand_status == "not_returned", na.rm = TRUE),
    n_computational_failure = sum(g$estimate_status %in% c("failed", "non_converged") |
      g$method_eval_status %in% c("failed", "non_converged"), na.rm = TRUE),
    availability_percent = 100 * sum(usable_estimate) / nrow(g),
    error_mean = if (length(e)) mean(e) else NA_real_,
    error_median = if (length(e)) stats::median(e) else NA_real_,
    error_sd = if (length(e) > 1L) stats::sd(e) else NA_real_,
    rmse = if (length(e)) sqrt(mean(e^2)) else NA_real_,
    estimate_sd = if (benchmark == "population" && length(estimate_values) > 1L)
      stats::sd(estimate_values) else NA_real_,
    mean_relative_error = if (length(relative_error)) mean(relative_error) else NA_real_,
    median_relative_error = if (length(relative_error)) stats::median(relative_error) else NA_real_,
    error_q025 = q[1L], error_q25 = q[2L], error_q75 = q[3L], error_q975 = q[4L],
    relative_q025 = rel_q[1L], relative_q25 = rel_q[2L],
    relative_q75 = rel_q[3L], relative_q975 = rel_q[4L]
  )
}

if (!dir.exists(results_dir)) fail(paste("Missing results directory:", results_dir))
if (!file.exists(audit_report) || !grepl("Suitable with reservations",
    paste(readLines(audit_report, warn = FALSE), collapse = " "), fixed = TRUE))
  fail("The data-quality audit report is missing or does not record the qualified audit outcome.")
head <- tryCatch(system2("git", "rev-parse HEAD", stdout = TRUE, stderr = FALSE)[1L],
  error = function(e) NA_character_)
assert(identical(head, expected_commit), paste("Unexpected HEAD:", head))
expected_scenarios <- data.table::as.data.table(expand.grid(
  sample_size = c(30L, 100L, 500L), loq_name = c("LOQ_10", "LOQ_40"),
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE))
expected_scenarios[, `:=`(scenario_index = seq_len(.N),
  scenario_id = paste0("n", sample_size, "_", loq_name),
  theoretical_censoring_probability = c(0.10, 0.40)[match(loq_name, c("LOQ_10", "LOQ_40"))])]
expected_scenarios[, loq := stats::qlnorm(theoretical_censoring_probability, 0, 1)]
scenario_order <- expected_scenarios$scenario_id
expected_rds <- c(paste0(expected_scenarios$scenario_id, ".rds"), "results.rds", "run_manifest.rds")
present_rds <- list.files(results_dir, pattern = "\\.rds$", full.names = FALSE)
assert(setequal(present_rds, expected_rds), "Expected six checkpoints, results.rds, and run_manifest.rds only.")
input_paths <- file.path(results_dir, expected_rds)
assert(all(file.exists(input_paths)), "One or more expected RDS inputs are missing.")
input_hashes_before <- sha256_files(input_paths)
manifest <- readRDS(file.path(results_dir, "run_manifest.rds"))
results <- data.table::as.data.table(readRDS(file.path(results_dir, "results.rds")))
config <- manifest$signature$config
assert(isTRUE(manifest$complete) && identical(as.integer(manifest$total_completed_replicates), 12000L),
  "Manifest does not mark all 12,000 scenario replicates complete.")
assert(identical(manifest$git_commit, expected_commit), "Manifest commit differs from expected commit.")
assert(identical(as.integer(manifest$seed), 20261122L) &&
  identical(as.integer(manifest$signature$seed), 20261122L), "Manifest seed is not 20261122L.")
assert(identical(config$distribution, "lognormal") && identical(config$meanlog, 0) &&
  identical(config$sdlog, 1) && identical(as.integer(config$sample_sizes), c(30L, 100L, 500L)) &&
  identical(as.integer(config$n_monte_carlo), 2000L) &&
  identical(manifest$methods, expected_methods) && identical(manifest$signature$inference, FALSE),
  "Manifest design, methods, or inference state differs from the approved design.")
assert(identical(manifest$source_md5, manifest$signature$source_md5) &&
  all(file.exists(names(manifest$source_md5))), "Manifest source hashes are missing or inconsistent.")
current_source_md5 <- tools::md5sum(names(manifest$source_md5))
assert(identical(unname(current_source_md5), unname(manifest$source_md5)),
  "Current simulation source hashes differ from the completed-run manifest.")
assert(all(manifest$completed_by_scenario[scenario_order] == 2000L),
  "Manifest does not report 2,000 replicates in each scenario.")

required_columns <- c(key_cols, "estimand_type", "probability", "estimate", "sample_truth",
  "population_target", "estimand_status", "sample_size", "loq_name", "loq",
  "theoretical_censoring_probability", "scenario_index", "replicate_seed", "n_censored",
  "observed_censoring_proportion", "error_vs_sample_truth", "error_vs_population",
  "estimate_status", "inference_enabled", "warnings", "error_message", "method_eval_status")
missing_columns <- setdiff(required_columns, names(results))
assert(!length(missing_columns), paste("Missing required result fields:", paste(missing_columns, collapse = ", ")))
assert(all(results$inference_enabled %in% FALSE), "Results include rows with inference enabled.")
assert(!anyDuplicated(results, by = key_cols), "Duplicate result keys were found.")

checkpoint_parts <- vector("list", length(scenario_order))
names(checkpoint_parts) <- scenario_order
for (i in seq_along(scenario_order)) {
  id <- scenario_order[i]
  cp <- readRDS(file.path(results_dir, paste0(id, ".rds")))
  scenario <- expected_scenarios[scenario_id == id]
  cp_data <- data.table::as.data.table(cp$results)
  assert(identical(cp$scenario_id, id) && identical(cp$signature, manifest$signature),
    paste("Checkpoint signature/scenario mismatch:", id))
  assert(nrow(cp_data) == 60000L && identical(sort(unique(as.integer(cp_data$replicate_id))), 1:2000),
    paste("Checkpoint replicate coverage mismatch:", id))
  assert(all(cp_data$scenario_id == id) && all(cp_data$scenario_index == scenario$scenario_index) &&
    all(cp_data$sample_size == scenario$sample_size) && all(cp_data$loq_name == scenario$loq_name) &&
    all(abs(cp_data$loq - scenario$loq) < 1e-12), paste("Checkpoint scenario fields mismatch:", id))
  checkpoint_parts[[id]] <- cp_data
}
checkpoint_results <- data.table::rbindlist(checkpoint_parts, use.names = TRUE, fill = TRUE)
data.table::setorderv(checkpoint_results, key_cols)
sorted_results <- data.table::copy(results)
data.table::setorderv(sorted_results, key_cols)
assert(nrow(results) == 360000L &&
  isTRUE(all.equal(as.data.frame(checkpoint_results), as.data.frame(sorted_results),
    check.attributes = FALSE, tolerance = 0)),
  "results.rds does not match checkpoint payloads record by record.")
assert(setequal(unique(results$scenario_id), scenario_order) &&
  setequal(unique(results$method), expected_methods), "Unexpected/missing scenarios or methods in results.")
expected_keys <- data.table::rbindlist(lapply(scenario_order, function(id) {
  data.table::rbindlist(lapply(expected_methods, function(method) {
    estimands <- if (method == "reverse_km") reverse_estimands else moment_estimands
    data.table::data.table(scenario_id = id, replicate_id = rep(1:2000, each = length(estimands)),
      method = method, estimand = rep(estimands, times = 2000L))
  }))
}))
observed_keys <- unique(results[, ..key_cols])
data.table::setkeyv(expected_keys, key_cols)
data.table::setkeyv(observed_keys, key_cols)
assert(nrow(data.table::fsetdiff(expected_keys, observed_keys)) == 0L &&
  nrow(data.table::fsetdiff(observed_keys, expected_keys)) == 0L,
  "Scenario/replicate/method/estimand combinations do not match the expected design.")

theoretical_p <- c(0.10, 0.40)
loq_by_name <- stats::setNames(stats::qlnorm(theoretical_p, 0, 1), c("LOQ_10", "LOQ_40"))
assert(all(abs(results$loq - unname(loq_by_name[results$loq_name])) < 1e-12),
  "LOQ values vary from the approved fixed scenario LOQs.")
assert(all(abs(results$observed_censoring_proportion - results$n_censored / results$sample_size) < 1e-12),
  "Censoring proportions do not match censored counts/sample sizes.")
assert(all(results$n_censored >= 0 & results$n_censored <= results$sample_size),
  "Censored counts outside [0, n] were found.")

# Stored error columns have already been validated by the data-quality audit;
# repeat the arithmetic guard before calculating metrics.
err_sample_calc <- results$estimate - results$sample_truth
err_population_calc <- results$estimate - results$population_target
ok_sample <- (is.na(err_sample_calc) & is.na(results$error_vs_sample_truth)) |
  (is.finite(err_sample_calc) & is.finite(results$error_vs_sample_truth) &
    abs(err_sample_calc - results$error_vs_sample_truth) <= 1e-12)
ok_population <- (is.na(err_population_calc) & is.na(results$error_vs_population)) |
  (is.finite(err_population_calc) & is.finite(results$error_vs_population) &
    abs(err_population_calc - results$error_vs_population) <= 1e-12)
assert(all(ok_sample) && all(ok_population), "Stored errors are inconsistent with their reference columns.")

results[, eligible_for_accuracy := !is.na(estimand_status) & estimand_status == "success" &
  estimate_status %in% c("success", "partial") & method_eval_status %in% c("success", "partial") &
  is.finite(estimate)]
# A finite reverse-KM number is not eligible unless the row-specific status is success.
results[method == "reverse_km" & estimand_status != "success", eligible_for_accuracy := FALSE]

group_key <- do.call(paste, c(lapply(results[, ..group_cols], as.character), sep = "\r"))
groups <- lapply(base::split(seq_len(nrow(results)), group_key), function(rows) results[rows])
metric_parts <- vector("list", length(groups) * 2L)
availability_parts <- vector("list", length(groups))
idx <- 0L
for (g in groups) {
  idx <- idx + 1L
  sample_metric <- summarize_errors(g, "error_vs_sample_truth", "sample_truth", "sample_truth")
  population_metric <- summarize_errors(g, "error_vs_population", "population_target", "population")
  metric_parts[[2L * idx - 1L]] <- sample_metric
  metric_parts[[2L * idx]] <- population_metric
  availability_parts[[idx]] <- data.table::data.table(
    scenario_id = g$scenario_id[1L], sample_size = g$sample_size[1L], loq_name = g$loq_name[1L],
    method = g$method[1L], estimand = g$estimand[1L], n_expected = nrow(g),
    n_finite = sum(is.finite(g$estimate)),
    n_usable = sum(g$eligible_for_accuracy),
    n_missing_nonfinite = sum(!is.finite(g$estimate)),
    n_missing_na = sum(is.na(g$estimate) & !is.nan(g$estimate)),
    n_nan = sum(is.nan(g$estimate)), n_infinite = sum(is.infinite(g$estimate)),
    n_finite_excluded_by_status = sum(is.finite(g$estimate) & !g$eligible_for_accuracy),
    n_theoretically_nonidentifiable = sum(g$estimand_status == "theoretically_non_identifiable", na.rm = TRUE),
    n_structurally_nonidentifiable = sum(g$estimand_status == "non_identifiable", na.rm = TRUE),
    n_not_estimable_in_sample = sum(g$estimand_status == "not_estimable_in_sample", na.rm = TRUE),
    n_not_returned = sum(g$estimand_status == "not_returned", na.rm = TRUE),
    n_computational_failure = sum(g$estimate_status %in% c("failed", "non_converged") |
      g$method_eval_status %in% c("failed", "non_converged"), na.rm = TRUE),
    availability_percent = 100 * sum(g$eligible_for_accuracy) / nrow(g),
    rel_n_sample = sample_metric$n_relative,
    rel_excl_sample_status = sample_metric$n_relative_excluded_status,
    rel_excl_sample_estimate = sample_metric$n_relative_excluded_estimate,
    rel_excl_sample_reference = sample_metric$n_relative_excluded_reference,
    rel_excl_sample_error = sample_metric$n_relative_excluded_error,
    rel_excl_sample_near_zero = sample_metric$n_relative_excluded_small_reference,
    rel_n_population = population_metric$n_relative,
    rel_excl_population_status = population_metric$n_relative_excluded_status,
    rel_excl_population_estimate = population_metric$n_relative_excluded_estimate,
    rel_excl_population_reference = population_metric$n_relative_excluded_reference,
    rel_excl_population_error = population_metric$n_relative_excluded_error,
    rel_excl_population_near_zero = population_metric$n_relative_excluded_small_reference
  )
}
metrics <- data.table::rbindlist(metric_parts, use.names = TRUE, fill = TRUE)
availability <- data.table::rbindlist(availability_parts, use.names = TRUE, fill = TRUE)
metrics[, scenario_sort := match(scenario_id, scenario_order)]
availability[, scenario_sort := match(scenario_id, scenario_order)]
data.table::setorderv(metrics, c("scenario_sort", "benchmark", "method", "estimand"))
data.table::setorderv(availability, c("scenario_sort", "method", "estimand"))
metrics[, scenario_sort := NULL]
availability[, scenario_sort := NULL]

# Verify the Reverse-KM exclusions directly against the known audit caveat.
reverse_metrics <- metrics[method == "reverse_km"]
assert(all(reverse_metrics[estimand == "mean", n_metric] == 0L) &&
  all(reverse_metrics[estimand %in% c("q0.025", "q0.1"), n_metric] == 0L) &&
  all(reverse_metrics[n_theoretically_nonidentifiable > 0L & estimand %in% c("q0.025", "q0.1"),
    n_metric] == 0L),
  "Theoretically non-identifiable Reverse-KM rows entered the primary accuracy metrics.")
assert(all(metrics$n_expected == 2000L), "Metric denominators are not 2,000 per scenario/method/estimand.")

# Sample-truth distribution uses one shared reference per replicate/statistic,
# not repeated method rows. Reverse-KM contributes its five quantile references.
truth_source <- data.table::rbindlist(list(
  results[method == "complete_case" & estimand %in% moment_estimands],
  results[method == "reverse_km" & estimand %in% reverse_estimands[-1L]]), use.names = TRUE)
truth_parts <- truth_source[, .(replicate_id, sample_truth, population_target),
  by = .(scenario_id, sample_size, loq_name, estimand)]
assert(!anyDuplicated(truth_parts, by = c("scenario_id", "replicate_id", "estimand")),
  "Sample-truth references are duplicated after selecting the single canonical method row.")
theoretical_moments <- c(mean = exp(0.5), median = 1,
  sd = sqrt((exp(1) - 1) * exp(1)), cv = sqrt(exp(1) - 1))
truth_summary <- truth_parts[, {
  x <- sample_truth[is.finite(sample_truth)]
  q <- row_quantiles(x)
  target <- population_target[which(is.finite(population_target))[1L]]
  .(n = length(x), sample_truth_mean = if (length(x)) mean(x) else NA_real_,
    sample_truth_sd = if (length(x) > 1L) stats::sd(x) else NA_real_,
    sample_truth_median = if (length(x)) stats::median(x) else NA_real_,
    truth_q025 = q[1L], truth_q25 = q[2L], truth_q75 = q[3L], truth_q975 = q[4L],
    population_target = target,
    mean_truth_minus_population = if (length(x) && is.finite(target)) mean(x) - target else NA_real_,
    rmse_truth_vs_population = if (length(x) && is.finite(target)) sqrt(mean((x - target)^2)) else NA_real_)
}, by = .(scenario_id, sample_size, loq_name, estimand)]
assert(all(truth_summary$n == 2000L), "Sample-truth distribution lacks 2,000 references in a group.")

censor_source <- unique(results[, .(scenario_id, sample_size, loq_name, replicate_id,
  theoretical_censoring_probability, loq, n_censored, observed_censoring_proportion)])
assert(nrow(censor_source) == 12000L && !anyDuplicated(censor_source, by = c("scenario_id", "replicate_id")),
  "Censoring metadata are not unique by scenario/replicate.")
censor_summary <- censor_source[, {
  p <- observed_censoring_proportion
  q <- row_quantiles(p)
  .(n = .N, theoretical = theoretical_censoring_probability[1L], loq = loq[1L],
    mean = mean(p), median = stats::median(p), sd = stats::sd(p), min = min(p), max = max(p),
    q025 = q[1L], q25 = q[2L], q75 = q[3L], q975 = q[4L],
    distinct_proportions = data.table::uniqueN(p))
}, by = .(scenario_id, sample_size, loq_name)]
censor_summary[, scenario_sort := match(scenario_id, scenario_order)]
data.table::setorderv(censor_summary, "scenario_sort")
censor_summary[, scenario_sort := NULL]

correlation_parts <- list()
ci <- 0L
for (g in groups) {
  for (bench in c("sample_truth", "population")) {
    err_col <- if (bench == "sample_truth") "error_vs_sample_truth" else "error_vs_population"
    eligible <- g[eligible_for_accuracy & is.finite(observed_censoring_proportion)]
    eligible <- eligible[is.finite(get(err_col))]
    e <- eligible[[err_col]]
    p <- eligible$observed_censoring_proportion
    pearson <- if (length(e) >= 3L && stats::sd(e) > 0 && stats::sd(p) > 0)
      stats::cor(e, p, method = "pearson") else NA_real_
    spearman <- if (length(e) >= 3L && stats::sd(e) > 0 && stats::sd(p) > 0)
      suppressWarnings(stats::cor(e, p, method = "spearman")) else NA_real_
    ci <- ci + 1L
    correlation_parts[[ci]] <- data.table::data.table(
      benchmark = bench, scenario_id = g$scenario_id[1L], sample_size = g$sample_size[1L],
      loq_name = g$loq_name[1L], method = g$method[1L], estimand = g$estimand[1L],
      n = length(e), pearson = pearson, spearman = spearman)
  }
}
censor_correlations <- data.table::rbindlist(correlation_parts, use.names = TRUE)
censor_correlations[, scenario_sort := match(scenario_id, scenario_order)]
data.table::setorderv(censor_correlations, c("scenario_sort", "benchmark", "method", "estimand"))
censor_correlations[, scenario_sort := NULL]

if (!dir.exists(figures_dir) && !dir.create(figures_dir, recursive = TRUE))
  fail(paste("Cannot create figure directory:", figures_dir))
method_colors <- c("#3B6FB6", "#D95F02", "#1B9E77", "#7570B3", "#E7298A", "#66A61E")
moment_methods <- expected_methods[expected_methods != "reverse_km"]
figure_paths <- character()
scenario_label <- function(id) {
  row <- expected_scenarios[scenario_id == id]
  paste0("n=", row$sample_size, ", theoretical censoring=", 100 * row$theoretical_censoring_probability, "%")
}
for (id in scenario_order) {
  for (bench in c("sample_truth", "population")) {
    err_col <- if (bench == "sample_truth") "error_vs_sample_truth" else "error_vs_population"
    ref_label <- if (bench == "sample_truth") "complete-sample truth" else "theoretical population parameter"
    path <- file.path(figures_dir, paste0("moment_errors_", bench, "_", id, ".png"))
    grDevices::png(path, width = 2200, height = 1600, res = 220)
    oldpar <- graphics::par(mfrow = c(2, 2), mar = c(7.5, 4, 3.5, 1), oma = c(1.2, 0, 2, 0))
    for (stat in moment_estimands) {
      d <- results[scenario_id == id & estimand == stat & method %in% moment_methods &
        eligible_for_accuracy & is.finite(get(err_col))]
      if (nrow(d)) {
        method_n <- d[, .N, by = method]
        labels <- paste0(moment_methods, "\nn=", method_n$N[match(moment_methods, method_n$method)])
        graphics::boxplot(d[[err_col]] ~ factor(d$method, levels = moment_methods),
          names = labels, col = method_colors, outline = TRUE,
          pars = list(outpch = 16, outcex = 0.35,
            outcol = grDevices::adjustcolor("black", alpha.f = 0.40)),
          las = 2, cex.axis = 0.65,
          main = stat, xlab = "Method", ylab = paste("Error vs", ref_label))
        graphics::abline(h = 0, lty = 2, col = "grey35")
      } else {
        graphics::plot.new(); graphics::title(main = stat)
        graphics::text(0.5, 0.5, "No eligible estimates")
      }
    }
    graphics::mtext(paste("Primary Monte Carlo:", scenario_label(id), "|", ref_label),
      outer = TRUE, cex = 1.05)
    graphics::mtext("Tukey outliers shown; no observations removed; full data range.",
      side = 1, outer = TRUE, line = 0.2, cex = 0.78)
    graphics::par(oldpar)
    invisible(grDevices::dev.off())
    figure_paths <- c(figure_paths, path)
  }

  # Descriptive censoring/error association for the four common moments.
  path <- file.path(figures_dir, paste0("censoring_vs_population_error_", id, ".png"))
  grDevices::png(path, width = 2200, height = 1600, res = 220)
  oldpar <- graphics::par(mfrow = c(2, 2), mar = c(4.5, 4, 3.5, 1), oma = c(0, 0, 3, 0))
  for (stat in moment_estimands) {
    d <- results[scenario_id == id & estimand == stat & method %in% moment_methods &
      eligible_for_accuracy & is.finite(error_vs_population)]
    if (nrow(d)) {
      graphics::plot(d$observed_censoring_proportion, d$error_vs_population,
        col = grDevices::adjustcolor(method_colors[match(d$method, moment_methods)], alpha.f = 0.32),
        pch = 16, cex = 0.55, xlab = "Empirical censored proportion",
        ylab = "Error vs theoretical parameter", main = stat)
      graphics::abline(h = 0, lty = 2, col = "grey35")
      graphics::legend("topright", legend = moment_methods, col = method_colors,
        pch = 16, cex = 0.55, bty = "n")
    } else {
      graphics::plot.new(); graphics::title(main = stat)
      graphics::text(0.5, 0.5, "No eligible estimates")
    }
  }
  graphics::mtext(paste("Empirical censoring and population error:", scenario_label(id),
    "(descriptive; n=2,000 per moment/method)"), outer = TRUE, cex = 0.95)
  graphics::par(oldpar)
  invisible(grDevices::dev.off())
  figure_paths <- c(figure_paths, path)
}

# Reverse-KM quantiles are shown only where the row-specific estimand status is success.
reverse_plot_path <- file.path(figures_dir, "reverse_km_quantile_errors_vs_population.png")
grDevices::png(reverse_plot_path, width = 2400, height = 1500, res = 220)
oldpar <- graphics::par(mfrow = c(2, 3), mar = c(6.5, 4, 3.5, 1), oma = c(1.2, 0, 2, 0))
identified_quantiles <- c("q0.5", "q0.9", "q0.975")
for (id in scenario_order) {
  d <- results[scenario_id == id & method == "reverse_km" & estimand %in% identified_quantiles &
    eligible_for_accuracy & is.finite(error_vs_population)]
  if (nrow(d)) {
    counts <- d[, .N, by = estimand]
    count_labels <- counts$N[match(identified_quantiles, counts$estimand)]
    count_labels[is.na(count_labels)] <- 0L
    labels <- paste0(identified_quantiles, "\nn=", count_labels)
    graphics::boxplot(d$error_vs_population ~ factor(d$estimand, levels = identified_quantiles),
      names = labels, col = c("#3B6FB6", "#1B9E77", "#D95F02"), outline = TRUE,
      pars = list(outpch = 16, outcex = 0.35,
        outcol = grDevices::adjustcolor("black", alpha.f = 0.40)),
      las = 1, cex.axis = 0.85, main = scenario_label(id), xlab = "Reverse-KM quantile",
      ylab = "Error vs theoretical quantile")
    graphics::abline(h = 0, lty = 2, col = "grey35")
  } else {
    graphics::plot.new(); graphics::title(main = scenario_label(id))
    graphics::text(0.5, 0.5, "No eligible quantiles")
  }
}
graphics::mtext("Reverse KM: conditional on usable estimates (n shown); theoretical reference",
  outer = TRUE, cex = 1.0)
graphics::mtext("Tukey outliers shown; no observations removed; full data range.",
  side = 1, outer = TRUE, line = 0.2, cex = 0.78)
graphics::par(oldpar)
invisible(grDevices::dev.off())
figure_paths <- c(figure_paths, reverse_plot_path)

# Render detailed tables and scenario-specific figures into the report.
md <- c(
  "# Primary Monte Carlo statistical analysis", "",
  "## Objective and scope", "",
  "This report analyzes the completed primary Monte Carlo without running new simulations, bootstrap, or inference. It separates error relative to the paired complete-sample statistic from error relative to the theoretical population parameter. All summaries remain stratified by the six scenarios, method, and estimand; no overall method ranking is produced.", "",
  "The underlying data-quality audit judged the results **suitable with reservations**. This analysis applies the audit's Reverse-KM status caveat: a finite value is included in accuracy metrics only when its row-specific `estimand_status` is `success`.", "",
  "## Design and inputs", "",
  paste0("The analysis used commit `", expected_commit, "` and `", nrow(results),
    "` rows from `", results_dir, "`. The manifest confirms a lognormal population (`meanlog = 0`, `sdlog = 1`), n = 30/100/500, fixed theoretical LOQs at 10% and 40%, 2,000 replicates per condition, seven methods, and inference disabled. The six checkpoint files were re-read and compared record by record with `results.rds`; the manifest and scenario signatures matched."), "",
  paste0("Input SHA-256 hashes were recorded before analysis and recalculated after report/figure generation; all eight primary RDS hashes remained unchanged during that analysis run. A later forensic check found that the current `results.rds` hash differs from the digest transcribed in the data-quality audit. Its current logical records match the six checkpoints exactly after key ordering; the historical byte provenance remains unverified. See [the hash investigation](primary_monte_carlo_hash_investigation.md) for the evidence and limits. This does not change the analysis-versus-checkpoint content reconciliation."), "",
  "## Definitions and eligibility", "",
  "For every scenario × method × estimand, the expected denominator is 2,000. `n_finite` counts finite numeric values regardless of status; `n_usable` counts finite estimates whose row-specific estimand status is `success` (and whose method status is not a failure). Availability is `n_usable / 2,000`. Missing values, finite-but-ineligible estimates, non-identifiability, sample non-estimability, not-returned values, and computational failures are reported separately.", "",
  "Relative error is `(estimate - reference) / reference`. References with absolute value at most `sqrt(.Machine$double.eps)` are treated as zero/near-zero and excluded rather than replaced. The tables report relative-error denominators and each exclusion reason. No such exclusion is silently converted to zero.", "",
  "For sample-truth comparisons, the stored `error_vs_sample_truth` is used. For population comparisons, the separately stored `error_vs_population` is used. Bias is mean population error; relative bias is the mean population-relative error; RMSE is `sqrt(mean(error^2))`. Error quantiles use type 7 sample quantiles. Estimator SD is computed across usable estimates. Sample-truth distributions are summarized separately from method errors.", "",
  "Complete case and substitution methods summarize detected-only or constructed samples. ROS reconstructs censored values using log-ROS: the method relates normal quantiles to detected log-concentrations and predicts below-LOQ values from that regression. Its lower-tail extrapolation and reconstructed summaries depend on the adequacy of the log-scale distributional relationship; ROS is not assumption-free or equivalent to multiple imputation. In this simulation, the lognormal generator is compatible with that working relationship, so the results do not establish robustness under misspecification.", "",
  "The lognormal MLE fits the censored lognormal likelihood: detected observations contribute their density and censored observations contribute the probability of falling below the LOQ. It estimates distribution parameters and derives moments from the fitted distribution; it does not reconstruct individual censored values. Because the generating distribution is lognormal here, the MLE specification is favorable by design. Its performance should not be generalized automatically to other distributions.", "",
  "For Reverse KM, `mean` is structurally non-identifiable and is not assigned an accuracy metric. Quantiles marked `theoretically_non_identifiable`, `not_estimable_in_sample`, or `not_returned` are excluded even when a numeric value exists. Not-returned values are not necessarily computational failures. Reverse-KM accuracy metrics are conditional on a finite estimate with row status `success`; read them together with their usable counts and availability, because the eligible replicates differ by quantile and scenario. They are not directly comparable as unconditional performance summaries with methods returning estimates in all replicates. The common-LOQ results describe only quantiles that are identifiable and returned under these data conditions; no nonparametric mean is identified without additional assumptions. Reverse-KM quantile results are not pooled or ranked against moment-estimating methods.", "",
  "## Data availability and reliability", "",
  "Counts use 2,000 expected replicates for every row. `finite` is numeric availability; `usable` is the analysis-eligible count; `finite/status-excluded` records finite outputs excluded by their estimand status. Relative-exclusion counts are provided separately for each reference.", ""
)
availability_columns <- c("method", "estimand", "n_expected", "n_finite", "n_usable",
  "n_missing_nonfinite", "n_missing_na", "n_nan", "n_infinite",
  "n_finite_excluded_by_status", "n_theoretically_nonidentifiable",
  "n_structurally_nonidentifiable", "n_not_estimable_in_sample", "n_not_returned",
  "n_computational_failure", "availability_percent")
availability_labels <- c("Method", "Estimand", "Expected", "Finite", "Usable", "Missing/nonfinite",
  "NA", "NaN", "Inf", "Finite, status-excluded", "Theoretical non-ID", "Structural non-ID", "Not estimable in sample",
  "Not returned", "Computational failure", "Available (%)")
for (id in scenario_order) {
  md <- c(md, paste0("### ", scenario_label(id)), "",
    md_table(availability[scenario_id == id], availability_columns, availability_labels, digits = 2L), "")
}
md <- c(md, "### Relative-error denominators and exclusions", "",
  paste0("Relative errors exclude rows with ineligible status, non-finite estimate/reference/error, or an absolute reference no larger than ",
    format(relative_reference_tolerance, scientific = TRUE), ". Counts below are mutually exclusive by exclusion reason and separate the two references."), "")
relative_exclusion_columns <- c("method", "estimand", "rel_n_sample", "rel_excl_sample_status",
  "rel_excl_sample_estimate", "rel_excl_sample_reference", "rel_excl_sample_error",
  "rel_excl_sample_near_zero", "rel_n_population", "rel_excl_population_status",
  "rel_excl_population_estimate", "rel_excl_population_reference", "rel_excl_population_error",
  "rel_excl_population_near_zero")
relative_exclusion_labels <- c("Method", "Estimand", "N sample", "Status", "Estimate", "Reference",
  "Error", "Near-zero ref", "N population", "Status", "Estimate", "Reference", "Error", "Near-zero ref")
for (id in scenario_order) {
  md <- c(md, paste0("#### ", scenario_label(id)), "",
    md_table(availability[scenario_id == id], relative_exclusion_columns,
      relative_exclusion_labels, digits = 0L), "")
}
md <- c(md,
  "## Comparison with the complete sample", "",
  "This paired comparison uses each method's estimate minus the complete-sample statistic from the same replicate. The signed mean error is the average method-minus-sample difference; it is not population bias. Each table shows the error distribution and relative-error summaries only among eligible, finite records.", "")
metric_columns <- c("method", "estimand", "n_metric", "n_relative", "error_mean", "error_median",
  "error_sd", "rmse", "mean_relative_error", "median_relative_error", "error_q025", "error_q25",
  "error_q75", "error_q975")
metric_labels <- c("Method", "Estimand", "N used", "Relative N", "Mean error", "Median error",
  "SD error", "RMSE", "Mean relative error", "Median relative error", "Error Q2.5%", "Error Q25%",
  "Error Q75%", "Error Q97.5%")
for (id in scenario_order) {
  md <- c(md, paste0("### ", scenario_label(id)), "",
    md_table(metrics[scenario_id == id & benchmark == "sample_truth"], metric_columns,
      metric_labels, digits = 5L), "")
}

md <- c(md,
  "## Comparison with theoretical population parameters", "",
  "These are separate, unconditional-to-the-sample-reference performance summaries: population error is estimate minus the configured population parameter. Bias, relative bias, RMSE, estimator SD, error SD, median error, and error quantiles are not interchangeable with the paired sample-truth errors above.", "")
pop_columns <- c("method", "estimand", "n_metric", "n_relative", "error_mean", "mean_relative_error",
  "rmse", "estimate_sd", "error_sd", "error_median", "error_q025", "error_q25", "error_q75", "error_q975")
pop_labels <- c("Method", "Estimand", "N used", "Relative N", "Bias", "Relative bias", "RMSE",
  "SD estimates", "SD error", "Median error", "Error Q2.5%", "Error Q25%", "Error Q75%", "Error Q97.5%")
for (id in scenario_order) {
  md <- c(md, paste0("### ", scenario_label(id)), "",
    md_table(metrics[scenario_id == id & benchmark == "population"], pop_columns,
      pop_labels, digits = 5L), "")
}

md <- c(md,
  "## Sampling variability in the complete samples", "",
  "The following distributions describe the complete-sample statistics before censoring. Their spread is sampling variability under the generating population; it must not be attributed to a censoring method. `truth_mean - population` and `truth RMSE` summarize those sample statistics against their population targets, not estimator performance after censoring.", "")
truth_columns <- c("scenario_id", "estimand", "n", "population_target", "sample_truth_mean",
  "sample_truth_sd", "sample_truth_median", "truth_q025", "truth_q25", "truth_q75", "truth_q975",
  "mean_truth_minus_population", "rmse_truth_vs_population")
truth_labels <- c("Scenario", "Statistic", "N", "Population target", "Mean truth", "SD truth",
  "Median truth", "Truth Q2.5%", "Truth Q25%", "Truth Q75%", "Truth Q97.5%",
  "Mean truth − population", "Truth RMSE vs population")
md <- c(md, md_table(truth_summary, truth_columns, truth_labels, digits = 5L), "")

md <- c(md,
  "## Empirical censoring and error association", "",
  "Censoring summaries are based on one record per scenario/replicate. The empirical proportion is a realized sample characteristic, not an assigned treatment. Correlations below are descriptive associations, not causal effects. Pearson and Spearman correlations use only status-eligible finite errors and are reported separately for the sample-truth and population references.", "",
  md_table(censor_summary, c("scenario_id", "theoretical", "loq", "n", "mean", "median", "sd",
    "min", "q025", "q25", "q75", "q975", "max", "distinct_proportions"),
    c("Scenario", "Theory", "Fixed LOQ", "Replicates", "Mean", "Median", "SD", "Min",
      "Q2.5%", "Q25%", "Q75%", "Q97.5%", "Max", "Distinct values"), digits = 5L), "",
  "Correlation table columns are scenario, method, statistic, number of usable replicate pairs, Pearson r, and Spearman rho. A missing correlation means fewer than three usable pairs or no variation in one of the variables.", "")
cor_columns <- c("scenario_id", "method", "estimand", "n", "pearson", "spearman")
cor_labels <- c("Scenario", "Method", "Statistic", "N pairs", "Pearson r", "Spearman rho")
for (bench in c("sample_truth", "population")) {
  md <- c(md, paste0("### Censoring vs error: ", if (bench == "sample_truth")
    "complete-sample reference" else "population reference"), "")
  for (id in scenario_order) {
    md <- c(md, paste0("#### ", scenario_label(id)), "",
      md_table(censor_correlations[benchmark == bench & scenario_id == id],
        cor_columns, cor_labels, digits = 3L), "")
  }
}

md <- c(md,
  "## Reverse Kaplan–Meier details", "",
  "All accuracy values below are conditional on a finite estimate with `estimand_status = success`; usable counts and availability identify how many of the 2,000 expected replicates contribute. Availability differs by quantile and scenario, so these conditional accuracy metrics are not unconditional performance summaries and should be read together with the availability table. The mean has no identified nonparametric target under the common LOQ. The low quantiles are theoretically non-identifiable and are excluded even when the implementation returned a finite number. Median non-estimability and missing 97.5th quantiles are retained as availability states, not imputed or automatically treated as computational failures.", "",
  md_table(availability[method == "reverse_km"], c("scenario_id", "estimand", "n_expected",
    "n_finite", "n_usable", "n_finite_excluded_by_status", "n_theoretically_nonidentifiable",
    "n_structurally_nonidentifiable", "n_not_estimable_in_sample", "n_not_returned",
    "n_computational_failure", "availability_percent"),
    c("Scenario", "Estimand", "Expected", "Finite", "Usable", "Finite, status-excluded",
      "Theoretical non-ID", "Structural non-ID", "Not estimable", "Not returned", "Failures", "Available (%)"),
    digits = 2L), "",
  "Across the six scenarios: mean non-identifiable = 12,000; q2.5% theoretically non-identifiable = 12,000; q10% theoretically non-identifiable = 12,000; median not estimable in sample = 415; q97.5% not returned = 4,000. In addition, finite values coexist with theoretical non-identifiability for q2.5% (76 records) and q10% (2,705 records); they were excluded from every primary accuracy metric.", "")

md <- c(md, "## Figures", "",
  "Moment-error plots include complete case, the three substitution methods, ROS, and lognormal MLE. Reverse KM is kept separate because it does not estimate those reconstructed-sample moments. All 13 boxplot figures show Tukey outlier points; no observations are removed and the y-axes span the full plotted data range, with scales set independently by panel. Thus the boxplots show complete tails, while central boxes may appear compressed where extremes are large. Censoring scatterplots already display all eligible points and are unchanged. Reverse-KM plots include only identified quantiles with success status; their accuracy is conditional on the usable n shown and must be interpreted with the availability table.", "")
for (id in scenario_order) {
  md <- c(md, paste0("### ", scenario_label(id)), "",
    paste0("**Moment errors vs complete-sample truth**  \n![Moment errors vs sample truth: ", id,
      "](figures/primary_monte_carlo/moment_errors_sample_truth_", id, ".png)"), "",
    paste0("**Moment errors vs theoretical population**  \n![Moment errors vs population: ", id,
      "](figures/primary_monte_carlo/moment_errors_population_", id, ".png)"), "",
    paste0("**Empirical censoring vs population error**  \n![Censoring and population error: ", id,
      "](figures/primary_monte_carlo/censoring_vs_population_error_", id, ".png)"), "")
}
md <- c(md, "### Reverse-KM identified quantiles", "",
  "Quantile-error distributions include only rows with a usable, status-success estimate. The n below each box is the contributing replicate count; accuracy is conditional on that availability. Outlier points are shown and no usable values are removed.", "",
  "![Reverse-KM quantile errors vs population](figures/primary_monte_carlo/reverse_km_quantile_errors_vs_population.png)", "")

md <- c(md,
  "## Quantitative reading across the reported strata", "",
  "The points below describe specific strata and are not an overall method ordering:", "")
truth_mean_small <- truth_summary[scenario_id == "n30_LOQ_10" & estimand == "mean"]
truth_mean_large <- truth_summary[scenario_id == "n500_LOQ_10" & estimand == "mean"]
small_mean_rmse <- metrics[scenario_id == "n30_LOQ_10" & benchmark == "population" &
  estimand == "mean" & method != "reverse_km"]
high_censor_large_mean <- metrics[scenario_id == "n500_LOQ_40" & benchmark == "population" &
  estimand == "mean" & method != "reverse_km"]
high_bias <- stats::setNames(high_censor_large_mean$error_mean, high_censor_large_mean$method)
high_rmse <- stats::setNames(high_censor_large_mean$rmse, high_censor_large_mean$method)
md <- c(md,
  paste0("- For the complete-sample mean at theoretical 10% censoring, the SD across generated samples was ",
    fmt(truth_mean_small$sample_truth_sd), " at n=30 and ", fmt(truth_mean_large$sample_truth_sd),
    " at n=500. This sampling spread is separate from each method's error."),
  paste0("- At n=30 and theoretical 10% censoring, population-reference RMSE for the six moment-output methods ranged from ",
    fmt(min(small_mean_rmse$rmse)), " to ", fmt(max(small_mean_rmse$rmse)),
    " for the mean; the complete-sample mean itself had SD ", fmt(truth_mean_small$sample_truth_sd),
    ". Read this beside the bias and paired-error tables rather than treating RMSE as censoring-only error."),
  paste0("- At n=500 and theoretical 40% censoring, population-reference mean biases (with RMSE in parentheses) were: complete case ",
    fmt(high_bias[["complete_case"]]), " (", fmt(high_rmse[["complete_case"]]), "); zero substitution ",
    fmt(high_bias[["substitute_zero"]]), " (", fmt(high_rmse[["substitute_zero"]]), "); LOQ/2 ",
    fmt(high_bias[["substitute_loq_half"]]), " (", fmt(high_rmse[["substitute_loq_half"]]), "); LOQ ",
    fmt(high_bias[["substitute_loq"]]), " (", fmt(high_rmse[["substitute_loq"]]), "); ROS ",
    fmt(high_bias[["ros_log"]]), " (", fmt(high_rmse[["ros_log"]]), "); and lognormal MLE ",
    fmt(high_bias[["mle_lognormal"]]), " (", fmt(high_rmse[["mle_lognormal"]]),"). These are one scenario/statistic only; MLE moments come from the fitted distribution, not imputation."),
  paste0("- Mean empirical censoring was ", fmt(100 * censor_summary[scenario_id == "n30_LOQ_10", mean]),
    "% (SD ", fmt(100 * censor_summary[scenario_id == "n30_LOQ_10", sd]), " percentage points) at n=30 and ",
    fmt(100 * censor_summary[scenario_id == "n500_LOQ_10", mean]),
    "% (SD ", fmt(100 * censor_summary[scenario_id == "n500_LOQ_10", sd]),
    " percentage points) at n=500; the realized proportion varies by replicate and its error associations are descriptive."), "")

md <- c(md,
  "## Limitations", "",
  "The audited result files do not retain the latent observations. Therefore it is not possible to independently regenerate the complete-sample `sample_truth` statistics for each replicate. Their stored values agree across the six applicable methods, the theoretical targets are consistent with the configured lognormal distribution, and the stored errors equal estimate minus their designated reference on every record; this verifies internal arithmetic but not the original latent-sample calculation.", "",
  "The audit classified the Monte Carlo outputs as **suitable with reservations**. Those reservations remain: Reverse-KM numeric availability is not equivalent to estimand identifiability, its accuracy summaries are conditional on usable estimates, and sample-truth provenance cannot be independently revalidated from the retained files. The current combined result content matches the checkpoint payloads, but the historical SHA-256 discrepancy is unresolved and does not establish byte-level provenance. Inference, confidence intervals, UCLs, and bootstrap were not performed. Results describe this single lognormal design and do not establish performance under other generating distributions.", "",
  "## Conclusions", "",
  "The stored results are complete and internally reconciled for all six design conditions. The two error references yield distinct questions and are reported separately. Tables show the metric values by scenario, method, and estimand; sampling distributions of the complete-sample statistics are also shown separately, so sample variability is not attributed to a censoring method. The availability tables and figures make clear where a method has fewer usable estimates.", "",
  "The analysis supports interpretation of these scenario-specific summaries with the stated reservations. Reverse-KM records marked non-identifiable, not estimable, or not returned are excluded; its accuracy is conditional on the remaining usable replicates and is not an unconditional cross-method comparison. No generic cross-method ranking is warranted. No further simulation or inferential claim is made.", "")

writeLines(md, analysis_report, useBytes = TRUE)
input_hashes_after <- sha256_files(input_paths)
assert(identical(input_hashes_before, input_hashes_after),
  "Primary RDS hashes changed during analysis; stop and preserve the current state.")

cat("Primary Monte Carlo analysis completed.\n")
cat("Scenarios analyzed:", paste(scenario_order, collapse = ", "), "\n")
cat("Metric rows:", nrow(metrics), "(two separate references)\n")
cat("Availability rows:", nrow(availability), "\n")
cat("Censoring correlation rows:", nrow(censor_correlations), "\n")
cat("Figures written:", length(figure_paths), "to", figures_dir, "\n")
cat("Report:", analysis_report, "\n")
cat("Primary input SHA-256 values unchanged:", identical(input_hashes_before, input_hashes_after), "\n")

invisible(list(metrics = metrics, availability = availability, truth_summary = truth_summary,
  censor_summary = censor_summary, censor_correlations = censor_correlations,
  figures = figure_paths, report = analysis_report))
