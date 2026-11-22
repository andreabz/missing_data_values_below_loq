# Read-only data-quality audit for the completed primary Monte Carlo.
# Run from the repository root with: Rscript R/13_validate_primary_monte_carlo.R

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("The data.table package is required for this audit.", call. = FALSE)
}

results_dir <- "results/primary_monte_carlo_lognormal_seed_20261122"
expected_commit <- "c29e92a3ea5b0116e46dbbdfb7465e60d96c88cf"
expected_methods <- c("complete_case", "substitute_zero", "substitute_loq_half",
  "substitute_loq", "ros_log", "mle_lognormal", "reverse_km")
moment_estimands <- c("mean", "median", "sd", "cv")
reverse_estimands <- c("mean", "q0.025", "q0.1", "q0.5", "q0.9", "q0.975")
key_cols <- c("scenario_id", "replicate_id", "method", "estimand")

checks <- list()
add_check <- function(name, status, details) {
  checks[[length(checks) + 1L]] <<- data.table::data.table(
    check = name, status = status, details = as.character(details))
}
check_status <- function(ok) if (isTRUE(ok)) "superato" else "anomalia"
close_enough <- function(x, y, tolerance = 1e-12) {
  length(x) == length(y) && all(is.finite(x) == is.finite(y)) &&
    all(is.na(x) == is.na(y)) &&
    (!any(is.finite(x)) || max(abs(x[is.finite(x)] - y[is.finite(y)])) <= tolerance)
}

if (!dir.exists(results_dir)) stop("Results directory is missing: ", results_dir, call. = FALSE)
expected_files <- c(paste0("n", rep(c(30, 100, 500), each = 2),
  "_LOQ_", rep(c(10, 40), times = 3), ".rds"), "results.rds", "run_manifest.rds")
present_files <- list.files(results_dir, pattern = "\\.rds$", full.names = FALSE)
add_check("Otto file RDS previsti", check_status(setequal(present_files, expected_files)),
  paste("presenti:", length(present_files), "attesi:", length(expected_files)))

git_head <- tryCatch(system2("git", "rev-parse HEAD", stdout = TRUE, stderr = FALSE)[1L],
  error = function(e) NA_character_)
add_check("Commit di riferimento", check_status(identical(git_head, expected_commit)),
  paste("HEAD =", git_head))

file_paths <- file.path(results_dir, expected_files)
if (!all(file.exists(file_paths))) {
  missing <- expected_files[!file.exists(file_paths)]
  add_check("Integrita file RDS", "anomalia", paste("mancanti:", paste(missing, collapse = ", ")))
  stop("Required RDS files are missing; see checks printed below.", call. = FALSE)
}
file_info <- data.table::data.table(
  file = expected_files,
  bytes = as.numeric(file.info(file_paths)$size),
  modified = format(file.info(file_paths)$mtime, "%Y-%m-%d %H:%M:%S %Z"),
  sha256 = NA_character_)
sha_output <- tryCatch(system2("sha256sum", file_paths, stdout = TRUE, stderr = FALSE),
  error = function(e) character())
if (length(sha_output) == length(file_paths)) {
  file_info[, sha256 := sub("  .*", "", sha_output)]
} else {
  add_check("Hash SHA-256", "non verificabile", "sha256sum non disponibile o output incompleto")
}

manifest <- readRDS(file.path(results_dir, "run_manifest.rds"))
combined <- readRDS(file.path(results_dir, "results.rds"))
if (!data.table::is.data.table(combined)) combined <- data.table::as.data.table(combined)
config <- manifest$signature$config
expected_scenarios <- data.table::as.data.table(expand.grid(
  sample_size = c(30L, 100L, 500L), loq_name = c("LOQ_10", "LOQ_40"),
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE))
expected_scenarios[, `:=`(
  scenario_index = seq_len(.N),
  scenario_id = paste0("n", sample_size, "_", loq_name),
  theoretical_censoring_probability = c(0.10, 0.40)[match(loq_name, c("LOQ_10", "LOQ_40"))])]
expected_scenarios[, loq := stats::qlnorm(theoretical_censoring_probability,
  meanlog = 0, sdlog = 1)]

add_check("Configurazione lognormale", check_status(
  identical(config$distribution, "lognormal") && identical(config$meanlog, 0) &&
    identical(config$sdlog, 1)), "lognormal; meanlog=0; sdlog=1")
add_check("Disegno campionario e LOQ", check_status(
  identical(as.integer(config$sample_sizes), c(30L, 100L, 500L)) &&
    identical(as.integer(manifest$signature$n_replicates), 2000L) &&
    identical(as.integer(config$n_monte_carlo), 2000L) &&
    close_enough(unname(config$loq_probabilities), c(0.10, 0.40)) &&
    close_enough(unname(config$loq_values), stats::qlnorm(c(0.10, 0.40), 0, 1))),
  "n=30/100/500; LOQ teorici 10%/40%; 2.000 repliche/condizione")
add_check("Seed", check_status(identical(as.integer(manifest$seed), 20261122L) &&
  identical(as.integer(manifest$signature$seed), 20261122L) &&
  grepl("20261122L", paste(readLines("simulations/run_primary_monte_carlo.R"), collapse = " ",
    useBytes = TRUE), fixed = TRUE)),
  paste("valore intero serializzato:", manifest$seed, "(notazione L nel launcher)"))
add_check("Metodi e inferenza", check_status(identical(manifest$methods, expected_methods) &&
  identical(manifest$signature$methods, expected_methods) &&
  identical(manifest$signature$inference, FALSE) &&
  all(combined$inference_enabled %in% FALSE)),
  paste("sette metodi previsti; inferenza disabilitata"))
add_check("Completamento manifest", check_status(isTRUE(manifest$complete) &&
  identical(as.integer(manifest$total_completed_replicates), 12000L) &&
  all(manifest$completed_by_scenario[expected_scenarios$scenario_id] == 2000L)),
  paste("complete=", manifest$complete, "; totale=", manifest$total_completed_replicates))
current_package_versions <- stats::setNames(vapply(names(manifest$package_versions), function(pkg)
  if (requireNamespace(pkg, quietly = TRUE)) as.character(utils::packageVersion(pkg)) else NA_character_,
  character(1)), names(manifest$package_versions))
versions_match <- identical(manifest$R_version, R.version.string) &&
  identical(unname(current_package_versions), unname(manifest$package_versions))
add_check("Versioni e firma", check_status(length(manifest$R_version) == 1L &&
  length(manifest$package_versions) >= 1L && versions_match && length(manifest$source_md5) >= 1L &&
  identical(manifest$source_md5, manifest$signature$source_md5) &&
  identical(manifest$R_version, manifest$signature$R_version)),
  paste("R:", manifest$R_version, "; pacchetti:",
    paste(names(manifest$package_versions), manifest$package_versions, collapse = ", "),
    "; versioni correnti coincidenti:", versions_match))
current_source_md5 <- if (all(file.exists(names(manifest$source_md5))))
  tools::md5sum(names(manifest$source_md5)) else NA
source_hash_ok <- !anyNA(current_source_md5) && identical(unname(current_source_md5),
  unname(manifest$source_md5))
add_check("Hash sorgenti", check_status(source_hash_ok),
  if (source_hash_ok) "hash sorgenti coincidenti con il manifest" else
    "file sorgente mancanti o hash differenti dal manifest")

scenario_tables <- list()
checkpoint_audit <- list()
for (i in seq_len(nrow(expected_scenarios))) {
  scenario <- expected_scenarios[i]
  scenario_id <- scenario$scenario_id[[1L]]
  checkpoint_path <- file.path(results_dir, paste0(scenario_id, ".rds"))
  checkpoint <- readRDS(checkpoint_path)
  dat <- data.table::as.data.table(checkpoint$results)
  scenario_tables[[scenario_id]] <- dat
  valid_ids <- identical(sort(unique(as.integer(dat$replicate_id))), 1:2000)
  unique_keys <- !anyDuplicated(dat, by = key_cols)
  metadata_by_rep <- dat[, .(
    rows = .N, n_censored_values = data.table::uniqueN(n_censored),
    censor_proportion_values = data.table::uniqueN(observed_censoring_proportion),
    loq_values = data.table::uniqueN(loq), size_values = data.table::uniqueN(sample_size),
    theoretical_probability_values = data.table::uniqueN(theoretical_censoring_probability),
    seed_values = data.table::uniqueN(replicate_seed)), by = replicate_id]
  expected_replicate_seed <- as.integer((20261122 + scenario$scenario_index * 1000003 +
    metadata_by_rep$replicate_id) %% (.Machine$integer.max - 1) + 1)
  metadata_consistent <- nrow(metadata_by_rep) == 2000L && all(metadata_by_rep$rows == 30L) &&
    all(metadata_by_rep$n_censored_values == 1L) &&
    all(metadata_by_rep$censor_proportion_values == 1L) &&
    all(metadata_by_rep$loq_values == 1L) && all(metadata_by_rep$size_values == 1L) &&
    all(metadata_by_rep$theoretical_probability_values == 1L) &&
    all(metadata_by_rep$seed_values == 1L) &&
    identical(sort(unique(as.integer(dat$replicate_seed))), sort(unique(expected_replicate_seed)))
  checkpoint_ok <- identical(checkpoint$scenario_id, scenario_id) &&
    identical(checkpoint$signature, manifest$signature) && nrow(dat) == 60000L &&
    valid_ids && unique_keys && metadata_consistent && all(dat$scenario_id == scenario_id) &&
    all(dat$scenario_index == scenario$scenario_index) &&
    all(dat$sample_size == scenario$sample_size) && all(dat$loq_name == scenario$loq_name) &&
    all(abs(dat$loq - scenario$loq) <= 1e-12) &&
    all(abs(dat$theoretical_censoring_probability - scenario$theoretical_censoring_probability) <= 1e-12)
  add_check(paste("Checkpoint", scenario_id), check_status(checkpoint_ok),
    paste("righe=", nrow(dat), "; repliche uniche=", data.table::uniqueN(dat$replicate_id),
      "; ID 1:2000=", valid_ids, "; chiavi uniche=", unique_keys,
      "; metadati/seed per replica coerenti=", metadata_consistent))
  checkpoint_audit[[scenario_id]] <- data.table::data.table(
    scenario_id = scenario_id, rows = nrow(dat), replicates = data.table::uniqueN(dat$replicate_id),
    min_replicate = if (nrow(dat)) min(dat$replicate_id) else NA_integer_,
    max_replicate = if (nrow(dat)) max(dat$replicate_id) else NA_integer_,
    unique_keys = unique_keys, signature_matches = identical(checkpoint$signature, manifest$signature))
}

scenario_summary <- data.table::rbindlist(lapply(expected_scenarios$scenario_id, function(id) {
  x <- scenario_tables[[id]]
  by_rep <- unique(x[, .(replicate_id, sample_size, loq, loq_name,
    theoretical_censoring_probability, n_censored, observed_censoring_proportion)])
  p <- by_rep$observed_censoring_proportion
  probs <- stats::quantile(p, c(0.025, 0.25, 0.75, 0.975), names = FALSE, type = 7)
  scenario <- expected_scenarios[scenario_id == id]
  data.table::data.table(scenario_id = id, sample_size = scenario$sample_size,
    loq_name = scenario$loq_name, loq = scenario$loq,
    theoretical_censoring_probability = scenario$theoretical_censoring_probability,
    mean = mean(p), median = stats::median(p), sd = stats::sd(p), min = min(p), max = max(p),
    q025 = probs[1L], q25 = probs[2L], q75 = probs[3L], q975 = probs[4L],
    distinct_proportions = data.table::uniqueN(p),
    censor_count_valid = all(by_rep$n_censored >= 0L & by_rep$n_censored <= by_rep$sample_size) &&
      all(abs(by_rep$observed_censoring_proportion - by_rep$n_censored / by_rep$sample_size) < 1e-12),
    fixed_loq = data.table::uniqueN(by_rep$loq) == 1L)
}), use.names = TRUE)
expected_rows <- data.table::rbindlist(lapply(expected_scenarios$scenario_id, function(id) {
  data.table::rbindlist(lapply(expected_methods, function(method) {
    estimands <- if (method == "reverse_km") reverse_estimands else moment_estimands
    data.table::data.table(scenario_id = id,
      replicate_id = rep(1:2000, each = length(estimands)), method = method,
      estimand = rep(estimands, times = 2000L))
  }))
}))
actual_keys <- unique(combined[, ..key_cols])
data.table::setkeyv(expected_rows, key_cols)
data.table::setkeyv(actual_keys, key_cols)
missing_keys <- data.table::fsetdiff(expected_rows, actual_keys)
extra_keys <- data.table::fsetdiff(actual_keys, expected_rows)
add_check("Completezza results.rds", check_status(nrow(combined) == 360000L &&
  data.table::uniqueN(combined$scenario_id) == 6L && nrow(missing_keys) == 0L &&
  nrow(extra_keys) == 0L), paste("righe=", nrow(combined), "; chiavi attese=", nrow(expected_rows),
  "; mancanti=", nrow(missing_keys), "; extra=", nrow(extra_keys)))
add_check("Unicita chiavi", check_status(!anyDuplicated(combined, by = key_cols)),
  paste("chiave scenario_id/replicate_id/method/estimand; duplicati=", anyDuplicated(combined, by = key_cols)))

checkpoint_combined <- data.table::rbindlist(scenario_tables, use.names = TRUE, fill = TRUE)
data.table::setorderv(checkpoint_combined, key_cols)
combined_sorted <- data.table::copy(combined)
data.table::setorderv(combined_sorted, key_cols)
same_payload <- nrow(checkpoint_combined) == nrow(combined_sorted) &&
  isTRUE(all.equal(as.data.frame(checkpoint_combined), as.data.frame(combined_sorted),
    check.attributes = FALSE, tolerance = 0))
add_check("Coerenza checkpoint/results.rds", check_status(same_payload),
  paste("confronto record per record su", nrow(combined), "righe"))

all_n <- all(combined$sample_size %in% c(30L, 100L, 500L))
row_counts <- combined[, .N, by = scenario_id]
add_check("Copertura scenario e repliche", check_status(all(row_counts$N == 60000L) && all_n),
  paste(paste(row_counts$scenario_id, row_counts$N, sep = "="), collapse = "; "))

expected_methods_per_estimand <- function(method, estimand) {
  if (method == "reverse_km") estimand %in% reverse_estimands else estimand %in% moment_estimands
}
method_coverage <- combined[, .(estimands = data.table::uniqueN(estimand), rows = .N), by = .(method)]
expected_method_counts <- c(complete_case = 48000L, substitute_zero = 48000L,
  substitute_loq_half = 48000L, substitute_loq = 48000L, ros_log = 48000L,
  mle_lognormal = 48000L, reverse_km = 72000L)
add_check("Copertura metodi e statistiche", check_status(
  setequal(unique(combined$method), expected_methods) &&
    all(method_coverage$rows[match(names(expected_method_counts), method_coverage$method)] ==
      expected_method_counts)),
  paste(paste(method_coverage$method, method_coverage$rows, sep = "="), collapse = "; "))

# Verify both error fields against their own, distinct reference columns.
expected_sample_error <- combined$estimate - combined$sample_truth
expected_population_error <- combined$estimate - combined$population_target
valid_sample_error <- is.finite(combined$estimate) & is.finite(combined$sample_truth)
valid_population_error <- is.finite(combined$estimate) & is.finite(combined$population_target)
sample_error_ok <- all(is.na(combined$error_vs_sample_truth[!valid_sample_error])) &&
  all(is.finite(combined$error_vs_sample_truth[valid_sample_error])) &&
  close_enough(combined$error_vs_sample_truth[valid_sample_error], expected_sample_error[valid_sample_error])
population_error_ok <- all(is.na(combined$error_vs_population[!valid_population_error])) &&
  all(is.finite(combined$error_vs_population[valid_population_error])) &&
  close_enough(combined$error_vs_population[valid_population_error], expected_population_error[valid_population_error])
add_check("Errori vs campione completo", check_status(sample_error_ok),
  "verificata identita estimate - sample_truth su tutti i record")
add_check("Errori vs parametri teorici", check_status(population_error_ok),
  "verificata identita estimate - population_target su tutti i record")
moment_rows <- combined[method != "reverse_km" & estimand %in% moment_estimands]
paired_reference_check <- moment_rows[, .(
  truth_values = data.table::uniqueN(sample_truth),
  population_values = data.table::uniqueN(population_target)),
  by = .(scenario_id, replicate_id, estimand)]
pair_refs_ok <- all(paired_reference_check$truth_values == 1L &
  paired_reference_check$population_values == 1L)
add_check("Benchmark appaiati tra metodi", check_status(pair_refs_ok),
  paste("un sample_truth e un target teorico per replica/statistica; gruppi=", nrow(paired_reference_check)))
theoretical_moments <- c(mean = exp(0 + 1 / 2), median = exp(0),
  sd = sqrt((exp(1) - 1) * exp(1)), cv = sqrt(exp(1) - 1))
moment_targets <- moment_rows[,
  .(target = population_target[1L]), by = estimand]
moment_target_ok <- all(vapply(seq_len(nrow(moment_targets)), function(i)
  isTRUE(all.equal(unname(moment_targets$target[i]),
    unname(theoretical_moments[moment_targets$estimand[i]]),
    tolerance = 1e-12)), logical(1)))
quantile_rows <- combined[method == "reverse_km" & estimand %in% reverse_estimands[-1L]]
quantile_targets_ok <- all(vapply(seq_len(nrow(quantile_rows)), function(i) {
  p <- as.numeric(sub("^q", "", quantile_rows$estimand[i]))
  isTRUE(all.equal(unname(quantile_rows$population_target[i]),
    unname(stats::qlnorm(p, 0, 1)), tolerance = 1e-12))
}, logical(1)))
add_check("Benchmark teorici", check_status(moment_target_ok && quantile_targets_ok),
  "momenti lognormali e quantili teorici coerenti con meanlog=0, sdlog=1")

finite_numeric_cols <- names(combined)[vapply(combined, is.numeric, logical(1))]
infinite_counts <- vapply(finite_numeric_cols, function(nm) sum(is.infinite(combined[[nm]])), integer(1))
nan_counts <- vapply(finite_numeric_cols, function(nm) sum(is.nan(combined[[nm]])), integer(1))
nonfinite_ok <- sum(infinite_counts) == 0L && sum(nan_counts) == 0L
add_check("Valori Inf/NaN", check_status(nonfinite_ok),
  paste("Inf=", sum(infinite_counts), "; NaN=", sum(nan_counts), "; NA sono valutati per stato"))

method_status <- combined[, .N, by = .(method, estimate_status, method_eval_status)]
estimand_status <- combined[, .N, by = .(method, estimand, estimand_status)]
combined[, estimate_availability := data.table::fifelse(is.finite(estimate), "estimate_available",
  data.table::fifelse(estimate_status %in% c("failed", "non_converged") |
    method_eval_status %in% c("failed", "non_converged"), "computational_failure",
  data.table::fifelse(estimand_status %in% c("non_identifiable", "theoretically_non_identifiable"),
    "structurally_non_identifiable",
  data.table::fifelse(estimand_status == "not_estimable_in_sample", "not_estimable_in_sample",
  data.table::fifelse(estimand_status == "not_returned", "not_returned_by_method",
  data.table::fifelse(estimate_status == "not_applicable", "not_applicable", "other_missing"))))))]
availability_status <- combined[, .N, by = .(method, estimand, estimate_availability)]
finite_but_nonidentified <- combined[is.finite(estimate) & estimand_status %in%
  c("non_identifiable", "theoretically_non_identifiable"),
  .N, by = .(scenario_id, method, estimand, estimand_status)]
missing_estimates <- combined[!is.finite(estimate)]
missing_status_consistent <- all(missing_estimates$estimate_availability != "other_missing")
add_check("Classificazione stime mancanti", check_status(missing_status_consistent),
  paste("mancanti=", nrow(missing_estimates), "; categorie=",
    paste(names(table(missing_estimates$estimate_availability)),
      as.integer(table(missing_estimates$estimate_availability)), collapse = ", ")))
error_rows <- combined[(!is.na(error_message) & nzchar(error_message)) |
  (!is.na(warnings) & nzchar(warnings))]
computational_failures <- combined[estimate_status %in% c("failed", "non_converged") |
  method_eval_status %in% c("failed", "non_converged")]
add_check("Errori/warning registrati", check_status(nrow(error_rows) == 0L),
  paste("record con messaggio o warning=", nrow(error_rows)))
add_check("Fallimenti computazionali", check_status(nrow(computational_failures) == 0L),
  paste("record failed/non_converged=", nrow(computational_failures)))
add_check("Inferenza disabilitata nei risultati", check_status(
  all(combined$inference_enabled %in% FALSE) && all(is.na(combined$interval_lower)) &&
  all(is.na(combined$interval_upper)) && all(is.na(combined$ucl95))),
  "nessun intervallo o UCL disponibile")

reverse <- combined[method == "reverse_km"]
reverse_expected <- list(
  mean_non_identifiable = 12000L,
  q0025_theoretically_non_identifiable = 12000L,
  q010_theoretically_non_identifiable = 12000L,
  median_not_estimable_in_sample = 415L,
  q0975_not_returned = 4000L)
reverse_observed <- c(
  mean_non_identifiable = reverse[estimand == "mean" & estimand_status == "non_identifiable", .N],
  q0025_theoretically_non_identifiable = reverse[estimand == "q0.025" & estimand_status == "theoretically_non_identifiable", .N],
  q010_theoretically_non_identifiable = reverse[estimand == "q0.1" & estimand_status == "theoretically_non_identifiable", .N],
  median_not_estimable_in_sample = reverse[estimand == "q0.5" & estimand_status == "not_estimable_in_sample", .N],
  q0975_not_returned = reverse[estimand == "q0.975" & estimand_status == "not_returned", .N])
add_check("Stati Reverse KM attesi", check_status(identical(as.integer(reverse_observed),
  as.integer(reverse_expected))), paste(paste(names(reverse_observed), reverse_observed, sep = "="), collapse = "; "))

summary_status <- data.table::rbindlist(checks, use.names = TRUE, fill = TRUE)
audit <- list(
  checks = summary_status,
  files = file_info,
  manifest = list(seed = manifest$seed, seed_scheme = manifest$seed_scheme,
    complete = manifest$complete, total_completed_replicates = manifest$total_completed_replicates,
    started_at_utc = manifest$started_at_utc, updated_at_utc = manifest$updated_at_utc,
    R_version = manifest$R_version, package_versions = manifest$package_versions,
    source_md5 = manifest$source_md5, git_commit = manifest$git_commit),
  checkpoints = data.table::rbindlist(checkpoint_audit),
  scenario_censoring = scenario_summary,
  method_status = method_status,
  estimand_status = estimand_status,
  estimate_availability = availability_status,
  finite_but_nonidentified = finite_but_nonidentified,
  finite_but_nonidentified = finite_but_nonidentified,
  reverse_km_expected_counts = data.table::data.table(
    category = names(reverse_expected), expected = as.integer(reverse_expected),
    observed = as.integer(reverse_observed)),
  errors_and_warnings = error_rows,
  computational_failures = computational_failures,
  row_count = nrow(combined),
  key_columns = key_cols,
  record_comparison_checkpoint_vs_results = same_payload
)

cat("\n=== Primary Monte Carlo data-quality audit ===\n")
print(audit$checks)
cat("\n=== File metadata and SHA-256 ===\n")
print(audit$files)
cat("\n=== Manifest ===\n")
print(audit$manifest)
cat("\n=== Checkpoints ===\n")
print(audit$checkpoints)
cat("\n=== Empirical censoring by condition ===\n")
print(audit$scenario_censoring)
cat("\n=== Method statuses ===\n")
print(audit$method_status)
cat("\n=== Estimand statuses ===\n")
print(audit$estimand_status)
cat("\n=== Estimate availability classifications ===\n")
print(audit$estimate_availability)
cat("\n=== Finite values with non-identifiable estimand status ===\n")
print(audit$finite_but_nonidentified)
cat("\n=== Reverse KM requested counts ===\n")
print(audit$reverse_km_expected_counts)
cat("\n=== Errors/warnings ===\n")
print(audit$errors_and_warnings)
cat("\n=== Computational failures ===\n")
print(audit$computational_failures)

if (any(summary_status$status == "anomalia")) {
  stop("Audit completed with anomalies; inspect the check table above.", call. = FALSE)
}
invisible(audit)
