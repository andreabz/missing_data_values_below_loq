# Reproducible, checkpointed primary Monte Carlo runner.
# Inference is deliberately disabled until method-specific procedures are
# selected; see the protocol and README for the current bootstrap limitation.

primary_mc_methods <- c("complete_case", "substitute_zero", "substitute_loq_half",
  "substitute_loq", "ros_log", "mle_lognormal", "reverse_km")

primary_mc_source_files <- c("R/00_config.R", "R/01_generate_sample.R",
  "R/02_method_results.R", "R/03_methods_simple.R", "R/04_method_ros.R",
  "R/05_method_mle.R", "R/06_method_reverse_km.R", "R/07_inference.R",
  "R/08_single_replication.R", "R/09_evaluation.R", "R/12_primary_monte_carlo.R",
  "simulations/run_primary_monte_carlo.R")

validate_primary_mc_config <- function(config) {
  required <- c("distribution", "meanlog", "sdlog", "sample_sizes",
    "loq_probabilities", "loq_values", "n_monte_carlo", "population_parameters")
  if (!is.list(config) || !all(required %in% names(config)) ||
      !identical(config$distribution, "lognormal"))
    stop("config must be a lognormal simulation_config with all required fields")
  if (!is.numeric(config$meanlog) || length(config$meanlog) != 1L ||
      !is.finite(config$meanlog) || !is.numeric(config$sdlog) ||
      length(config$sdlog) != 1L || !is.finite(config$sdlog) || config$sdlog <= 0)
    stop("config meanlog and sdlog are invalid")
  if (!is.numeric(config$sample_sizes) || !length(config$sample_sizes) ||
      any(!is.finite(config$sample_sizes)) || any(config$sample_sizes < 2) ||
      any(config$sample_sizes != floor(config$sample_sizes)))
    stop("config sample_sizes must contain integer sample sizes >= 2")
  loq_names <- names(config$loq_values)
  if (!is.numeric(config$loq_probabilities) || is.null(names(config$loq_probabilities)) ||
      !length(loq_names) || !identical(names(config$loq_probabilities), loq_names) ||
      any(!is.finite(config$loq_probabilities)) ||
      any(config$loq_probabilities <= 0 | config$loq_probabilities >= 1) ||
      !is.numeric(config$loq_values) || any(!is.finite(config$loq_values)) ||
      any(config$loq_values <= 0))
    stop("config loq_probabilities and loq_values must be matching named vectors")
  if (length(config$n_monte_carlo) != 1L || !is.finite(config$n_monte_carlo) ||
      config$n_monte_carlo < 1 || config$n_monte_carlo != floor(config$n_monte_carlo))
    stop("config n_monte_carlo must be a positive integer")
  invisible(TRUE)
}

primary_mc_package_versions <- function() {
  packages <- c("data.table", "NADA", "NADA2", "survival")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing required packages: ", paste(missing, collapse = ", "))
  stats::setNames(vapply(packages, function(pkg)
    as.character(utils::packageVersion(pkg)), character(1)), packages)
}

primary_mc_signature <- function(config, n_replicates, seed, methods, inference) {
  absent <- primary_mc_source_files[!file.exists(primary_mc_source_files)]
  if (length(absent)) stop("Cannot fingerprint missing source files: ", paste(absent, collapse = ", "))
  list(config = config, n_replicates = as.integer(n_replicates), seed = as.integer(seed),
    seed_scheme = "set.seed((master_seed + scenario_index*1000003 + replicate_id) mod (.Machine$integer.max - 1) + 1)",
    methods = methods, inference = inference, package_versions = primary_mc_package_versions(),
    R_version = R.version.string, RNGkind = RNGkind(),
    source_md5 = tools::md5sum(primary_mc_source_files))
}

primary_mc_scenarios <- function(config) {
  grid <- expand.grid(sample_size = config$sample_sizes,
    loq_name = names(config$loq_values), KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE)
  data.table::as.data.table(grid)[, `:=`(
    scenario_index = seq_len(.N),
    scenario_id = paste0("n", sample_size, "_", loq_name),
    loq = unname(config$loq_values[loq_name]),
    theoretical_censoring_probability = unname(config$loq_probabilities[loq_name]))]
}

primary_mc_replicate_seed <- function(master_seed, scenario_index, replicate_id) {
  modulus <- .Machine$integer.max - 1
  as.integer((as.double(master_seed) + as.double(scenario_index) * 1000003 +
    as.double(replicate_id)) %% modulus) + 1L
}

primary_mc_safe_error <- function(estimate, target) {
  if (length(estimate) != 1L || length(target) != 1L ||
      !is.finite(estimate) || !is.finite(target)) return(NA_real_)
  estimate - target
}

primary_mc_method_rows <- function(replication, evaluation, sample, scenario,
                                   replicate_id, replicate_seed, config, inference_enabled) {
  population <- config$population_parameters
  sample_stats <- evaluation$sample_truth
  source_by_method <- c(complete_case = "observed_uncensored_subset",
    substitute_zero = "observed_value_substitution",
    substitute_loq_half = "observed_value_substitution",
    substitute_loq = "observed_value_substitution",
    ros_log = "ROS_reconstructed_distribution",
    mle_lognormal = "censored_lognormal_model",
    reverse_km = "reverse_KM_empirical_distribution")
  rows <- list()
  for (method in names(replication$methods)) {
    result <- replication$methods[[method]]
    method_rows <- list()
    if (method == "reverse_km") {
      method_rows[[length(method_rows) + 1L]] <- list(estimand = "mean", estimand_type = "moment",
        probability = NA_real_, estimate = NA_real_, sample_truth = NA_real_,
        population_target = NA_real_, estimand_status = "non_identifiable")
      quantiles <- result$estimates$quantiles
      if (length(quantiles)) for (qname in names(quantiles)) {
        probability <- as.numeric(sub("^q", "", qname))
        estimate <- unname(quantiles[[qname]])
        q_status <- if (probability <= scenario$theoretical_censoring_probability) {
          "theoretically_non_identifiable"
        } else if (probability <= sample$metadata$observed_censoring_proportion) {
          "not_estimable_in_sample"
        } else if (is.finite(estimate)) "success" else "not_returned"
        sample_q <- unname(stats::quantile(sample$truth, probability, names = FALSE, type = 7))
        population_q <- stats::qlnorm(probability, meanlog = config$meanlog, sdlog = config$sdlog)
        method_rows[[length(method_rows) + 1L]] <- list(estimand = qname, estimand_type = "quantile",
          probability = probability, estimate = estimate, sample_truth = sample_q,
          population_target = population_q, estimand_status = q_status)
      }
    } else {
      for (statistic in c("mean", "median", "sd", "cv")) {
        estimate <- result$estimates[[statistic]]
        if (is.null(estimate) || length(estimate) != 1L) estimate <- NA_real_
        sample_target <- sample_stats[[statistic]]
        population_target <- population[[statistic]]
        method_rows[[length(method_rows) + 1L]] <- list(estimand = statistic, estimand_type = "moment",
          probability = NA_real_, estimate = estimate, sample_truth = sample_target,
          population_target = population_target,
          estimand_status = if (is.finite(estimate)) result$estimate_status else
            if (result$estimate_status == "success") "not_returned" else result$estimate_status)
      }
    }

    method_eval <- evaluation$methods[[method]]
    completed_method_rows <- lapply(method_rows, function(row) {
      ci <- result$inference$ci_mean_95
      ci_available <- if (!inference_enabled) NA else
        length(ci) >= 2L && all(is.finite(ci)) && ci[[1L]] <= ci[[2L]]
      ucl <- result$inference$ucl95
      ucl_available <- if (!inference_enabled) NA else is.numeric(ucl) &&
        length(ucl) == 1L && is.finite(ucl)
      q_interval <- NULL
      q_intervals <- result$inference$quantile_intervals
      if (is.list(q_intervals) && row$estimand %in% names(q_intervals))
        q_interval <- q_intervals[[row$estimand]]
      lower <- if (identical(row$estimand, "mean") && length(ci) >= 2L) ci[[1L]] else
        if (!is.null(q_interval)) q_interval$ci[[1L]] else NA_real_
      upper <- if (identical(row$estimand, "mean") && length(ci) >= 2L) ci[[2L]] else
        if (!is.null(q_interval)) q_interval$ci[[2L]] else NA_real_
      c(row, list(
        scenario_id = scenario$scenario_id, scenario_index = scenario$scenario_index,
        sample_size = scenario$sample_size, loq_name = scenario$loq_name,
        loq = scenario$loq,
        theoretical_censoring_probability = scenario$theoretical_censoring_probability,
        replicate_id = as.integer(replicate_id), replicate_seed = replicate_seed,
        n_censored = sample$metadata$n_censored,
        observed_censoring_proportion = sample$metadata$observed_censoring_proportion,
        method = method, estimate_source = unname(source_by_method[[method]]),
        error_vs_sample_truth = primary_mc_safe_error(row$estimate, row$sample_truth),
        error_vs_population = primary_mc_safe_error(row$estimate, row$population_target),
        estimate_status = result$estimate_status,
        inference_enabled = inference_enabled, inference_status = result$inference$status,
        interval_lower = lower, interval_upper = upper, interval_available = ci_available,
        ucl95 = if (inference_enabled && identical(row$estimand, "mean")) ucl else NA_real_,
        ucl_available = if (inference_enabled && identical(row$estimand, "mean")) ucl_available else NA,
        warnings = if (length(result$warnings)) paste(result$warnings, collapse = " | ") else NA_character_,
        error_message = result$error, method_eval_status = method_eval$estimate_status))
    })
    rows <- c(rows, completed_method_rows)
  }
  data.table::rbindlist(rows, use.names = TRUE, fill = TRUE)
}

primary_mc_atomic_save <- function(object, path) {
  directory <- dirname(path)
  if (!dir.exists(directory) && !dir.create(directory, recursive = TRUE))
    stop("Cannot create checkpoint directory: ", directory)
  temporary <- tempfile(pattern = paste0(basename(path), "."), tmpdir = directory)
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  saveRDS(object, temporary)
  if (!file.rename(temporary, path)) stop("Atomic checkpoint rename failed for: ", path)
  invisible(path)
}

run_primary_monte_carlo <- function(seed, config = simulation_config,
                                    n_replicates = config$n_monte_carlo,
                                    methods = primary_mc_methods,
                                    inference = FALSE, checkpoint_dir = NULL,
                                    checkpoint_every = 100L, resume = TRUE,
                                    max_new_replicates = Inf) {
  validate_primary_mc_config(config)
  if (length(seed) != 1L || !is.numeric(seed) || !is.finite(seed) ||
      seed != floor(seed) || seed < 0 || seed > .Machine$integer.max)
    stop("seed must be one explicit integer in [0, .Machine$integer.max]")
  if (length(n_replicates) != 1L || !is.numeric(n_replicates) ||
      !is.finite(n_replicates) || n_replicates < 1 || n_replicates != floor(n_replicates))
    stop("n_replicates must be a positive integer")
  if (!is.character(methods) || !length(methods) || anyDuplicated(methods) ||
      any(!methods %in% primary_mc_methods)) stop("methods contains invalid or duplicate methods")
  if (!is.logical(inference) || length(inference) != 1L || is.na(inference))
    stop("inference must be TRUE or FALSE")
  if (!is.logical(resume) || length(resume) != 1L || is.na(resume))
    stop("resume must be TRUE or FALSE")
  if (isTRUE(inference))
    stop(paste0("Primary simulation currently requires inference=FALSE: run_single_replication(inference=TRUE) ",
      "activates MLE and ROS bootstrap together. Select procedures and method-specific bootstrap controls ",
      "after scientific review before enabling inference."))
  if (length(checkpoint_every) != 1L || !is.finite(checkpoint_every) ||
      checkpoint_every < 1 || checkpoint_every != floor(checkpoint_every))
    stop("checkpoint_every must be a positive integer")
  if (!(length(max_new_replicates) == 1L && is.numeric(max_new_replicates) &&
      ((is.finite(max_new_replicates) && max_new_replicates >= 0 &&
        max_new_replicates == floor(max_new_replicates)) ||
       (is.infinite(max_new_replicates) && max_new_replicates > 0))))
    stop("max_new_replicates must be a non-negative integer or Inf")
  if (!is.null(checkpoint_dir) && (!is.character(checkpoint_dir) ||
      length(checkpoint_dir) != 1L || is.na(checkpoint_dir) || !nzchar(checkpoint_dir)))
    stop("checkpoint_dir must be NULL or one non-empty path")

  old_rng_kind <- RNGkind()
  had_rng_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_rng_seed) old_rng_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    do.call(RNGkind, as.list(old_rng_kind))
    if (had_rng_seed) assign(".Random.seed", old_rng_seed, envir = .GlobalEnv) else
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)

  seed <- as.integer(seed)
  n_replicates <- as.integer(n_replicates)
  checkpoint_every <- as.integer(checkpoint_every)
  scenarios <- primary_mc_scenarios(config)
  signature <- primary_mc_signature(config, n_replicates, seed, methods, inference)
  started_at <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  run_manifest_path <- if (is.null(checkpoint_dir)) NULL else
    file.path(checkpoint_dir, "run_manifest.rds")
  if (!is.null(checkpoint_dir) && !dir.exists(checkpoint_dir) &&
      !dir.create(checkpoint_dir, recursive = TRUE)) stop("Cannot create checkpoint_dir")

  run_manifest <- list(signature = signature, config = config, methods = methods,
    seed = seed, seed_scheme = signature$seed_scheme,
    package_versions = signature$package_versions, R_version = signature$R_version,
    RNGkind = signature$RNGkind, source_md5 = signature$source_md5,
    git_commit = tryCatch(system2("git", c("rev-parse", "HEAD"), stdout = TRUE,
      stderr = FALSE)[[1L]], error = function(e) NA_character_),
    started_at_utc = started_at, updated_at_utc = started_at,
    completed_by_scenario = stats::setNames(integer(nrow(scenarios)), scenarios$scenario_id),
    complete = FALSE)

  if (!is.null(run_manifest_path) && file.exists(run_manifest_path)) {
    if (!isTRUE(resume)) stop("Checkpoint manifest exists; set resume=TRUE or choose a new directory")
    previous <- readRDS(run_manifest_path)
    if (!identical(previous$signature, signature))
      stop("Checkpoint manifest does not match configuration, seed, methods, packages, or source code")
    run_manifest$started_at_utc <- previous$started_at_utc
    run_manifest$git_commit <- previous$git_commit
    run_manifest$completed_by_scenario <- previous$completed_by_scenario
  } else if (!is.null(run_manifest_path)) {
    # Without a manifest, only an empty directory or its reservation marker is safe.
    existing <- list.files(checkpoint_dir, all.files = TRUE, no.. = TRUE)
    existing <- setdiff(existing, ".gitkeep")
    if (length(existing))
      stop("Checkpoint directory is not empty and has no resumable run manifest")
  }

  if (!is.null(run_manifest_path)) primary_mc_atomic_save(run_manifest, run_manifest_path)
  started_call <- proc.time()[["elapsed"]]
  new_replicates <- 0L
  scenario_results <- vector("list", nrow(scenarios))
  names(scenario_results) <- scenarios$scenario_id
  stopped_early <- FALSE

  for (scenario_row in seq_len(nrow(scenarios))) {
    scenario <- scenarios[scenario_row]
    scenario_id <- scenario$scenario_id[[1L]]
    checkpoint_path <- if (is.null(checkpoint_dir)) NULL else
      file.path(checkpoint_dir, paste0(scenario_id, ".rds"))
    results <- data.table::data.table()
    if (!is.null(checkpoint_path) && file.exists(checkpoint_path)) {
      if (!isTRUE(resume)) stop("Scenario checkpoint exists; set resume=TRUE or choose a new directory")
      saved <- readRDS(checkpoint_path)
      if (!identical(saved$signature, signature) || !identical(saved$scenario_id, scenario_id))
        stop("Scenario checkpoint signature mismatch: ", scenario_id)
      results <- saved$results
    }
    if (nrow(results)) {
      duplicate_keys <- results[, .N, by = .(replicate_id, method, estimand)][N > 1L]
      if (nrow(duplicate_keys)) stop("Checkpoint contains duplicate replicate/method/estimand rows")
      completed <- sort(unique(results$replicate_id))
      if (any(completed < 1L | completed > n_replicates))
        stop("Checkpoint contains replicate IDs outside requested range")
      methods_per_id <- results[, data.table::uniqueN(method), by = replicate_id][[2L]]
      if (any(methods_per_id != length(methods)))
        stop("Checkpoint contains incomplete method sets for a replicate")
    } else completed <- integer()

    pending <- setdiff(seq_len(n_replicates), completed)
    pending_new <- 0L
    for (replicate_id in pending) {
      if (new_replicates >= max_new_replicates) {
        stopped_early <- TRUE
        break
      }
      replicate_seed <- primary_mc_replicate_seed(seed, scenario$scenario_index[[1L]], replicate_id)
      set.seed(replicate_seed)
      sample <- generate_left_censored_sample(n = scenario$sample_size[[1L]],
        loq = scenario$loq[[1L]], meanlog = config$meanlog, sdlog = config$sdlog)
      replication <- run_single_replication(sample, replicate_id = replicate_id,
        methods = methods, inference = FALSE)
      evaluation <- evaluate_replication(replication, sample$truth,
        population_parameters = config$population_parameters)
      rows <- primary_mc_method_rows(replication, evaluation, sample, scenario,
        replicate_id, replicate_seed, config, inference_enabled = FALSE)
      results <- data.table::rbindlist(list(results, rows), use.names = TRUE, fill = TRUE)
      completed <- c(completed, replicate_id)
      pending_new <- pending_new + 1L
      new_replicates <- new_replicates + 1L
      if (!is.null(checkpoint_path) && pending_new %% checkpoint_every == 0L) {
        primary_mc_atomic_save(list(signature = signature, scenario_id = scenario_id,
          results = results, saved_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)),
          checkpoint_path)
      }
      if (new_replicates >= max_new_replicates) {
        stopped_early <- TRUE
        break
      }
    }
    if (!is.null(checkpoint_path) && (pending_new > 0L || !file.exists(checkpoint_path)))
      primary_mc_atomic_save(list(signature = signature, scenario_id = scenario_id,
        results = results, saved_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)),
        checkpoint_path)
    scenario_results[[scenario_id]] <- results
    run_manifest$completed_by_scenario[[scenario_id]] <- length(unique(completed))
    run_manifest$updated_at_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
    if (!is.null(run_manifest_path)) primary_mc_atomic_save(run_manifest, run_manifest_path)
    if (stopped_early) break
  }

  if (!is.null(checkpoint_dir)) {
    for (scenario_row in seq_len(nrow(scenarios))) {
      scenario_id <- scenarios$scenario_id[[scenario_row]]
      if (!is.null(scenario_results[[scenario_id]])) next
      checkpoint_path <- file.path(checkpoint_dir, paste0(scenario_id, ".rds"))
      if (file.exists(checkpoint_path)) scenario_results[[scenario_id]] <- readRDS(checkpoint_path)$results
    }
  }
  results <- data.table::rbindlist(scenario_results, use.names = TRUE, fill = TRUE)
  if (is.null(results)) results <- data.table::data.table()
  if (nrow(results)) data.table::setorderv(results,
    c("scenario_index", "replicate_id", "method", "estimand"))
  completed_by_scenario <- vapply(scenarios$scenario_id, function(id) {
    x <- results[scenario_id == id]
    if (nrow(x)) data.table::uniqueN(x$replicate_id) else 0L
  }, integer(1))
  pending_rows <- lapply(seq_len(nrow(scenarios)), function(i) {
    id <- scenarios$scenario_id[[i]]
    done <- if (nrow(results)) unique(results[scenario_id == id, replicate_id]) else integer()
    todo <- setdiff(seq_len(n_replicates), done)
    if (!length(todo)) return(NULL)
    data.table::data.table(scenario_id = id, replicate_id = todo)
  })
  pending <- data.table::rbindlist(pending_rows, use.names = TRUE, fill = TRUE)
  if (is.null(pending)) pending <- data.table::data.table(scenario_id = character(),
    replicate_id = integer())
  complete <- nrow(pending) == 0L
  run_manifest$completed_by_scenario <- completed_by_scenario
  run_manifest$complete <- complete
  run_manifest$updated_at_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  run_manifest$elapsed_seconds_this_call <- unname(proc.time()[["elapsed"]] - started_call)
  run_manifest$total_completed_replicates <- sum(completed_by_scenario)
  if (!is.null(run_manifest_path)) {
    primary_mc_atomic_save(run_manifest, run_manifest_path)
    if (complete) primary_mc_atomic_save(results, file.path(checkpoint_dir, "results.rds"))
    else primary_mc_atomic_save(results, file.path(checkpoint_dir, "partial_results.rds"))
  }
  structure(list(results = results, manifest = run_manifest, pending = pending,
    complete = complete), class = "primary_monte_carlo_run")
}
