# Small functional check only: no inferential bootstrap or scientific run.
root <- normalizePath(".", mustWork = TRUE)
for (file in c("R/00_config.R", "R/01_generate_sample.R", "R/02_method_results.R",
               "R/03_methods_simple.R", "R/04_method_ros.R", "R/05_method_mle.R",
               "R/06_method_reverse_km.R", "R/07_inference.R", "R/08_single_replication.R",
               "R/09_evaluation.R", "R/12_primary_monte_carlo.R")) source(file)

smoke_config <- simulation_config
smoke_config$n_monte_carlo <- 1L
smoke_seed <- 20261121L
temporary_root <- tempfile("primary-mc-smoke-")
dir.create(temporary_root)
checkpoint_dir <- file.path(temporary_root, "checkpoint")

first <- run_primary_monte_carlo(seed = smoke_seed, config = smoke_config,
  checkpoint_dir = checkpoint_dir, checkpoint_every = 1L,
  max_new_replicates = 3L)
stopifnot(!first$complete, sum(first$manifest$completed_by_scenario) == 3L,
  nrow(first$pending) == 3L)
signature_mismatch <- try(run_primary_monte_carlo(seed = smoke_seed + 1L,
  config = smoke_config, checkpoint_dir = checkpoint_dir), silent = TRUE)
stopifnot(inherits(signature_mismatch, "try-error"),
  grepl("does not match configuration", as.character(signature_mismatch), fixed = TRUE))

resumed <- run_primary_monte_carlo(seed = smoke_seed, config = smoke_config,
  checkpoint_dir = checkpoint_dir, checkpoint_every = 1L, resume = TRUE)
stopifnot(resumed$complete, nrow(resumed$pending) == 0L,
  length(unique(resumed$results$scenario_id)) == 6L,
  all(resumed$manifest$completed_by_scenario == 1L))

baseline <- run_primary_monte_carlo(seed = smoke_seed, config = smoke_config)
stopifnot(baseline$complete,
  isTRUE(all.equal(as.data.frame(resumed$results), as.data.frame(baseline$results),
    check.attributes = FALSE)))

stopifnot(!any(resumed$results$inference_enabled),
  all(is.na(resumed$results$interval_available)),
  all(is.na(resumed$results$interval_lower)), all(is.na(resumed$results$interval_upper)),
  all(is.na(resumed$results$ucl95)),
  !any(grepl("bootstrap", resumed$results$warnings, ignore.case = TRUE), na.rm = TRUE))
keys <- resumed$results[, .N, by = .(scenario_id, replicate_id, method, estimand)]
stopifnot(all(keys$N == 1L))
stopifnot(identical(primary_mc_replicate_seed(smoke_seed, 2L, 1L),
  primary_mc_replicate_seed(smoke_seed, 2L, 1L)),
  primary_mc_replicate_seed(smoke_seed, 2L, 1L) !=
    primary_mc_replicate_seed(smoke_seed, 2L, 2L))
unsafe_dir <- file.path(temporary_root, "unrecognized-files")
dir.create(unsafe_dir)
writeLines("preserve", file.path(unsafe_dir, "keep-me.txt"))
unsafe_error <- try(run_primary_monte_carlo(seed = smoke_seed, config = smoke_config,
  checkpoint_dir = unsafe_dir), silent = TRUE)
stopifnot(inherits(unsafe_error, "try-error"),
  grepl("not empty", as.character(unsafe_error), fixed = TRUE),
  file.exists(file.path(unsafe_dir, "keep-me.txt")))
inference_error <- try(run_primary_monte_carlo(seed = smoke_seed, config = smoke_config,
  inference = TRUE), silent = TRUE)
stopifnot(inherits(inference_error, "try-error"),
  grepl("requires inference=FALSE", as.character(inference_error), fixed = TRUE))

cat("Primary Monte Carlo smoke check passed (6 scenarios x 1 replicate; no bootstrap).\n")
unlink(temporary_root, recursive = TRUE)
