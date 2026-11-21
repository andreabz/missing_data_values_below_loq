# Definitive first lognormal point-estimator Monte Carlo.
# Run only after review of the prepared design: Rscript simulations/run_primary_monte_carlo.R
if (!file.exists("R/00_config.R"))
  stop("Run this script from the repository root")

for (file in c("R/00_config.R", "R/01_generate_sample.R", "R/02_method_results.R",
               "R/03_methods_simple.R", "R/04_method_ros.R", "R/05_method_mle.R",
               "R/06_method_reverse_km.R", "R/07_inference.R", "R/08_single_replication.R",
               "R/09_evaluation.R", "R/12_primary_monte_carlo.R")) source(file)

primary_run <- list(
  seed = 20261122L,
  n_replicates = 2000L,
  methods = primary_mc_methods,
  inference = FALSE,
  checkpoint_dir = "results/primary_monte_carlo_lognormal_seed_20261122"
)
stopifnot(identical(simulation_config$distribution, "lognormal"),
  identical(simulation_config$meanlog, 0), identical(simulation_config$sdlog, 1),
  identical(simulation_config$sample_sizes, c(30L, 100L, 500L)),
  identical(unname(simulation_config$loq_probabilities), c(0.10, 0.40)),
  identical(simulation_config$n_monte_carlo, 2000L),
  identical(primary_run$n_replicates, simulation_config$n_monte_carlo),
  length(primary_mc_scenarios(simulation_config)$scenario_id) == 6L,
  identical(primary_run$methods, primary_mc_methods),
  identical(primary_run$inference, FALSE))

if (dir.exists(primary_run$checkpoint_dir)) {
  entries <- setdiff(list.files(primary_run$checkpoint_dir, all.files = TRUE, no.. = TRUE), ".gitkeep")
  if (length(entries) && !"run_manifest.rds" %in% entries)
    stop("Results directory contains files but no run manifest; refusing to continue")
} else if (!dir.create(primary_run$checkpoint_dir, recursive = TRUE)) {
  stop("Could not create the dedicated results directory")
}

run <- run_primary_monte_carlo(
  seed = primary_run$seed,
  config = simulation_config,
  n_replicates = primary_run$n_replicates,
  methods = primary_run$methods,
  inference = primary_run$inference,
  checkpoint_dir = primary_run$checkpoint_dir,
  checkpoint_every = 100L,
  resume = TRUE
)

message("Completed: ", run$complete,
  "; total completed scenario replicates: ", run$manifest$total_completed_replicates)
