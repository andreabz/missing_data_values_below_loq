# Values below the limit of quantification (LOQ)

Reproducible R study of methods for left-censored environmental and chemical concentration data. The repository now contains implemented method wrappers, validation reports, and a checkpointed runner prepared for the first lognormal Monte Carlo study. The definitive simulation has **not** been run; scientific review is required before launch.

## Study design

The first simulation uses `Lognormal(meanlog = 0, sdlog = 1)`, sample sizes 30, 100, and 500, and fixed LOQs at the generating distribution's 10th and 40th percentiles. The current configuration specifies 2,000 replicates per condition, for six conditions. Each sample is passed to all applicable methods, with theoretical censoring probability and realized sample censoring proportion retained separately.

The seven methods are complete case, zero substitution, LOQ/2 substitution, LOQ substitution, log-ROS, censored lognormal MLE, and reverse Kaplan–Meier (KM). No additional method is implemented for the primary design.

## Current scientific status

Method wrappers and deterministic functional checks exist. Independent estimator diagnostics and inferential coverage studies are reported separately; passing functional checks does not establish statistical validity.

- [Estimator validation](reports/estimator_validation.md) evaluates point-estimator behaviour in targeted diagnostics. Its results have reservations and do not validate confidence intervals.
- [Inference validation](reports/inference_validation.md) reports decisions by condition, method, estimand, and interval type. Some intervals are `not validated`; others are `inconclusive`.
- [Bootstrap stability supplement](reports/inference_bootstrap_stability.md) found material B=99 to B=499 endpoint changes in a targeted subset. Bootstrap procedures remain preliminary and must not be represented as validated.
- Reverse KM is used for identifiable quantiles, not a nonparametric mean under a common LOQ. Structural non-identifiability, sample non-estimability, and a quantile not returned by the fitted curve remain distinct.
- [Inference methods review](reports/inference_methods_review.md) records the method assumptions and unresolved choices.

The simulation must keep point-estimator availability, estimation accuracy, inferential coverage, interval availability, failures, and non-identifiability as separate outcomes.

## Implemented simulation architecture

`R/12_primary_monte_carlo.R` implements the primary simulation loop around `generate_left_censored_sample()`, `run_single_replication()`, and `evaluate_replication()`.

- A single observed sample is sent to all selected methods; latent values are used only for evaluation.
- Each output row identifies the scenario, replicate, method, estimand, estimate source, point-estimator status, inference status, errors against sample truth and population targets, warnings, and errors.
- Reverse-KM quantiles are stored as quantile estimands. Its mean is recorded as non-identifiable, with no mean estimate supplied.
- Results use `data.table` and are stored as RDS. Optional per-scenario checkpoints are written atomically; a run manifest records configuration, explicit master seed, derived replicate seed scheme, R and package versions, RNG kind, source-file hashes, progress, and completion state.
- Resuming checks the run signature, skips completed replicate IDs, and rejects duplicate replicate/method/estimand rows.

The runner requires an explicit seed; no definitive seed has been selected. Its current seed scheme derives a reproducible seed from the master seed, scenario index, and replicate ID. The runner defaults to the 2,000-replicate configuration but permits small smoke runs.

## Smoke check

Run the small, deterministic flow check from the repository root:

```sh
Rscript tests/manual/smoke_primary_monte_carlo.R
```

It runs one replicate in each of the six scenarios, interrupts and resumes from checkpoints, and compares the recovered output with an uninterrupted smoke run. It writes temporary artifacts and does not produce scientific results.

## Preparing a primary run

After scientific review, load the project files and pass the approved seed and a new output directory explicitly:

```r
source("R/00_config.R")
for (file in c("R/01_generate_sample.R", "R/02_method_results.R",
               "R/03_methods_simple.R", "R/04_method_ros.R", "R/05_method_mle.R",
               "R/06_method_reverse_km.R", "R/07_inference.R",
               "R/08_single_replication.R", "R/09_evaluation.R",
               "R/12_primary_monte_carlo.R")) source(file)

run <- run_primary_monte_carlo(
  seed = seed_approved_for_run,
  checkpoint_dir = run_directory
)
```

`seed_approved_for_run` and `run_directory` are placeholders to replace after review. The runner currently rejects `inference = TRUE`: `run_single_replication(inference = TRUE)` activates MLE and ROS bootstrap procedures together, although the bootstrap results are preliminary. A minimal follow-up is to add method- and procedure-specific inference selection, including a way to run the MLE log-Wald candidate without activating its bootstrap. Which procedures belong in the primary comparison remains a scientific decision.

The configuration's `n_bootstrap_pilot = 2000` is a pilot setting, not a final bootstrap choice. No bootstrap is run by the prepared Monte Carlo runner.

## Existing checks

Run functional checks from the repository root, for example:

```sh
Rscript tests/manual/test_sample_generator.R
Rscript tests/manual/test_censoring_methods.R
Rscript tests/manual/test_single_replication.R
Rscript tests/manual/validate_inference.R
```

`validate_estimators.R`, `validate_inference_simulation.R`, and `validate_inference_stability.R` perform simulations or bootstrap refits; they are validation activities, not smoke tests, and should not be run as part of a primary-run smoke check.

## Open decisions before definitive execution

- Approve the primary-run seed.
- Decide whether the planned 2,000 replicates suffice or whether the protocol's Monte Carlo precision review at 5,000 and 10,000 is needed.
- Select which inferential procedures to include, if any, given the per-condition status and bootstrap stability evidence in the reports.
- Set final bootstrap replication counts only after selecting a procedure that merits further inference.
- Design the later Gamma misspecification scenario separately.
- Decide whether any additional statistics belong in the study.

Do not describe `not validated` or `inconclusive` procedures as validated. Do not start the definitive Monte Carlo or Gamma scenario before the shared scientific review.
