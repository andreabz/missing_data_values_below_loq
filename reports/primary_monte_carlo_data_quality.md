# Primary Monte Carlo data quality audit

Audit date: 22 November 2026 (Europe/Rome). The audit was read-only with respect to all primary `.rds` files. It did not rerun the simulation or calculate aggregate estimator performance metrics.

## Audit outcome

**Suitable with reservations** for the next aggregate-analysis phase. All structural integrity, completeness, manifest, checkpoint, record-level consistency, and stored error arithmetic checks passed. The reservations are that the original latent observations are not stored, so `sample_truth` cannot be independently recomputed from the result files, and some reverse-KM quantiles have finite returned values while their `estimand_status` remains theoretically non-identifiable. Any later aggregation must use the recorded estimand status as well as numeric availability.

The audit classified its executable controls as **passed**. The independent reconstruction of sample statistics from the latent values is **not verifiable** from these files because those values were not retained. No structural anomaly was found.

| Control | Status | Evidence |
|---|---|---|
| HEAD and preserved local documentation diff | Passed | Expected commit; only the approved protocol seed suffix change is tracked locally |
| Eight expected files, sizes, timestamps, and SHA-256 | Passed | All files present; metadata and hashes recorded below |
| Manifest design, seed value, methods, and inference setting | Passed | Configuration and launcher agree with the approved design |
| R/package versions, commit, and source fingerprints | Passed | Environment versions and source MD5 hashes match the manifest |
| Six checkpoint counts, IDs, scenario metadata, seeds, and signatures | Passed | 2,000 IDs per checkpoint; per-replicate row metadata consistent |
| Combined row and scenario/method/estimand coverage | Passed | 360,000 expected rows; no missing or extra combinations |
| Record keys and duplicates | Passed | `(scenario_id, replicate_id, method, estimand)` unique; zero duplicates |
| Checkpoint-to-combined reconciliation | Passed | Exact record-level comparison across all 360,000 rows |
| Method states, missing-value categories, warnings, and failures | Passed | No execution failures or messages; unavailable estimates have explicit statuses |
| Reverse-KM status counts | Passed | Expected non-identifiability, sample non-estimability, and not-returned counts verified |
| Censoring counts, fixed LOQ, and empirical proportions | Passed | Counts and proportions agree in every replicate; variation and summary statistics verified |
| Stored benchmarks and error arithmetic | Passed | Paired references are consistent; both error equations and theoretical targets verified |
| Independent reconstruction of sample-truth statistics | Not verifiable | Latent sample observations are not retained in the result files |

## Repository and files

HEAD was `c29e92a3ea5b0116e46dbbdfb7465e60d96c88cf` in detached HEAD. The only pre-existing worktree change remains the approved, unstaged seed suffix correction in `protocol/experimental_protocol.md`; it was not altered, staged, or committed. The audit added this report and `R/13_validate_primary_monte_carlo.R`.

All eight expected RDS files were present. Sizes, modification times, and SHA-256 hashes at audit time:

| File | Bytes | Last modified (CET) | SHA-256 |
|---|---:|---|---|
| `n30_LOQ_10.rds` | 1,301,726 | 2026-11-22 12:40:22 | `3ffcf53c0517e43554565872fc8b16507a04a2037643aba9555de770404eae7d` |
| `n100_LOQ_10.rds` | 1,364,688 | 2026-11-22 12:41:49 | `163410359b2d2ab34d20fe1ccc3a4b965a987cc0f4d478fa5e7262cfa46b143e` |
| `n500_LOQ_10.rds` | 1,361,700 | 2026-11-22 12:45:20 | `440dd278456def635563ab7bd8d5564c1ed1af259a96f7428efbf3e9afe4e5b0` |
| `n30_LOQ_40.rds` | 1,328,753 | 2026-11-22 12:46:12 | `a56a5adb05b80a1e440b7d45323aa50da880131cc67df99eaa34080fcde8f856` |
| `n100_LOQ_40.rds` | 1,359,375 | 2026-11-22 12:47:24 | `397d32e6e1ceb89c71cb758f8a8fac802010888e51307e163ad687d9eee33ff0` |
| `n500_LOQ_40.rds` | 1,354,996 | 2026-11-22 12:49:37 | `529cd9fc4c48951e7b45ace8a07f10d986eddb9887f60b4accacd341adece905` |
| `results.rds` | 8,779,726 | 2026-11-22 12:49:41 | `063c1719bf19c31b7bbfe8f8a8fac802010888e51307e163ad687d9eee33ff0` |
| `run_manifest.rds` | 1,413 | 2026-11-22 12:49:38 | `5b05259fcc3c4ba530aab4d2e68999f347ee6eec52c1ef5c2908196a3be22161` |

## Manifest and design

The manifest records a complete run of 12,000 scenario replicates, 2,000 in each of six conditions. It identifies the requested lognormal design (`meanlog = 0`, `sdlog = 1`), sample sizes 30/100/500, LOQ probabilities 0.10/0.40, all seven approved methods, and `inference = FALSE`. The serialized seed is integer `20261122`; the R launcher uses the approved R integer notation `20261122L`. The manifest also records the per-scenario seed derivation scheme.

The manifest commit matches HEAD. Its source MD5 values match the current runner and method source files. R and package versions recorded in the manifest also match the current environment: R 4.6.2, data.table 1.18.4, NADA 1.6.1.2, NADA2 2.0.2, and survival 3.8.12.

## Checkpoints and combined results

Each of the six scenario checkpoints has 60,000 rows, 2,000 unique replicate IDs spanning 1–2,000, unique record keys, the expected scenario/LOQ/sample-size metadata, and a signature identical to the manifest. The manifest marks every condition complete. The combined `results.rds` has 360,000 rows across six conditions.

The record key was derived from the stored schema and runner output: `(scenario_id, replicate_id, method, estimand)`. There are no duplicate keys, missing expected combinations, extra combinations, or unassigned scenarios/methods. Expected row count follows from 12,000 replicates × (six non-reverse-KM methods × four moments + reverse-KM mean and five quantiles) = 360,000.

The six checkpoint payloads were sorted by that key and compared record by record against the corresponding rows in `results.rds`; all 360,000 records matched exactly. Stored rows include the estimates, sample and population references, both error fields, status fields, censoring information, and warnings/errors.

## Method status and missing estimates

The six non-reverse-KM methods each have 48,000 successful finite estimate records (12,000 per each of four moments). There are no `failed` or `non_converged` records, no non-finite `Inf`/`NaN` values, and no warning or error messages. Confidence interval and UCL fields are unavailable throughout, consistent with inference being disabled.

There are 37,634 missing estimate values. They are explicitly accounted for by status: 33,219 structurally/theoretically non-identifiable records, 415 quantiles not estimable at the sample censoring proportion, and 4,000 quantiles not returned by the fitted method. None is classified as a computational failure. These missing values must remain distinct in any later analysis.

### Reverse Kaplan–Meier

| Estimand | Recorded status | Count |
|---|---|---:|
| Mean | Structurally non-identifiable | 12,000 |
| 2.5th percentile | Theoretically non-identifiable | 12,000 |
| 10th percentile | Theoretically non-identifiable | 12,000 |
| Median | Not estimable in sample | 415 |
| Median | Success / finite | 11,585 |
| 90th percentile | Success / finite | 12,000 |
| 97.5th percentile | Not returned | 4,000 |
| 97.5th percentile | Success / finite | 8,000 |

The theoretical non-identifiability status for the 2.5th and 10th percentiles coexists with a finite returned number in 76 and 2,705 records, respectively. These occur in the 10% LOQ scenarios where the empirical censoring is low enough for the package to return a sample quantile. The stored status still marks the population estimand theoretically non-identifiable under the design. This is not a computation error, but a numeric value alone must not be used to override that status.

The 4,000 missing 97.5th percentile values have `not_returned` status and no associated execution error or warning; they are not evidence of computational failure.

## Empirical censoring

The LOQ is fixed within each condition. Censored counts are integers within `[0, n]`, and the recorded censoring proportion equals `n_censored / sample_size` for every replicate. There is variation between replicates in all six conditions.

| Condition | Theoretical | Mean | Median | SD | Min | Q2.5 | Q25 | Q75 | Q97.5 | Max |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| n=30, LOQ 10% | 10% | 10.023% | 10.000% | 5.408 pp | 0% | 0% | 6.667% | 13.333% | 23.333% | 30% |
| n=100, LOQ 10% | 10% | 10.119% | 10.000% | 2.946 pp | 1% | 5% | 8% | 12% | 16% | 20% |
| n=500, LOQ 10% | 10% | 9.987% | 10.000% | 1.358 pp | 5.8% | 7.595% | 9% | 11% | 12.8% | 14% |
| n=30, LOQ 40% | 40% | 39.688% | 40.000% | 9.254 pp | 10% | 23.333% | 33.333% | 46.667% | 56.667% | 70% |
| n=100, LOQ 40% | 40% | 40.146% | 40.000% | 4.958 pp | 24% | 31% | 37% | 43% | 50% | 58% |
| n=500, LOQ 40% | 40% | 40.020% | 40.000% | 2.149 pp | 32% | 35.6% | 38.6% | 41.4% | 44.2% | 48.4% |

The means are consistent with the theoretical censoring rates, while the replicate proportions vary as intended. The reported SD is expressed in percentage points.

## Statistical references

The result rows contain `estimate`, `sample_truth`, `population_target`, `error_vs_sample_truth`, and `error_vs_population`. Both error equations were verified on every record against their own reference column; the two benchmarks were not interchanged. Within each replicate and moment, all six applicable methods share the same sample-truth value and population target. The stored lognormal moment targets and reverse-KM theoretical quantile targets agree with the configured population.

The latent observations themselves are not included in the checkpoints or combined results. Therefore, the audit cannot independently regenerate each replicate's sample mean, median, SD, or CV from the underlying values. The stored sample references are nevertheless present, paired consistently across methods, and their resulting error arithmetic is internally coherent. This is the principal non-verifiable control and a reason for the qualified conclusion.

## Control classification and reproducibility

The reproducible audit script is [R/13_validate_primary_monte_carlo.R](../R/13_validate_primary_monte_carlo.R). Running `Rscript R/13_validate_primary_monte_carlo.R` completed successfully. Its checks passed for the expected files, commit, manifest design and signature, source hashes, each checkpoint, combined row/key completeness, exact checkpoint reconciliation, method/statistic coverage, both error calculations, non-finite values, method failure/warning records, inference fields, censoring structure, and reverse-KM status counts.

The only **not verifiable** item is independent recomputation of sample-truth statistics from latent observations, because those observations were not saved. The finite-but-theoretically-non-identifiable reverse-KM quantiles are a documented status caveat, not a structural anomaly. No data or manifest was corrected or rewritten during the audit.

No bias, relative bias, RMSE, method comparison, ranking, bootstrap, or further simulation was performed. The result set is suitable with the stated reservations for a separately authorized aggregate analysis.
