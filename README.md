# Values below the limit of quantification (LOQ)

Reproducible R study of methods for left-censored environmental and chemical concentration data. The prepared primary Monte Carlo evaluates point estimators under a lognormal generating model. It has **not** been run.

## Definitive primary design

- Distribution: `Lognormal(meanlog = 0, sdlog = 1)`.
- Sample sizes: `n = 30`, `100`, and `500`.
- Fixed LOQs: the generating distribution's 10th and 40th percentiles. LOQs do not change between replicates or adapt to observed samples.
- Replicates: 2,000 in each of six sample-size/LOQ conditions, or 12,000 generated samples total.
- Approved master seed: `20261122L`.
- Comparison: paired. Every method receives the same censored sample; estimates are compared both with statistics from that replicate's complete latent sample and, separately, with the theoretical population parameter.
- Censoring: each latent sample is generated first, then the fixed LOQ is applied. The theoretical censoring probability and realized sample proportion are both retained; the latter varies across replicates.
- Methods: complete case, zero substitution, LOQ/2 substitution, LOQ substitution, log-ROS, censored lognormal MLE, and reverse Kaplan–Meier (KM).
- Scope: point-estimator assessment only. Inference, confidence intervals, UCLs, bootstrap, and coverage simulation are disabled for this run.

Reverse KM does not return a nonparametric mean under the common LOQ. Quantiles remain subject to identifiability and sample-estimability statuses. The methods do not all reconstruct a completed sample, so the output preserves each method's estimand and status.

## Scientific evidence

Functional checks establish that the code paths run; they do not establish statistical validity. See [estimator validation](reports/estimator_validation.md), [inference validation](reports/inference_validation.md), [bootstrap stability](reports/inference_bootstrap_stability.md), and [inference methods review](reports/inference_methods_review.md). Previous inferential validation results are retained unchanged and are not part of the primary point-estimator run.

## Runner and reproducibility

`R/12_primary_monte_carlo.R` generates one sample per replicate, evaluates the same sample with all seven methods, and stores paired errors against sample truth and population targets in separate fields. It supports atomic scenario checkpoints and signature-checked resumption. The manifest records the seed and seed derivation scheme, scenario configuration, replicate count, R and package versions, RNG kind, source hashes, progress, and completion state. An unrecognized nonempty output directory without a manifest is rejected.

The dedicated output directory is `results/primary_monte_carlo_lognormal_seed_20261122/`. It was checked to contain no prior results and reserved with a `.gitkeep` marker. The launcher pins the approved seed, 2,000 replicates, seven methods, and `inference = FALSE` in `simulations/run_primary_monte_carlo.R`.

The exact launch command, **after review**, is:

```sh
Rscript simulations/run_primary_monte_carlo.R
```

Do not run it as part of preparatory checks. The smoke test uses temporary files and verifies checkpoint/resume reproducibility without producing scientific outputs:

```sh
Rscript tests/manual/smoke_primary_monte_carlo.R
```

Other functional checks:

```sh
Rscript tests/manual/test_sample_generator.R
Rscript tests/manual/test_censoring_methods.R
Rscript tests/manual/test_single_replication.R
```

Do not run `validate_estimators.R`, `validate_inference_simulation.R`, or `validate_inference_stability.R` as part of the primary-run preparation; those perform simulations or bootstrap refits.
