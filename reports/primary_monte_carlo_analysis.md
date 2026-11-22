# Primary Monte Carlo statistical analysis

## Objective and scope

This report analyzes the completed primary Monte Carlo without running new simulations, bootstrap, or inference. It separates error relative to the paired complete-sample statistic from error relative to the theoretical population parameter. All summaries remain stratified by the six scenarios, method, and estimand; no overall method ranking is produced.

The underlying data-quality audit judged the results **suitable with reservations**. This analysis applies the audit's Reverse-KM status caveat: a finite value is included in accuracy metrics only when its row-specific `estimand_status` is `success`.

## Design and inputs

The analysis used commit `c29e92a3ea5b0116e46dbbdfb7465e60d96c88cf` and `360000` rows from `results/primary_monte_carlo_lognormal_seed_20261122`. The manifest confirms a lognormal population (`meanlog = 0`, `sdlog = 1`), n = 30/100/500, fixed theoretical LOQs at 10% and 40%, 2,000 replicates per condition, seven methods, and inference disabled. The six checkpoint files were re-read and compared record by record with `results.rds`; the manifest and scenario signatures matched.

Input SHA-256 hashes were recorded before analysis and recalculated after report/figure generation; all eight primary RDS hashes remained unchanged during that analysis run. A later forensic check found that the current `results.rds` hash differs from the digest transcribed in the data-quality audit. Its current logical records match the six checkpoints exactly after key ordering; the historical byte provenance remains unverified. See [the hash investigation](primary_monte_carlo_hash_investigation.md) for the evidence and limits. This does not change the analysis-versus-checkpoint content reconciliation.

## Definitions and eligibility

For every scenario × method × estimand, the expected denominator is 2,000. `n_finite` counts finite numeric values regardless of status; `n_usable` counts finite estimates whose row-specific estimand status is `success` (and whose method status is not a failure). Availability is `n_usable / 2,000`. Missing values, finite-but-ineligible estimates, non-identifiability, sample non-estimability, not-returned values, and computational failures are reported separately.

Relative error is `(estimate - reference) / reference`. References with absolute value at most `sqrt(.Machine$double.eps)` are treated as zero/near-zero and excluded rather than replaced. The tables report relative-error denominators and each exclusion reason. No such exclusion is silently converted to zero.

For sample-truth comparisons, the stored `error_vs_sample_truth` is used. For population comparisons, the separately stored `error_vs_population` is used. Bias is mean population error; relative bias is the mean population-relative error; RMSE is `sqrt(mean(error^2))`. Error quantiles use type 7 sample quantiles. Estimator SD is computed across usable estimates. Sample-truth distributions are summarized separately from method errors.

Complete case and substitution methods summarize detected-only or constructed samples. ROS reconstructs censored values using log-ROS: the method relates normal quantiles to detected log-concentrations and predicts below-LOQ values from that regression. Its lower-tail extrapolation and reconstructed summaries depend on the adequacy of the log-scale distributional relationship; ROS is not assumption-free or equivalent to multiple imputation. In this simulation, the lognormal generator is compatible with that working relationship, so the results do not establish robustness under misspecification.

The lognormal MLE fits the censored lognormal likelihood: detected observations contribute their density and censored observations contribute the probability of falling below the LOQ. It estimates distribution parameters and derives moments from the fitted distribution; it does not reconstruct individual censored values. Because the generating distribution is lognormal here, the MLE specification is favorable by design. Its performance should not be generalized automatically to other distributions.

For Reverse KM, `mean` is structurally non-identifiable and is not assigned an accuracy metric. Quantiles marked `theoretically_non_identifiable`, `not_estimable_in_sample`, or `not_returned` are excluded even when a numeric value exists. Not-returned values are not necessarily computational failures. Reverse-KM accuracy metrics are conditional on a finite estimate with row status `success`; read them together with their usable counts and availability, because the eligible replicates differ by quantile and scenario. They are not directly comparable as unconditional performance summaries with methods returning estimates in all replicates. The common-LOQ results describe only quantiles that are identifiable and returned under these data conditions; no nonparametric mean is identified without additional assumptions. Reverse-KM quantile results are not pooled or ranked against moment-estimating methods.

## Data availability and reliability

Counts use 2,000 expected replicates for every row. `finite` is numeric availability; `usable` is the analysis-eligible count; `finite/status-excluded` records finite outputs excluded by their estimand status. Relative-exclusion counts are provided separately for each reference.

### n=30, theoretical censoring=10%

| Method | Estimand | Expected | Finite | Usable | Missing/nonfinite | NA | NaN | Inf | Finite, status-excluded | Theoretical non-ID | Structural non-ID | Not estimable in sample | Not returned | Computational failure | Available (%) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| complete_case | cv | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| complete_case | mean | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| complete_case | median | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| complete_case | sd | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| mle_lognormal | cv | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| mle_lognormal | mean | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| mle_lognormal | median | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| mle_lognormal | sd | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| reverse_km | mean | 2000 | 0 | 0 | 2000 | 2000 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0.00 |
| reverse_km | q0.025 | 2000 | 72 | 0 | 1928 | 1928 | 0 | 0 | 72 | 2000 | 0 | 0 | 0 | 0 | 0.00 |
| reverse_km | q0.1 | 2000 | 831 | 0 | 1169 | 1169 | 0 | 0 | 831 | 2000 | 0 | 0 | 0 | 0 | 0.00 |
| reverse_km | q0.5 | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| reverse_km | q0.9 | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| reverse_km | q0.975 | 2000 | 0 | 0 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0.00 |
| ros_log | cv | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| ros_log | mean | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| ros_log | median | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| ros_log | sd | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_loq | cv | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_loq | mean | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_loq | median | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_loq | sd | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_loq_half | cv | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_loq_half | mean | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_loq_half | median | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_loq_half | sd | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_zero | cv | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_zero | mean | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_zero | median | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_zero | sd | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |

### n=100, theoretical censoring=10%

| Method | Estimand | Expected | Finite | Usable | Missing/nonfinite | NA | NaN | Inf | Finite, status-excluded | Theoretical non-ID | Structural non-ID | Not estimable in sample | Not returned | Computational failure | Available (%) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| complete_case | cv | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| complete_case | mean | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| complete_case | median | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| complete_case | sd | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| mle_lognormal | cv | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| mle_lognormal | mean | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| mle_lognormal | median | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| mle_lognormal | sd | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| reverse_km | mean | 2000 | 0 | 0 | 2000 | 2000 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0.00 |
| reverse_km | q0.025 | 2000 | 4 | 0 | 1996 | 1996 | 0 | 0 | 4 | 2000 | 0 | 0 | 0 | 0 | 0.00 |
| reverse_km | q0.1 | 2000 | 901 | 0 | 1099 | 1099 | 0 | 0 | 901 | 2000 | 0 | 0 | 0 | 0 | 0.00 |
| reverse_km | q0.5 | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| reverse_km | q0.9 | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| reverse_km | q0.975 | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| ros_log | cv | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| ros_log | mean | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| ros_log | median | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| ros_log | sd | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_loq | cv | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_loq | mean | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_loq | median | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_loq | sd | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_loq_half | cv | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_loq_half | mean | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_loq_half | median | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_loq_half | sd | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_zero | cv | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_zero | mean | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_zero | median | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_zero | sd | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |

### n=500, theoretical censoring=10%

| Method | Estimand | Expected | Finite | Usable | Missing/nonfinite | NA | NaN | Inf | Finite, status-excluded | Theoretical non-ID | Structural non-ID | Not estimable in sample | Not returned | Computational failure | Available (%) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| complete_case | cv | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| complete_case | mean | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| complete_case | median | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| complete_case | sd | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| mle_lognormal | cv | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| mle_lognormal | mean | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| mle_lognormal | median | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| mle_lognormal | sd | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| reverse_km | mean | 2000 | 0 | 0 | 2000 | 2000 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0.00 |
| reverse_km | q0.025 | 2000 | 0 | 0 | 2000 | 2000 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0.00 |
| reverse_km | q0.1 | 2000 | 973 | 0 | 1027 | 1027 | 0 | 0 | 973 | 2000 | 0 | 0 | 0 | 0 | 0.00 |
| reverse_km | q0.5 | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| reverse_km | q0.9 | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| reverse_km | q0.975 | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| ros_log | cv | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| ros_log | mean | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| ros_log | median | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| ros_log | sd | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_loq | cv | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_loq | mean | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_loq | median | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_loq | sd | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_loq_half | cv | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_loq_half | mean | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_loq_half | median | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_loq_half | sd | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_zero | cv | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_zero | mean | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_zero | median | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_zero | sd | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |

### n=30, theoretical censoring=40%

| Method | Estimand | Expected | Finite | Usable | Missing/nonfinite | NA | NaN | Inf | Finite, status-excluded | Theoretical non-ID | Structural non-ID | Not estimable in sample | Not returned | Computational failure | Available (%) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| complete_case | cv | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| complete_case | mean | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| complete_case | median | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| complete_case | sd | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| mle_lognormal | cv | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| mle_lognormal | mean | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| mle_lognormal | median | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| mle_lognormal | sd | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| reverse_km | mean | 2000 | 0 | 0 | 2000 | 2000 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0.00 |
| reverse_km | q0.025 | 2000 | 0 | 0 | 2000 | 2000 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0.00 |
| reverse_km | q0.1 | 2000 | 0 | 0 | 2000 | 2000 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0.00 |
| reverse_km | q0.5 | 2000 | 1654 | 1654 | 346 | 346 | 0 | 0 | 0 | 0 | 0 | 346 | 0 | 0 | 82.70 |
| reverse_km | q0.9 | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| reverse_km | q0.975 | 2000 | 0 | 0 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0.00 |
| ros_log | cv | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| ros_log | mean | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| ros_log | median | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| ros_log | sd | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_loq | cv | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_loq | mean | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_loq | median | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_loq | sd | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_loq_half | cv | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_loq_half | mean | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_loq_half | median | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_loq_half | sd | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_zero | cv | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_zero | mean | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_zero | median | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_zero | sd | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |

### n=100, theoretical censoring=40%

| Method | Estimand | Expected | Finite | Usable | Missing/nonfinite | NA | NaN | Inf | Finite, status-excluded | Theoretical non-ID | Structural non-ID | Not estimable in sample | Not returned | Computational failure | Available (%) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| complete_case | cv | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| complete_case | mean | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| complete_case | median | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| complete_case | sd | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| mle_lognormal | cv | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| mle_lognormal | mean | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| mle_lognormal | median | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| mle_lognormal | sd | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| reverse_km | mean | 2000 | 0 | 0 | 2000 | 2000 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0.00 |
| reverse_km | q0.025 | 2000 | 0 | 0 | 2000 | 2000 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0.00 |
| reverse_km | q0.1 | 2000 | 0 | 0 | 2000 | 2000 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0.00 |
| reverse_km | q0.5 | 2000 | 1931 | 1931 | 69 | 69 | 0 | 0 | 0 | 0 | 0 | 69 | 0 | 0 | 96.55 |
| reverse_km | q0.9 | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| reverse_km | q0.975 | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| ros_log | cv | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| ros_log | mean | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| ros_log | median | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| ros_log | sd | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_loq | cv | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_loq | mean | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_loq | median | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_loq | sd | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_loq_half | cv | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_loq_half | mean | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_loq_half | median | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_loq_half | sd | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_zero | cv | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_zero | mean | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_zero | median | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_zero | sd | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |

### n=500, theoretical censoring=40%

| Method | Estimand | Expected | Finite | Usable | Missing/nonfinite | NA | NaN | Inf | Finite, status-excluded | Theoretical non-ID | Structural non-ID | Not estimable in sample | Not returned | Computational failure | Available (%) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| complete_case | cv | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| complete_case | mean | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| complete_case | median | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| complete_case | sd | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| mle_lognormal | cv | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| mle_lognormal | mean | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| mle_lognormal | median | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| mle_lognormal | sd | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| reverse_km | mean | 2000 | 0 | 0 | 2000 | 2000 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0.00 |
| reverse_km | q0.025 | 2000 | 0 | 0 | 2000 | 2000 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0.00 |
| reverse_km | q0.1 | 2000 | 0 | 0 | 2000 | 2000 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0.00 |
| reverse_km | q0.5 | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| reverse_km | q0.9 | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| reverse_km | q0.975 | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| ros_log | cv | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| ros_log | mean | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| ros_log | median | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| ros_log | sd | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_loq | cv | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_loq | mean | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_loq | median | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_loq | sd | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_loq_half | cv | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_loq_half | mean | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_loq_half | median | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_loq_half | sd | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_zero | cv | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_zero | mean | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_zero | median | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| substitute_zero | sd | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |

### Relative-error denominators and exclusions

Relative errors exclude rows with ineligible status, non-finite estimate/reference/error, or an absolute reference no larger than 1.490116e-08. Counts below are mutually exclusive by exclusion reason and separate the two references.

#### n=30, theoretical censoring=10%

| Method | Estimand | N sample | Status | Estimate | Reference | Error | Near-zero ref | N population | Status | Estimate | Reference | Error | Near-zero ref |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| complete_case | cv | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| complete_case | mean | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| complete_case | median | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| complete_case | sd | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| mle_lognormal | cv | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| mle_lognormal | mean | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| mle_lognormal | median | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| mle_lognormal | sd | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| reverse_km | mean | 0 | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 |
| reverse_km | q0.025 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 |
| reverse_km | q0.1 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 |
| reverse_km | q0.5 | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| reverse_km | q0.9 | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| reverse_km | q0.975 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 |
| ros_log | cv | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| ros_log | mean | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| ros_log | median | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| ros_log | sd | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_loq | cv | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_loq | mean | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_loq | median | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_loq | sd | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_loq_half | cv | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_loq_half | mean | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_loq_half | median | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_loq_half | sd | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_zero | cv | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_zero | mean | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_zero | median | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_zero | sd | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |

#### n=100, theoretical censoring=10%

| Method | Estimand | N sample | Status | Estimate | Reference | Error | Near-zero ref | N population | Status | Estimate | Reference | Error | Near-zero ref |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| complete_case | cv | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| complete_case | mean | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| complete_case | median | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| complete_case | sd | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| mle_lognormal | cv | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| mle_lognormal | mean | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| mle_lognormal | median | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| mle_lognormal | sd | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| reverse_km | mean | 0 | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 |
| reverse_km | q0.025 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 |
| reverse_km | q0.1 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 |
| reverse_km | q0.5 | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| reverse_km | q0.9 | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| reverse_km | q0.975 | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| ros_log | cv | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| ros_log | mean | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| ros_log | median | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| ros_log | sd | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_loq | cv | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_loq | mean | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_loq | median | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_loq | sd | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_loq_half | cv | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_loq_half | mean | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_loq_half | median | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_loq_half | sd | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_zero | cv | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_zero | mean | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_zero | median | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_zero | sd | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |

#### n=500, theoretical censoring=10%

| Method | Estimand | N sample | Status | Estimate | Reference | Error | Near-zero ref | N population | Status | Estimate | Reference | Error | Near-zero ref |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| complete_case | cv | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| complete_case | mean | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| complete_case | median | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| complete_case | sd | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| mle_lognormal | cv | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| mle_lognormal | mean | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| mle_lognormal | median | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| mle_lognormal | sd | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| reverse_km | mean | 0 | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 |
| reverse_km | q0.025 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 |
| reverse_km | q0.1 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 |
| reverse_km | q0.5 | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| reverse_km | q0.9 | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| reverse_km | q0.975 | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| ros_log | cv | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| ros_log | mean | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| ros_log | median | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| ros_log | sd | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_loq | cv | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_loq | mean | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_loq | median | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_loq | sd | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_loq_half | cv | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_loq_half | mean | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_loq_half | median | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_loq_half | sd | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_zero | cv | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_zero | mean | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_zero | median | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_zero | sd | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |

#### n=30, theoretical censoring=40%

| Method | Estimand | N sample | Status | Estimate | Reference | Error | Near-zero ref | N population | Status | Estimate | Reference | Error | Near-zero ref |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| complete_case | cv | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| complete_case | mean | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| complete_case | median | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| complete_case | sd | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| mle_lognormal | cv | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| mle_lognormal | mean | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| mle_lognormal | median | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| mle_lognormal | sd | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| reverse_km | mean | 0 | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 |
| reverse_km | q0.025 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 |
| reverse_km | q0.1 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 |
| reverse_km | q0.5 | 1654 | 346 | 0 | 0 | 0 | 0 | 1654 | 346 | 0 | 0 | 0 | 0 |
| reverse_km | q0.9 | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| reverse_km | q0.975 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 |
| ros_log | cv | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| ros_log | mean | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| ros_log | median | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| ros_log | sd | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_loq | cv | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_loq | mean | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_loq | median | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_loq | sd | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_loq_half | cv | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_loq_half | mean | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_loq_half | median | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_loq_half | sd | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_zero | cv | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_zero | mean | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_zero | median | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_zero | sd | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |

#### n=100, theoretical censoring=40%

| Method | Estimand | N sample | Status | Estimate | Reference | Error | Near-zero ref | N population | Status | Estimate | Reference | Error | Near-zero ref |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| complete_case | cv | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| complete_case | mean | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| complete_case | median | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| complete_case | sd | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| mle_lognormal | cv | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| mle_lognormal | mean | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| mle_lognormal | median | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| mle_lognormal | sd | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| reverse_km | mean | 0 | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 |
| reverse_km | q0.025 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 |
| reverse_km | q0.1 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 |
| reverse_km | q0.5 | 1931 | 69 | 0 | 0 | 0 | 0 | 1931 | 69 | 0 | 0 | 0 | 0 |
| reverse_km | q0.9 | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| reverse_km | q0.975 | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| ros_log | cv | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| ros_log | mean | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| ros_log | median | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| ros_log | sd | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_loq | cv | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_loq | mean | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_loq | median | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_loq | sd | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_loq_half | cv | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_loq_half | mean | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_loq_half | median | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_loq_half | sd | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_zero | cv | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_zero | mean | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_zero | median | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_zero | sd | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |

#### n=500, theoretical censoring=40%

| Method | Estimand | N sample | Status | Estimate | Reference | Error | Near-zero ref | N population | Status | Estimate | Reference | Error | Near-zero ref |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| complete_case | cv | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| complete_case | mean | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| complete_case | median | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| complete_case | sd | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| mle_lognormal | cv | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| mle_lognormal | mean | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| mle_lognormal | median | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| mle_lognormal | sd | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| reverse_km | mean | 0 | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 |
| reverse_km | q0.025 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 |
| reverse_km | q0.1 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 |
| reverse_km | q0.5 | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| reverse_km | q0.9 | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| reverse_km | q0.975 | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| ros_log | cv | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| ros_log | mean | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| ros_log | median | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| ros_log | sd | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_loq | cv | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_loq | mean | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_loq | median | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_loq | sd | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_loq_half | cv | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_loq_half | mean | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_loq_half | median | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_loq_half | sd | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_zero | cv | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_zero | mean | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_zero | median | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |
| substitute_zero | sd | 2000 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0 |

## Comparison with the complete sample

This paired comparison uses each method's estimate minus the complete-sample statistic from the same replicate. The signed mean error is the average method-minus-sample difference; it is not population bias. Each table shows the error distribution and relative-error summaries only among eligible, finite records.

### n=30, theoretical censoring=10%

| Method | Estimand | N used | Relative N | Mean error | Median error | SD error | RMSE | Mean relative error | Median relative error | Error Q2.5% | Error Q25% | Error Q75% | Error Q97.5% |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| complete_case | cv | 2000 | 2000 | -0.08597 | -0.08215 | 0.04924 | 0.09907 | -0.07954 | -0.07394 | -0.19886 | -0.11588 | -0.05261 | 0.00000 |
| complete_case | mean | 2000 | 2000 | 0.16411 | 0.14973 | 0.10173 | 0.19307 | 0.10140 | 0.09723 | 0.00000 | 0.09173 | 0.21699 | 0.41587 |
| complete_case | median | 2000 | 2000 | 0.14184 | 0.10376 | 0.13206 | 0.19378 | 0.14826 | 0.10018 | 0.00000 | 0.04914 | 0.19878 | 0.49997 |
| complete_case | sd | 2000 | 2000 | 0.03476 | 0.01488 | 0.06911 | 0.07734 | 0.01133 | 0.00869 | -0.03589 | -0.00069 | 0.05065 | 0.22261 |
| mle_lognormal | cv | 2000 | 2000 | 0.19601 | 0.19983 | 0.28625 | 0.34687 | 0.21009 | 0.20160 | -0.44730 | 0.05503 | 0.35699 | 0.74736 |
| mle_lognormal | mean | 2000 | 2000 | 0.02818 | 0.03388 | 0.11943 | 0.12268 | 0.01930 | 0.02317 | -0.25281 | -0.01484 | 0.08386 | 0.24303 |
| mle_lognormal | median | 2000 | 2000 | -0.01591 | -0.00730 | 0.13794 | 0.13882 | 0.00175 | -0.00702 | -0.31891 | -0.09175 | 0.07721 | 0.23330 |
| mle_lognormal | sd | 2000 | 2000 | 0.37218 | 0.33424 | 0.74651 | 0.83397 | 0.24691 | 0.22510 | -1.09540 | 0.06294 | 0.69142 | 1.95112 |
| reverse_km | mean | 0 | 0 | — | — | — | — | — | — | — | — | — | — |
| reverse_km | q0.025 | 0 | 0 | — | — | — | — | — | — | — | — | — | — |
| reverse_km | q0.1 | 0 | 0 | — | — | — | — | — | — | — | — | — | — |
| reverse_km | q0.5 | 2000 | 2000 | -0.04284 | -0.02858 | 0.04428 | 0.06160 | -0.04140 | -0.02933 | -0.16053 | -0.05847 | -0.01160 | -0.00095 |
| reverse_km | q0.9 | 2000 | 2000 | -0.07371 | -0.04927 | 0.08164 | 0.10998 | -0.02128 | -0.01454 | -0.29971 | -0.09596 | -0.01934 | -0.00153 |
| reverse_km | q0.975 | 0 | 0 | — | — | — | — | — | — | — | — | — | — |
| ros_log | cv | 2000 | 2000 | 0.00006 | 0.00064 | 0.00877 | 0.00877 | -0.00069 | 0.00050 | -0.01877 | -0.00336 | 0.00443 | 0.01599 |
| ros_log | mean | 2000 | 2000 | -0.00012 | -0.00051 | 0.00669 | 0.00669 | 0.00005 | -0.00032 | -0.01222 | -0.00381 | 0.00273 | 0.01477 |
| ros_log | median | 2000 | 2000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 |
| ros_log | sd | 2000 | 2000 | -0.00037 | 0.00038 | 0.00590 | 0.00591 | -0.00068 | 0.00019 | -0.01412 | -0.00252 | 0.00304 | 0.00874 |
| substitute_loq | cv | 2000 | 2000 | -0.01163 | -0.00999 | 0.00869 | 0.01452 | -0.01104 | -0.00916 | -0.03365 | -0.01632 | -0.00532 | 0.00000 |
| substitute_loq | mean | 2000 | 2000 | 0.00916 | 0.00837 | 0.00606 | 0.01098 | 0.00601 | 0.00507 | 0.00000 | 0.00468 | 0.01260 | 0.02361 |
| substitute_loq | median | 2000 | 2000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 |
| substitute_loq | sd | 2000 | 2000 | -0.00765 | -0.00684 | 0.00517 | 0.00923 | -0.00514 | -0.00385 | -0.01974 | -0.01070 | -0.00384 | 0.00000 |
| substitute_loq_half | cv | 2000 | 2000 | 0.00611 | 0.00522 | 0.00604 | 0.00859 | 0.00587 | 0.00482 | -0.00280 | 0.00195 | 0.00922 | 0.02048 |
| substitute_loq_half | mean | 2000 | 2000 | -0.00476 | -0.00431 | 0.00433 | 0.00643 | -0.00312 | -0.00263 | -0.01447 | -0.00736 | -0.00164 | 0.00227 |
| substitute_loq_half | median | 2000 | 2000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 |
| substitute_loq_half | sd | 2000 | 2000 | 0.00399 | 0.00367 | 0.00372 | 0.00546 | 0.00271 | 0.00196 | -0.00211 | 0.00139 | 0.00621 | 0.01247 |
| substitute_zero | cv | 2000 | 2000 | 0.02501 | 0.02220 | 0.01700 | 0.03024 | 0.02393 | 0.02045 | 0.00000 | 0.01292 | 0.03345 | 0.06760 |
| substitute_zero | mean | 2000 | 2000 | -0.01867 | -0.01751 | 0.01066 | 0.02150 | -0.01225 | -0.01078 | -0.04201 | -0.02529 | -0.01083 | 0.00000 |
| substitute_zero | median | 2000 | 2000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 |
| substitute_zero | sd | 2000 | 2000 | 0.01659 | 0.01549 | 0.00953 | 0.01913 | 0.01125 | 0.00892 | 0.00000 | 0.00935 | 0.02247 | 0.03800 |

### n=100, theoretical censoring=10%

| Method | Estimand | N used | Relative N | Mean error | Median error | SD error | RMSE | Mean relative error | Median relative error | Error Q2.5% | Error Q25% | Error Q75% | Error Q97.5% |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| complete_case | cv | 2000 | 2000 | -0.08658 | -0.08434 | 0.02759 | 0.09087 | -0.07279 | -0.07100 | -0.14617 | -0.10356 | -0.06696 | -0.03964 |
| complete_case | mean | 2000 | 2000 | 0.16368 | 0.15863 | 0.05466 | 0.17256 | 0.10053 | 0.09823 | 0.07298 | 0.12566 | 0.19711 | 0.28295 |
| complete_case | median | 2000 | 2000 | 0.13608 | 0.12464 | 0.07300 | 0.15442 | 0.13864 | 0.12229 | 0.03231 | 0.08247 | 0.17650 | 0.31203 |
| complete_case | sd | 2000 | 2000 | 0.04641 | 0.03453 | 0.04720 | 0.06619 | 0.01976 | 0.01825 | -0.00581 | 0.01600 | 0.06332 | 0.17142 |
| mle_lognormal | cv | 2000 | 2000 | 0.09199 | 0.13731 | 0.25265 | 0.26882 | 0.10615 | 0.12424 | -0.55852 | -0.00019 | 0.24056 | 0.44107 |
| mle_lognormal | mean | 2000 | 2000 | 0.00855 | 0.01745 | 0.07028 | 0.07078 | 0.00651 | 0.01117 | -0.15393 | -0.01981 | 0.04942 | 0.11925 |
| mle_lognormal | median | 2000 | 2000 | -0.00569 | -0.00285 | 0.07234 | 0.07255 | -0.00048 | -0.00263 | -0.15510 | -0.05169 | 0.04340 | 0.12529 |
| mle_lognormal | sd | 2000 | 2000 | 0.15302 | 0.23658 | 0.55545 | 0.57601 | 0.11981 | 0.13435 | -1.11732 | -0.01874 | 0.44431 | 0.93472 |
| reverse_km | mean | 0 | 0 | — | — | — | — | — | — | — | — | — | — |
| reverse_km | q0.025 | 0 | 0 | — | — | — | — | — | — | — | — | — | — |
| reverse_km | q0.1 | 0 | 0 | — | — | — | — | — | — | — | — | — | — |
| reverse_km | q0.5 | 2000 | 2000 | -0.01264 | -0.00841 | 0.01285 | 0.01802 | -0.01252 | -0.00835 | -0.04601 | -0.01732 | -0.00351 | -0.00026 |
| reverse_km | q0.9 | 2000 | 2000 | -0.02134 | -0.01428 | 0.02213 | 0.03074 | -0.00609 | -0.00424 | -0.08069 | -0.02934 | -0.00602 | -0.00052 |
| reverse_km | q0.975 | 2000 | 2000 | 0.47707 | 0.30063 | 0.53245 | 0.71481 | 0.06714 | 0.04879 | 0.00995 | 0.12544 | 0.63118 | 1.91814 |
| ros_log | cv | 2000 | 2000 | 0.00016 | 0.00057 | 0.00413 | 0.00413 | -0.00011 | 0.00045 | -0.00890 | -0.00203 | 0.00283 | 0.00732 |
| ros_log | mean | 2000 | 2000 | -0.00015 | -0.00046 | 0.00335 | 0.00335 | -0.00004 | -0.00029 | -0.00602 | -0.00229 | 0.00168 | 0.00743 |
| ros_log | median | 2000 | 2000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 |
| ros_log | sd | 2000 | 2000 | -0.00007 | 0.00031 | 0.00262 | 0.00262 | -0.00016 | 0.00015 | -0.00611 | -0.00134 | 0.00162 | 0.00410 |
| substitute_loq | cv | 2000 | 2000 | -0.01126 | -0.01072 | 0.00447 | 0.01211 | -0.00957 | -0.00912 | -0.02140 | -0.01394 | -0.00818 | -0.00378 |
| substitute_loq | mean | 2000 | 2000 | 0.00924 | 0.00900 | 0.00330 | 0.00981 | 0.00578 | 0.00553 | 0.00336 | 0.00694 | 0.01134 | 0.01637 |
| substitute_loq | median | 2000 | 2000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 |
| substitute_loq | sd | 2000 | 2000 | -0.00691 | -0.00663 | 0.00267 | 0.00740 | -0.00386 | -0.00352 | -0.01290 | -0.00855 | -0.00507 | -0.00239 |
| substitute_loq_half | cv | 2000 | 2000 | 0.00591 | 0.00567 | 0.00303 | 0.00664 | 0.00504 | 0.00470 | 0.00079 | 0.00378 | 0.00770 | 0.01277 |
| substitute_loq_half | mean | 2000 | 2000 | -0.00481 | -0.00467 | 0.00232 | 0.00534 | -0.00300 | -0.00287 | -0.00985 | -0.00625 | -0.00319 | -0.00074 |
| substitute_loq_half | median | 2000 | 2000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 |
| substitute_loq_half | sd | 2000 | 2000 | 0.00360 | 0.00342 | 0.00184 | 0.00404 | 0.00201 | 0.00178 | 0.00047 | 0.00232 | 0.00472 | 0.00764 |
| substitute_zero | cv | 2000 | 2000 | 0.02399 | 0.02283 | 0.00845 | 0.02543 | 0.02044 | 0.01947 | 0.00990 | 0.01817 | 0.02940 | 0.04236 |
| substitute_zero | mean | 2000 | 2000 | -0.01886 | -0.01847 | 0.00577 | 0.01972 | -0.01179 | -0.01130 | -0.03089 | -0.02264 | -0.01473 | -0.00853 |
| substitute_zero | median | 2000 | 2000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 |
| substitute_zero | sd | 2000 | 2000 | 0.01496 | 0.01466 | 0.00499 | 0.01577 | 0.00838 | 0.00775 | 0.00627 | 0.01130 | 0.01812 | 0.02567 |

### n=500, theoretical censoring=10%

| Method | Estimand | N used | Relative N | Mean error | Median error | SD error | RMSE | Mean relative error | Median relative error | Error Q2.5% | Error Q25% | Error Q75% | Error Q97.5% |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| complete_case | cv | 2000 | 2000 | -0.08566 | -0.08497 | 0.01305 | 0.08664 | -0.06749 | -0.06706 | -0.11280 | -0.09412 | -0.07617 | -0.06314 |
| complete_case | mean | 2000 | 2000 | 0.16221 | 0.16132 | 0.02470 | 0.16408 | 0.09861 | 0.09790 | 0.11830 | 0.14475 | 0.17849 | 0.21412 |
| complete_case | median | 2000 | 2000 | 0.13305 | 0.13018 | 0.03251 | 0.13697 | 0.13344 | 0.13064 | 0.07703 | 0.10974 | 0.15504 | 0.20334 |
| complete_case | sd | 2000 | 2000 | 0.05399 | 0.04853 | 0.02770 | 0.06068 | 0.02433 | 0.02357 | 0.01803 | 0.03581 | 0.06606 | 0.11651 |
| mle_lognormal | cv | 2000 | 2000 | 0.02621 | 0.06231 | 0.18388 | 0.18569 | 0.03613 | 0.05036 | -0.42158 | -0.03212 | 0.13521 | 0.24584 |
| mle_lognormal | mean | 2000 | 2000 | 0.00159 | 0.00552 | 0.03226 | 0.03229 | 0.00131 | 0.00332 | -0.07166 | -0.01540 | 0.02299 | 0.05513 |
| mle_lognormal | median | 2000 | 2000 | -0.00179 | -0.00167 | 0.03289 | 0.03293 | -0.00070 | -0.00169 | -0.06764 | -0.02276 | 0.02044 | 0.06194 |
| mle_lognormal | sd | 2000 | 2000 | 0.04047 | 0.10673 | 0.35484 | 0.35705 | 0.03956 | 0.05307 | -0.80509 | -0.07100 | 0.24633 | 0.47932 |
| reverse_km | mean | 0 | 0 | — | — | — | — | — | — | — | — | — | — |
| reverse_km | q0.025 | 0 | 0 | — | — | — | — | — | — | — | — | — | — |
| reverse_km | q0.1 | 0 | 0 | — | — | — | — | — | — | — | — | — | — |
| reverse_km | q0.5 | 2000 | 2000 | -0.00250 | -0.00180 | 0.00251 | 0.00354 | -0.00249 | -0.00181 | -0.00905 | -0.00341 | -0.00075 | -0.00007 |
| reverse_km | q0.9 | 2000 | 2000 | -0.00414 | -0.00283 | 0.00414 | 0.00585 | -0.00116 | -0.00079 | -0.01590 | -0.00573 | -0.00119 | -0.00012 |
| reverse_km | q0.975 | 2000 | 2000 | 0.11397 | 0.07487 | 0.11875 | 0.16457 | 0.01607 | 0.01090 | 0.00240 | 0.02979 | 0.15587 | 0.44388 |
| ros_log | cv | 2000 | 2000 | 0.00010 | 0.00018 | 0.00169 | 0.00169 | 0.00003 | 0.00014 | -0.00348 | -0.00099 | 0.00125 | 0.00327 |
| ros_log | mean | 2000 | 2000 | -0.00008 | -0.00015 | 0.00140 | 0.00140 | -0.00004 | -0.00009 | -0.00264 | -0.00104 | 0.00082 | 0.00282 |
| ros_log | median | 2000 | 2000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 |
| ros_log | sd | 2000 | 2000 | 0.00002 | 0.00011 | 0.00099 | 0.00099 | -0.00001 | 0.00005 | -0.00216 | -0.00059 | 0.00069 | 0.00177 |
| substitute_loq | cv | 2000 | 2000 | -0.01101 | -0.01086 | 0.00205 | 0.01120 | -0.00869 | -0.00861 | -0.01540 | -0.01234 | -0.00950 | -0.00741 |
| substitute_loq | mean | 2000 | 2000 | 0.00918 | 0.00908 | 0.00153 | 0.00931 | 0.00560 | 0.00554 | 0.00647 | 0.00810 | 0.01021 | 0.01238 |
| substitute_loq | median | 2000 | 2000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 |
| substitute_loq | sd | 2000 | 2000 | -0.00637 | -0.00633 | 0.00126 | 0.00650 | -0.00314 | -0.00309 | -0.00894 | -0.00716 | -0.00556 | -0.00395 |
| substitute_loq_half | cv | 2000 | 2000 | 0.00566 | 0.00561 | 0.00138 | 0.00582 | 0.00447 | 0.00439 | 0.00318 | 0.00468 | 0.00657 | 0.00856 |
| substitute_loq_half | mean | 2000 | 2000 | -0.00468 | -0.00465 | 0.00108 | 0.00480 | -0.00286 | -0.00283 | -0.00701 | -0.00536 | -0.00394 | -0.00267 |
| substitute_loq_half | median | 2000 | 2000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 |
| substitute_loq_half | sd | 2000 | 2000 | 0.00325 | 0.00322 | 0.00085 | 0.00336 | 0.00161 | 0.00155 | 0.00171 | 0.00265 | 0.00379 | 0.00506 |
| substitute_zero | cv | 2000 | 2000 | 0.02310 | 0.02288 | 0.00383 | 0.02341 | 0.01827 | 0.01808 | 0.01641 | 0.02032 | 0.02561 | 0.03107 |
| substitute_zero | mean | 2000 | 2000 | -0.01854 | -0.01849 | 0.00267 | 0.01873 | -0.01131 | -0.01125 | -0.02396 | -0.02036 | -0.01655 | -0.01373 |
| substitute_zero | median | 2000 | 2000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 |
| substitute_zero | sd | 2000 | 2000 | 0.01366 | 0.01359 | 0.00246 | 0.01388 | 0.00674 | 0.00659 | 0.00887 | 0.01204 | 0.01535 | 0.01850 |

### n=30, theoretical censoring=40%

| Method | Estimand | N used | Relative N | Mean error | Median error | SD error | RMSE | Mean relative error | Median relative error | Error Q2.5% | Error Q25% | Error Q75% | Error Q97.5% |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| complete_case | cv | 2000 | 2000 | -0.30785 | -0.29805 | 0.09539 | 0.32228 | -0.28331 | -0.27688 | -0.51560 | -0.36348 | -0.24290 | -0.15374 |
| complete_case | mean | 2000 | 2000 | 0.80500 | 0.74929 | 0.31786 | 0.86545 | 0.49910 | 0.47423 | 0.36846 | 0.59000 | 0.94552 | 1.57620 |
| complete_case | median | 2000 | 2000 | 0.71189 | 0.63080 | 0.39554 | 0.81435 | 0.76574 | 0.62048 | 0.17264 | 0.43180 | 0.91385 | 1.69238 |
| complete_case | sd | 2000 | 2000 | 0.19024 | 0.10414 | 0.33794 | 0.38774 | 0.06548 | 0.06152 | -0.14517 | -0.00821 | 0.27846 | 1.01290 |
| mle_lognormal | cv | 2000 | 2000 | 0.20772 | 0.16290 | 0.35028 | 0.40716 | 0.19583 | 0.16470 | -0.38939 | 0.03139 | 0.34472 | 1.00388 |
| mle_lognormal | mean | 2000 | 2000 | 0.02933 | 0.03513 | 0.10666 | 0.11060 | 0.01932 | 0.02384 | -0.19967 | -0.00445 | 0.07125 | 0.20818 |
| mle_lognormal | median | 2000 | 2000 | -0.01603 | -0.00512 | 0.10059 | 0.10184 | -0.00791 | -0.00574 | -0.24702 | -0.06510 | 0.04517 | 0.16642 |
| mle_lognormal | sd | 2000 | 2000 | 0.40736 | 0.27829 | 0.89621 | 0.98424 | 0.22854 | 0.18612 | -1.03690 | 0.05216 | 0.65073 | 2.44450 |
| reverse_km | mean | 0 | 0 | — | — | — | — | — | — | — | — | — | — |
| reverse_km | q0.025 | 0 | 0 | — | — | — | — | — | — | — | — | — | — |
| reverse_km | q0.1 | 0 | 0 | — | — | — | — | — | — | — | — | — | — |
| reverse_km | q0.5 | 1654 | 1654 | -0.04348 | -0.03020 | 0.04404 | 0.06188 | -0.03900 | -0.02853 | -0.16920 | -0.05899 | -0.01220 | -0.00094 |
| reverse_km | q0.9 | 2000 | 2000 | -0.07454 | -0.04945 | 0.08151 | 0.11044 | -0.02180 | -0.01513 | -0.28662 | -0.09761 | -0.01937 | -0.00214 |
| reverse_km | q0.975 | 0 | 0 | — | — | — | — | — | — | — | — | — | — |
| ros_log | cv | 2000 | 2000 | 0.01233 | 0.00824 | 0.05987 | 0.06111 | 0.00364 | 0.00772 | -0.09337 | -0.02344 | 0.04330 | 0.13834 |
| ros_log | mean | 2000 | 2000 | -0.00920 | -0.00845 | 0.04696 | 0.04784 | -0.00464 | -0.00487 | -0.10634 | -0.03763 | 0.02137 | 0.07930 |
| ros_log | median | 2000 | 2000 | -0.00440 | 0.00000 | 0.04791 | 0.04810 | -0.00647 | 0.00000 | -0.13380 | 0.00000 | 0.00000 | 0.08256 |
| ros_log | sd | 2000 | 2000 | 0.00149 | 0.00532 | 0.03174 | 0.03177 | -0.00277 | 0.00255 | -0.06856 | -0.01663 | 0.02325 | 0.05354 |
| substitute_loq | cv | 2000 | 2000 | -0.14616 | -0.13922 | 0.05701 | 0.15688 | -0.13770 | -0.12705 | -0.27832 | -0.17740 | -0.10492 | -0.05582 |
| substitute_loq | mean | 2000 | 2000 | 0.13695 | 0.13486 | 0.03918 | 0.14244 | 0.09010 | 0.08367 | 0.06428 | 0.10921 | 0.16295 | 0.21824 |
| substitute_loq | median | 2000 | 2000 | 0.01175 | 0.00000 | 0.03808 | 0.03985 | 0.01852 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.13754 |
| substitute_loq | sd | 2000 | 2000 | -0.09334 | -0.09151 | 0.02741 | 0.09728 | -0.06232 | -0.05400 | -0.15230 | -0.11086 | -0.07377 | -0.04425 |
| substitute_loq_half | cv | 2000 | 2000 | 0.01667 | 0.01549 | 0.02708 | 0.03179 | 0.01619 | 0.01438 | -0.03744 | -0.00022 | 0.03290 | 0.07327 |
| substitute_loq_half | mean | 2000 | 2000 | -0.01708 | -0.01657 | 0.02239 | 0.02815 | -0.01087 | -0.00986 | -0.06320 | -0.03185 | -0.00230 | 0.02734 |
| substitute_loq_half | median | 2000 | 2000 | -0.03958 | 0.00000 | 0.09188 | 0.10002 | -0.05521 | 0.00000 | -0.32887 | 0.00000 | 0.00000 | 0.00000 |
| substitute_loq_half | sd | 2000 | 2000 | 0.00751 | 0.00731 | 0.01585 | 0.01754 | 0.00474 | 0.00364 | -0.02421 | -0.00252 | 0.01742 | 0.03927 |
| substitute_zero | cv | 2000 | 2000 | 0.23406 | 0.21538 | 0.10409 | 0.25615 | 0.22328 | 0.19778 | 0.08818 | 0.16107 | 0.28704 | 0.49233 |
| substitute_zero | mean | 2000 | 2000 | -0.17111 | -0.16950 | 0.04525 | 0.17699 | -0.11185 | -0.10427 | -0.26825 | -0.20022 | -0.14129 | -0.08442 |
| substitute_zero | median | 2000 | 2000 | -0.09090 | 0.00000 | 0.20948 | 0.22830 | -0.12893 | 0.00000 | -0.71697 | 0.00000 | 0.00000 | 0.00000 |
| substitute_zero | sd | 2000 | 2000 | 0.12220 | 0.12198 | 0.03055 | 0.12596 | 0.08125 | 0.07154 | 0.06277 | 0.10154 | 0.14298 | 0.18222 |

### n=100, theoretical censoring=40%

| Method | Estimand | N used | Relative N | Mean error | Median error | SD error | RMSE | Mean relative error | Median relative error | Error Q2.5% | Error Q25% | Error Q75% | Error Q97.5% |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| complete_case | cv | 2000 | 2000 | -0.31073 | -0.30491 | 0.05868 | 0.31622 | -0.25796 | -0.25639 | -0.43242 | -0.34167 | -0.27327 | -0.21349 |
| complete_case | mean | 2000 | 2000 | 0.81963 | 0.79979 | 0.16758 | 0.83658 | 0.49842 | 0.48961 | 0.54967 | 0.70272 | 0.91739 | 1.20095 |
| complete_case | median | 2000 | 2000 | 0.70673 | 0.68425 | 0.20687 | 0.73637 | 0.72525 | 0.67922 | 0.37701 | 0.55585 | 0.83031 | 1.19117 |
| complete_case | sd | 2000 | 2000 | 0.26006 | 0.20421 | 0.25795 | 0.36625 | 0.10943 | 0.10579 | -0.02042 | 0.11222 | 0.33858 | 0.84712 |
| mle_lognormal | cv | 2000 | 2000 | 0.09679 | 0.11645 | 0.25987 | 0.27724 | 0.10266 | 0.10473 | -0.46255 | -0.00641 | 0.23726 | 0.49559 |
| mle_lognormal | mean | 2000 | 2000 | 0.00691 | 0.01722 | 0.06662 | 0.06696 | 0.00546 | 0.01053 | -0.13586 | -0.01653 | 0.04277 | 0.09913 |
| mle_lognormal | median | 2000 | 2000 | -0.00390 | 0.00026 | 0.05318 | 0.05331 | -0.00137 | 0.00027 | -0.12138 | -0.03490 | 0.02973 | 0.09498 |
| mle_lognormal | sd | 2000 | 2000 | 0.15942 | 0.20654 | 0.59890 | 0.61961 | 0.11356 | 0.11815 | -1.01067 | -0.02512 | 0.43299 | 1.02581 |
| reverse_km | mean | 0 | 0 | — | — | — | — | — | — | — | — | — | — |
| reverse_km | q0.025 | 0 | 0 | — | — | — | — | — | — | — | — | — | — |
| reverse_km | q0.1 | 0 | 0 | — | — | — | — | — | — | — | — | — | — |
| reverse_km | q0.5 | 1931 | 1931 | -0.01305 | -0.00900 | 0.01343 | 0.01872 | -0.01284 | -0.00903 | -0.04866 | -0.01778 | -0.00369 | -0.00031 |
| reverse_km | q0.9 | 2000 | 2000 | -0.02079 | -0.01382 | 0.02147 | 0.02988 | -0.00585 | -0.00395 | -0.07487 | -0.02877 | -0.00582 | -0.00044 |
| reverse_km | q0.975 | 2000 | 2000 | 0.47646 | 0.32303 | 0.50986 | 0.69774 | 0.06618 | 0.04894 | 0.01152 | 0.12201 | 0.64546 | 1.87980 |
| ros_log | cv | 2000 | 2000 | 0.00768 | 0.00689 | 0.02870 | 0.02970 | 0.00383 | 0.00563 | -0.04702 | -0.01116 | 0.02430 | 0.06972 |
| ros_log | mean | 2000 | 2000 | -0.00606 | -0.00621 | 0.02435 | 0.02508 | -0.00331 | -0.00371 | -0.05441 | -0.02154 | 0.00985 | 0.04214 |
| ros_log | median | 2000 | 2000 | 0.00002 | 0.00000 | 0.00980 | 0.00980 | -0.00003 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00440 |
| ros_log | sd | 2000 | 2000 | 0.00216 | 0.00374 | 0.01481 | 0.01496 | 0.00016 | 0.00173 | -0.03038 | -0.00693 | 0.01213 | 0.02766 |
| substitute_loq | cv | 2000 | 2000 | -0.14338 | -0.14002 | 0.03106 | 0.14670 | -0.12018 | -0.11694 | -0.21093 | -0.16215 | -0.12225 | -0.09092 |
| substitute_loq | mean | 2000 | 2000 | 0.13761 | 0.13724 | 0.02074 | 0.13916 | 0.08525 | 0.08304 | 0.09851 | 0.12375 | 0.15132 | 0.17991 |
| substitute_loq | median | 2000 | 2000 | 0.00089 | 0.00000 | 0.00700 | 0.00706 | 0.00124 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00542 |
| substitute_loq | sd | 2000 | 2000 | -0.08369 | -0.08396 | 0.01676 | 0.08535 | -0.04575 | -0.04335 | -0.11744 | -0.09518 | -0.07181 | -0.05171 |
| substitute_loq_half | cv | 2000 | 2000 | 0.01849 | 0.01785 | 0.01402 | 0.02320 | 0.01561 | 0.01456 | -0.00757 | 0.00904 | 0.02732 | 0.04806 |
| substitute_loq_half | mean | 2000 | 2000 | -0.01820 | -0.01797 | 0.01198 | 0.02179 | -0.01123 | -0.01079 | -0.04229 | -0.02594 | -0.01026 | 0.00435 |
| substitute_loq_half | median | 2000 | 2000 | -0.00978 | 0.00000 | 0.05402 | 0.05489 | -0.01305 | 0.00000 | -0.19138 | 0.00000 | 0.00000 | 0.00000 |
| substitute_loq_half | sd | 2000 | 2000 | 0.00760 | 0.00731 | 0.00750 | 0.01067 | 0.00412 | 0.00352 | -0.00644 | 0.00251 | 0.01256 | 0.02368 |
| substitute_zero | cv | 2000 | 2000 | 0.22658 | 0.22007 | 0.05440 | 0.23302 | 0.19092 | 0.18388 | 0.14000 | 0.18774 | 0.25773 | 0.35440 |
| substitute_zero | mean | 2000 | 2000 | -0.17401 | -0.17341 | 0.02445 | 0.17571 | -0.10771 | -0.10557 | -0.22476 | -0.19030 | -0.15658 | -0.12833 |
| substitute_zero | median | 2000 | 2000 | -0.02046 | 0.00000 | 0.11311 | 0.11491 | -0.02733 | 0.00000 | -0.38543 | 0.00000 | 0.00000 | 0.00000 |
| substitute_zero | sd | 2000 | 2000 | 0.11228 | 0.11306 | 0.02019 | 0.11408 | 0.06136 | 0.05875 | 0.07019 | 0.09910 | 0.12645 | 0.14980 |

### n=500, theoretical censoring=40%

| Method | Estimand | N used | Relative N | Mean error | Median error | SD error | RMSE | Mean relative error | Median relative error | Error Q2.5% | Error Q25% | Error Q75% | Error Q97.5% |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| complete_case | cv | 2000 | 2000 | -0.30719 | -0.30537 | 0.02605 | 0.30829 | -0.24082 | -0.24168 | -0.36350 | -0.32240 | -0.28946 | -0.26341 |
| complete_case | mean | 2000 | 2000 | 0.81079 | 0.80608 | 0.07315 | 0.81409 | 0.49266 | 0.49155 | 0.68131 | 0.75996 | 0.85976 | 0.96173 |
| complete_case | median | 2000 | 2000 | 0.68925 | 0.68596 | 0.08824 | 0.69487 | 0.69226 | 0.68636 | 0.52737 | 0.62848 | 0.74604 | 0.87440 |
| complete_case | sd | 2000 | 2000 | 0.29536 | 0.26783 | 0.13977 | 0.32675 | 0.13282 | 0.12942 | 0.10652 | 0.20056 | 0.35966 | 0.63180 |
| mle_lognormal | cv | 2000 | 2000 | 0.02176 | 0.05427 | 0.16114 | 0.16256 | 0.02953 | 0.04361 | -0.37903 | -0.03999 | 0.12447 | 0.23401 |
| mle_lognormal | mean | 2000 | 2000 | 0.00087 | 0.00427 | 0.02728 | 0.02728 | 0.00081 | 0.00260 | -0.06132 | -0.01368 | 0.01932 | 0.04360 |
| mle_lognormal | median | 2000 | 2000 | -0.00159 | -0.00058 | 0.02397 | 0.02402 | -0.00099 | -0.00058 | -0.05079 | -0.01796 | 0.01599 | 0.04157 |
| mle_lognormal | sd | 2000 | 2000 | 0.03312 | 0.09348 | 0.30687 | 0.30857 | 0.03189 | 0.04748 | -0.71512 | -0.08301 | 0.22564 | 0.43722 |
| reverse_km | mean | 0 | 0 | — | — | — | — | — | — | — | — | — | — |
| reverse_km | q0.025 | 0 | 0 | — | — | — | — | — | — | — | — | — | — |
| reverse_km | q0.1 | 0 | 0 | — | — | — | — | — | — | — | — | — | — |
| reverse_km | q0.5 | 2000 | 2000 | -0.00261 | -0.00180 | 0.00251 | 0.00362 | -0.00261 | -0.00182 | -0.00917 | -0.00369 | -0.00079 | -0.00008 |
| reverse_km | q0.9 | 2000 | 2000 | -0.00413 | -0.00293 | 0.00410 | 0.00582 | -0.00116 | -0.00083 | -0.01490 | -0.00585 | -0.00117 | -0.00012 |
| reverse_km | q0.975 | 2000 | 2000 | 0.11372 | 0.07311 | 0.12037 | 0.16557 | 0.01595 | 0.01048 | 0.00266 | 0.03087 | 0.15557 | 0.46393 |
| ros_log | cv | 2000 | 2000 | 0.00207 | 0.00202 | 0.01234 | 0.01251 | 0.00101 | 0.00162 | -0.02157 | -0.00652 | 0.01037 | 0.02694 |
| ros_log | mean | 2000 | 2000 | -0.00164 | -0.00184 | 0.01079 | 0.01091 | -0.00091 | -0.00114 | -0.02261 | -0.00919 | 0.00570 | 0.01991 |
| ros_log | median | 2000 | 2000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 |
| ros_log | sd | 2000 | 2000 | 0.00055 | 0.00098 | 0.00620 | 0.00623 | 0.00004 | 0.00042 | -0.01286 | -0.00342 | 0.00493 | 0.01125 |
| substitute_loq | cv | 2000 | 2000 | -0.14342 | -0.14236 | 0.01465 | 0.14417 | -0.11255 | -0.11245 | -0.17616 | -0.15184 | -0.13372 | -0.11704 |
| substitute_loq | mean | 2000 | 2000 | 0.13750 | 0.13741 | 0.00913 | 0.13781 | 0.08385 | 0.08350 | 0.11944 | 0.13123 | 0.14338 | 0.15548 |
| substitute_loq | median | 2000 | 2000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 |
| substitute_loq | sd | 2000 | 2000 | -0.07804 | -0.07873 | 0.00956 | 0.07863 | -0.03823 | -0.03806 | -0.09480 | -0.08474 | -0.07221 | -0.05762 |
| substitute_loq_half | cv | 2000 | 2000 | 0.01836 | 0.01818 | 0.00611 | 0.01935 | 0.01443 | 0.01428 | 0.00651 | 0.01421 | 0.02240 | 0.03041 |
| substitute_loq_half | mean | 2000 | 2000 | -0.01781 | -0.01768 | 0.00521 | 0.01856 | -0.01086 | -0.01076 | -0.02797 | -0.02133 | -0.01436 | -0.00745 |
| substitute_loq_half | median | 2000 | 2000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 |
| substitute_loq_half | sd | 2000 | 2000 | 0.00695 | 0.00695 | 0.00311 | 0.00761 | 0.00340 | 0.00327 | 0.00090 | 0.00490 | 0.00898 | 0.01321 |
| substitute_zero | cv | 2000 | 2000 | 0.22371 | 0.22196 | 0.02395 | 0.22499 | 0.17581 | 0.17417 | 0.18172 | 0.20764 | 0.23818 | 0.27450 |
| substitute_zero | mean | 2000 | 2000 | -0.17313 | -0.17290 | 0.01049 | 0.17345 | -0.10556 | -0.10503 | -0.19363 | -0.18011 | -0.16644 | -0.15182 |
| substitute_zero | median | 2000 | 2000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 | 0.00000 |
| substitute_zero | sd | 2000 | 2000 | 0.10508 | 0.10636 | 0.01199 | 0.10576 | 0.05147 | 0.05145 | 0.07868 | 0.09788 | 0.11349 | 0.12540 |

## Comparison with theoretical population parameters

These are separate, unconditional-to-the-sample-reference performance summaries: population error is estimate minus the configured population parameter. Bias, relative bias, RMSE, estimator SD, error SD, median error, and error quantiles are not interchangeable with the paired sample-truth errors above.

### n=30, theoretical censoring=10%

| Method | Estimand | N used | Relative N | Bias | Relative bias | RMSE | SD estimates | SD error | Median error | Error Q2.5% | Error Q25% | Error Q75% | Error Q97.5% |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| complete_case | cv | 2000 | 2000 | -0.28686 | -0.21884 | 0.41367 | 0.29813 | 0.29813 | -0.35445 | -0.68630 | -0.48404 | -0.17004 | 0.49322 |
| complete_case | mean | 2000 | 2000 | 0.17218 | 0.10443 | 0.47273 | 0.44037 | 0.44037 | 0.08573 | -0.47759 | -0.13029 | 0.41509 | 1.21330 |
| complete_case | median | 2000 | 2000 | 0.17657 | 0.17657 | 0.31476 | 0.26063 | 0.26063 | 0.14183 | -0.24087 | -0.01931 | 0.33736 | 0.77282 |
| complete_case | sd | 2000 | 2000 | -0.22193 | -0.10269 | 1.02688 | 1.00286 | 1.00286 | -0.47413 | -1.35259 | -0.86720 | 0.09874 | 2.37453 |
| mle_lognormal | cv | 2000 | 2000 | -0.00488 | -0.00372 | 0.31008 | 0.31012 | 0.31012 | -0.05259 | -0.49812 | -0.21965 | 0.15476 | 0.71769 |
| mle_lognormal | mean | 2000 | 2000 | 0.03625 | 0.02199 | 0.40740 | 0.40588 | 0.40588 | -0.03764 | -0.56632 | -0.24777 | 0.26866 | 0.93111 |
| mle_lognormal | median | 2000 | 2000 | 0.01882 | 0.01882 | 0.19189 | 0.19101 | 0.19101 | -0.00328 | -0.30539 | -0.11399 | 0.13287 | 0.44419 |
| mle_lognormal | sd | 2000 | 2000 | 0.11549 | 0.05344 | 1.04618 | 1.04005 | 1.04005 | -0.15059 | -1.18559 | -0.58023 | 0.54196 | 2.82552 |
| reverse_km | mean | 0 | 0 | — | — | — | — | — | — | — | — | — | — |
| reverse_km | q0.025 | 0 | 0 | — | — | — | — | — | — | — | — | — | — |
| reverse_km | q0.1 | 0 | 0 | — | — | — | — | — | — | — | — | — | — |
| reverse_km | q0.5 | 2000 | 2000 | -0.00811 | -0.00811 | 0.23495 | 0.23487 | 0.23487 | -0.03694 | -0.38846 | -0.17747 | 0.12778 | 0.52721 |
| reverse_km | q0.9 | 2000 | 2000 | -0.16003 | -0.04443 | 1.11903 | 1.10780 | 1.10780 | -0.35112 | -1.72396 | -0.93527 | 0.37824 | 2.53117 |
| reverse_km | q0.975 | 0 | 0 | — | — | — | — | — | — | — | — | — | — |
| ros_log | cv | 2000 | 2000 | -0.20083 | -0.15321 | 0.36975 | 0.31053 | 0.31053 | -0.26546 | -0.61936 | -0.40949 | -0.07804 | 0.60000 |
| ros_log | mean | 2000 | 2000 | 0.00794 | 0.00482 | 0.40342 | 0.40344 | 0.40344 | -0.06567 | -0.59049 | -0.27124 | 0.24042 | 0.98900 |
| ros_log | median | 2000 | 2000 | 0.03473 | 0.03473 | 0.24287 | 0.24044 | 0.24044 | 0.00573 | -0.35280 | -0.13821 | 0.17190 | 0.59391 |
| ros_log | sd | 2000 | 2000 | -0.25707 | -0.11895 | 0.98394 | 0.95001 | 0.95001 | -0.49358 | -1.33835 | -0.86825 | 0.05144 | 2.22665 |
| substitute_loq | cv | 2000 | 2000 | -0.21252 | -0.16213 | 0.37346 | 0.30717 | 0.30717 | -0.27871 | -0.62210 | -0.42012 | -0.09318 | 0.56970 |
| substitute_loq | mean | 2000 | 2000 | 0.01722 | 0.01045 | 0.40365 | 0.40338 | 0.40338 | -0.05591 | -0.57809 | -0.26132 | 0.24732 | 0.99650 |
| substitute_loq | median | 2000 | 2000 | 0.03473 | 0.03473 | 0.24287 | 0.24044 | 0.24044 | 0.00573 | -0.35280 | -0.13821 | 0.17190 | 0.59391 |
| substitute_loq | sd | 2000 | 2000 | -0.26434 | -0.12231 | 0.98574 | 0.94987 | 0.94987 | -0.50091 | -1.34119 | -0.87399 | 0.04223 | 2.22336 |
| substitute_loq_half | cv | 2000 | 2000 | -0.19478 | -0.14859 | 0.36396 | 0.30753 | 0.30753 | -0.26083 | -0.60747 | -0.40157 | -0.07211 | 0.59260 |
| substitute_loq_half | mean | 2000 | 2000 | 0.00331 | 0.00201 | 0.40481 | 0.40490 | 0.40490 | -0.06994 | -0.59878 | -0.27674 | 0.23497 | 0.98133 |
| substitute_loq_half | median | 2000 | 2000 | 0.03473 | 0.03473 | 0.24287 | 0.24044 | 0.24044 | 0.00573 | -0.35280 | -0.13821 | 0.17190 | 0.59391 |
| substitute_loq_half | sd | 2000 | 2000 | -0.25270 | -0.11692 | 0.98075 | 0.94787 | 0.94787 | -0.49072 | -1.33071 | -0.86263 | 0.05271 | 2.22592 |
| substitute_zero | cv | 2000 | 2000 | -0.17588 | -0.13417 | 0.35475 | 0.30816 | 0.30816 | -0.24088 | -0.59514 | -0.38187 | -0.05138 | 0.60913 |
| substitute_zero | mean | 2000 | 2000 | -0.01060 | -0.00643 | 0.40658 | 0.40655 | 0.40655 | -0.08402 | -0.62417 | -0.29268 | 0.22432 | 0.96756 |
| substitute_zero | median | 2000 | 2000 | 0.03473 | 0.03473 | 0.24287 | 0.24044 | 0.24044 | 0.00573 | -0.35280 | -0.13821 | 0.17190 | 0.59391 |
| substitute_zero | sd | 2000 | 2000 | -0.24010 | -0.11110 | 0.97539 | 0.94562 | 0.94562 | -0.47968 | -1.30982 | -0.84942 | 0.06190 | 2.23168 |

### n=100, theoretical censoring=10%

| Method | Estimand | N used | Relative N | Bias | Relative bias | RMSE | SD estimates | SD error | Median error | Error Q2.5% | Error Q25% | Error Q75% | Error Q97.5% |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| complete_case | cv | 2000 | 2000 | -0.18180 | -0.13869 | 0.32438 | 0.26871 | 0.26871 | -0.23806 | -0.51446 | -0.35669 | -0.08320 | 0.46553 |
| complete_case | mean | 2000 | 2000 | 0.15258 | 0.09255 | 0.28007 | 0.23491 | 0.23491 | 0.12841 | -0.23999 | -0.00478 | 0.28414 | 0.68574 |
| complete_case | median | 2000 | 2000 | 0.14015 | 0.14015 | 0.19585 | 0.13684 | 0.13684 | 0.13137 | -0.09833 | 0.04387 | 0.22564 | 0.42516 |
| complete_case | sd | 2000 | 2000 | -0.09458 | -0.04376 | 0.72751 | 0.72152 | 0.72152 | -0.24311 | -0.95582 | -0.55184 | 0.17814 | 1.64619 |
| mle_lognormal | cv | 2000 | 2000 | -0.00323 | -0.00246 | 0.16299 | 0.16300 | 0.16300 | -0.01572 | -0.29271 | -0.12044 | 0.09913 | 0.34493 |
| mle_lognormal | mean | 2000 | 2000 | -0.00255 | -0.00155 | 0.20853 | 0.20856 | 0.20856 | -0.02344 | -0.35585 | -0.14440 | 0.11813 | 0.47501 |
| mle_lognormal | median | 2000 | 2000 | -0.00163 | -0.00163 | 0.09975 | 0.09976 | 0.09976 | -0.00712 | -0.17929 | -0.07230 | 0.06272 | 0.20663 |
| mle_lognormal | sd | 2000 | 2000 | 0.01203 | 0.00556 | 0.49894 | 0.49892 | 0.49892 | -0.06275 | -0.80263 | -0.33147 | 0.29781 | 1.20269 |
| reverse_km | mean | 0 | 0 | — | — | — | — | — | — | — | — | — | — |
| reverse_km | q0.025 | 0 | 0 | — | — | — | — | — | — | — | — | — | — |
| reverse_km | q0.1 | 0 | 0 | — | — | — | — | — | — | — | — | — | — |
| reverse_km | q0.5 | 2000 | 2000 | -0.00857 | -0.00857 | 0.12295 | 0.12268 | 0.12268 | -0.01575 | -0.22006 | -0.09256 | 0.06693 | 0.25409 |
| reverse_km | q0.9 | 2000 | 2000 | -0.08292 | -0.02302 | 0.61444 | 0.60897 | 0.60897 | -0.15882 | -1.07239 | -0.50556 | 0.27466 | 1.35106 |
| reverse_km | q0.975 | 2000 | 2000 | 0.07529 | 0.01061 | 2.02632 | 2.02543 | 2.02543 | -0.28606 | -2.78904 | -1.33622 | 1.03393 | 4.88101 |
| ros_log | cv | 2000 | 2000 | -0.09506 | -0.07252 | 0.29347 | 0.27772 | 0.27772 | -0.15542 | -0.43249 | -0.27566 | 0.01260 | 0.57123 |
| ros_log | mean | 2000 | 2000 | -0.01125 | -0.00682 | 0.21440 | 0.21416 | 0.21416 | -0.03234 | -0.36556 | -0.16085 | 0.11347 | 0.48765 |
| ros_log | median | 2000 | 2000 | 0.00407 | 0.00407 | 0.12409 | 0.12406 | 0.12406 | -0.00248 | -0.21266 | -0.08316 | 0.07882 | 0.27291 |
| ros_log | sd | 2000 | 2000 | -0.14106 | -0.06527 | 0.69502 | 0.68072 | 0.68072 | -0.28190 | -0.95665 | -0.57529 | 0.12334 | 1.46177 |
| substitute_loq | cv | 2000 | 2000 | -0.10648 | -0.08123 | 0.29564 | 0.27587 | 0.27587 | -0.16585 | -0.44226 | -0.28562 | -0.00246 | 0.55916 |
| substitute_loq | mean | 2000 | 2000 | -0.00186 | -0.00113 | 0.21424 | 0.21428 | 0.21428 | -0.02339 | -0.35762 | -0.15261 | 0.12183 | 0.49796 |
| substitute_loq | median | 2000 | 2000 | 0.00407 | 0.00407 | 0.12409 | 0.12406 | 0.12406 | -0.00248 | -0.21266 | -0.08316 | 0.07882 | 0.27291 |
| substitute_loq | sd | 2000 | 2000 | -0.14790 | -0.06843 | 0.69674 | 0.68103 | 0.68103 | -0.28812 | -0.96188 | -0.58291 | 0.11620 | 1.45431 |
| substitute_loq_half | cv | 2000 | 2000 | -0.08931 | -0.06813 | 0.29049 | 0.27649 | 0.27649 | -0.14765 | -0.42515 | -0.26969 | 0.01662 | 0.57520 |
| substitute_loq_half | mean | 2000 | 2000 | -0.01591 | -0.00965 | 0.21556 | 0.21503 | 0.21503 | -0.03643 | -0.37235 | -0.16660 | 0.10943 | 0.48348 |
| substitute_loq_half | median | 2000 | 2000 | 0.00407 | 0.00407 | 0.12409 | 0.12406 | 0.12406 | -0.00248 | -0.21266 | -0.08316 | 0.07882 | 0.27291 |
| substitute_loq_half | sd | 2000 | 2000 | -0.13739 | -0.06357 | 0.69308 | 0.67949 | 0.67949 | -0.27832 | -0.95014 | -0.56934 | 0.12537 | 1.46279 |
| substitute_zero | cv | 2000 | 2000 | -0.07123 | -0.05434 | 0.28605 | 0.27711 | 0.27711 | -0.13014 | -0.41233 | -0.25092 | 0.03685 | 0.60076 |
| substitute_zero | mean | 2000 | 2000 | -0.02995 | -0.01817 | 0.21786 | 0.21585 | 0.21585 | -0.05040 | -0.38535 | -0.18067 | 0.09597 | 0.46845 |
| substitute_zero | median | 2000 | 2000 | 0.00407 | 0.00407 | 0.12409 | 0.12406 | 0.12406 | -0.00248 | -0.21266 | -0.08316 | 0.07882 | 0.27291 |
| substitute_zero | sd | 2000 | 2000 | -0.12603 | -0.05831 | 0.68922 | 0.67777 | 0.67777 | -0.26649 | -0.94012 | -0.55694 | 0.13454 | 1.47187 |

### n=500, theoretical censoring=10%

| Method | Estimand | N used | Relative N | Bias | Relative bias | RMSE | SD estimates | SD error | Median error | Error Q2.5% | Error Q25% | Error Q75% | Error Q97.5% |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| complete_case | cv | 2000 | 2000 | -0.11313 | -0.08630 | 0.21961 | 0.18827 | 0.18827 | -0.15063 | -0.34392 | -0.22845 | -0.04321 | 0.31819 |
| complete_case | mean | 2000 | 2000 | 0.16055 | 0.09738 | 0.19008 | 0.10179 | 0.10179 | 0.15550 | -0.02203 | 0.08994 | 0.22452 | 0.37746 |
| complete_case | median | 2000 | 2000 | 0.13522 | 0.13522 | 0.14799 | 0.06017 | 0.06017 | 0.13237 | 0.02467 | 0.09427 | 0.17439 | 0.25837 |
| complete_case | sd | 2000 | 2000 | 0.01497 | 0.00693 | 0.42822 | 0.42807 | 0.42807 | -0.06214 | -0.53221 | -0.25399 | 0.19469 | 1.00768 |
| mle_lognormal | cv | 2000 | 2000 | -0.00126 | -0.00096 | 0.07235 | 0.07236 | 0.07236 | -0.00378 | -0.13306 | -0.05384 | 0.04780 | 0.13798 |
| mle_lognormal | mean | 2000 | 2000 | -0.00007 | -0.00004 | 0.08849 | 0.08851 | 0.08851 | -0.00380 | -0.16075 | -0.05817 | 0.05378 | 0.18904 |
| mle_lognormal | median | 2000 | 2000 | 0.00037 | 0.00037 | 0.04442 | 0.04443 | 0.04443 | 0.00043 | -0.08335 | -0.03089 | 0.02861 | 0.08978 |
| mle_lognormal | sd | 2000 | 2000 | 0.00145 | 0.00067 | 0.20899 | 0.20904 | 0.20904 | -0.01335 | -0.36339 | -0.14609 | 0.13582 | 0.46008 |
| reverse_km | mean | 0 | 0 | — | — | — | — | — | — | — | — | — | — |
| reverse_km | q0.025 | 0 | 0 | — | — | — | — | — | — | — | — | — | — |
| reverse_km | q0.1 | 0 | 0 | — | — | — | — | — | — | — | — | — | — |
| reverse_km | q0.5 | 2000 | 2000 | -0.00034 | -0.00034 | 0.05551 | 0.05552 | 0.05552 | -0.00301 | -0.10223 | -0.03811 | 0.03730 | 0.11215 |
| reverse_km | q0.9 | 2000 | 2000 | -0.01424 | -0.00395 | 0.27512 | 0.27482 | 0.27482 | -0.02495 | -0.51800 | -0.20253 | 0.15602 | 0.60769 |
| reverse_km | q0.975 | 2000 | 2000 | 0.02454 | 0.00346 | 0.83232 | 0.83216 | 0.83216 | -0.03491 | -1.41349 | -0.56696 | 0.51370 | 1.83949 |
| ros_log | cv | 2000 | 2000 | -0.02737 | -0.02088 | 0.19502 | 0.19314 | 0.19314 | -0.06607 | -0.26591 | -0.14690 | 0.04289 | 0.40697 |
| ros_log | mean | 2000 | 2000 | -0.00174 | -0.00106 | 0.09317 | 0.09318 | 0.09318 | -0.00665 | -0.16776 | -0.06512 | 0.05378 | 0.19825 |
| ros_log | median | 2000 | 2000 | 0.00216 | 0.00216 | 0.05566 | 0.05563 | 0.05563 | -0.00045 | -0.09969 | -0.03570 | 0.04051 | 0.11473 |
| ros_log | sd | 2000 | 2000 | -0.03900 | -0.01804 | 0.40438 | 0.40260 | 0.40260 | -0.11035 | -0.55513 | -0.29091 | 0.13327 | 0.88065 |
| substitute_loq | cv | 2000 | 2000 | -0.03848 | -0.02935 | 0.19607 | 0.19231 | 0.19231 | -0.07718 | -0.27532 | -0.15721 | 0.03095 | 0.39207 |
| substitute_loq | mean | 2000 | 2000 | 0.00752 | 0.00456 | 0.09347 | 0.09319 | 0.09319 | 0.00248 | -0.15868 | -0.05610 | 0.06292 | 0.20840 |
| substitute_loq | median | 2000 | 2000 | 0.00216 | 0.00216 | 0.05566 | 0.05563 | 0.05563 | -0.00045 | -0.09969 | -0.03570 | 0.04051 | 0.11473 |
| substitute_loq | sd | 2000 | 2000 | -0.04539 | -0.02100 | 0.40551 | 0.40306 | 0.40306 | -0.11630 | -0.56247 | -0.29790 | 0.12678 | 0.87642 |
| substitute_loq_half | cv | 2000 | 2000 | -0.02182 | -0.01664 | 0.19405 | 0.19286 | 0.19286 | -0.06074 | -0.25778 | -0.14154 | 0.04801 | 0.41210 |
| substitute_loq_half | mean | 2000 | 2000 | -0.00634 | -0.00384 | 0.09375 | 0.09356 | 0.09356 | -0.01106 | -0.17212 | -0.06982 | 0.04940 | 0.19490 |
| substitute_loq_half | median | 2000 | 2000 | 0.00216 | 0.00216 | 0.05566 | 0.05563 | 0.05563 | -0.00045 | -0.09969 | -0.03570 | 0.04051 | 0.11473 |
| substitute_loq_half | sd | 2000 | 2000 | -0.03577 | -0.01655 | 0.40348 | 0.40199 | 0.40199 | -0.10756 | -0.55207 | -0.28780 | 0.13601 | 0.88407 |
| substitute_zero | cv | 2000 | 2000 | -0.00437 | -0.00334 | 0.19339 | 0.19339 | 0.19339 | -0.04468 | -0.24025 | -0.12407 | 0.06519 | 0.42800 |
| substitute_zero | mean | 2000 | 2000 | -0.02020 | -0.01225 | 0.09609 | 0.09396 | 0.09396 | -0.02422 | -0.18805 | -0.08354 | 0.03665 | 0.18297 |
| substitute_zero | median | 2000 | 2000 | 0.00216 | 0.00216 | 0.05566 | 0.05563 | 0.05563 | -0.00045 | -0.09969 | -0.03570 | 0.04051 | 0.11473 |
| substitute_zero | sd | 2000 | 2000 | -0.02536 | -0.01173 | 0.40152 | 0.40082 | 0.40082 | -0.09754 | -0.54067 | -0.27668 | 0.14581 | 0.89233 |

### n=30, theoretical censoring=40%

| Method | Estimand | N used | Relative N | Bias | Relative bias | RMSE | SD estimates | SD error | Median error | Error Q2.5% | Error Q25% | Error Q75% | Error Q97.5% |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| complete_case | cv | 2000 | 2000 | -0.49899 | -0.38067 | 0.57405 | 0.28387 | 0.28387 | -0.55415 | -0.88807 | -0.68866 | -0.37772 | 0.20165 |
| complete_case | mean | 2000 | 2000 | 0.80721 | 0.48959 | 0.99871 | 0.58823 | 0.58823 | 0.73448 | -0.02927 | 0.39074 | 1.11663 | 2.16446 |
| complete_case | median | 2000 | 2000 | 0.74609 | 0.74609 | 0.83460 | 0.37413 | 0.37413 | 0.68379 | 0.19638 | 0.48875 | 0.94524 | 1.56500 |
| complete_case | sd | 2000 | 2000 | -0.05857 | -0.02710 | 1.28889 | 1.28788 | 1.28788 | -0.34179 | -1.42268 | -0.84765 | 0.33250 | 2.97570 |
| mle_lognormal | cv | 2000 | 2000 | 0.01657 | 0.01264 | 0.46048 | 0.46029 | 0.46029 | -0.07495 | -0.57735 | -0.28895 | 0.20799 | 1.17378 |
| mle_lognormal | mean | 2000 | 2000 | 0.03153 | 0.01912 | 0.39349 | 0.39233 | 0.39233 | -0.01549 | -0.58227 | -0.24448 | 0.25550 | 0.89390 |
| mle_lognormal | median | 2000 | 2000 | 0.01817 | 0.01817 | 0.21516 | 0.21445 | 0.21445 | 0.00448 | -0.37273 | -0.12543 | 0.16278 | 0.46640 |
| mle_lognormal | sd | 2000 | 2000 | 0.15854 | 0.07336 | 1.31598 | 1.30672 | 1.30672 | -0.13041 | -1.25495 | -0.67217 | 0.61547 | 3.31852 |
| reverse_km | mean | 0 | 0 | — | — | — | — | — | — | — | — | — | — |
| reverse_km | q0.025 | 0 | 0 | — | — | — | — | — | — | — | — | — | — |
| reverse_km | q0.1 | 0 | 0 | — | — | — | — | — | — | — | — | — | — |
| reverse_km | q0.5 | 1654 | 1654 | 0.05453 | 0.05453 | 0.21125 | 0.20415 | 0.20415 | 0.01236 | -0.20620 | -0.10004 | 0.16916 | 0.54854 |
| reverse_km | q0.9 | 2000 | 2000 | -0.22097 | -0.06134 | 1.04802 | 1.02472 | 1.02472 | -0.40950 | -1.73807 | -0.96349 | 0.38399 | 1.97466 |
| reverse_km | q0.975 | 0 | 0 | — | — | — | — | — | — | — | — | — | — |
| ros_log | cv | 2000 | 2000 | -0.17881 | -0.13641 | 0.39994 | 0.35783 | 0.35783 | -0.24609 | -0.67949 | -0.42116 | -0.02773 | 0.70958 |
| ros_log | mean | 2000 | 2000 | -0.00699 | -0.00424 | 0.38471 | 0.38474 | 0.38474 | -0.04908 | -0.60280 | -0.27392 | 0.21231 | 0.81480 |
| ros_log | median | 2000 | 2000 | 0.02980 | 0.02980 | 0.25255 | 0.25085 | 0.25085 | 0.00654 | -0.45671 | -0.12857 | 0.18276 | 0.57581 |
| ros_log | sd | 2000 | 2000 | -0.24732 | -0.11444 | 1.03126 | 1.00142 | 1.00142 | -0.45115 | -1.36395 | -0.86123 | 0.08231 | 2.19002 |
| substitute_loq | cv | 2000 | 2000 | -0.33731 | -0.25732 | 0.46333 | 0.31772 | 0.31772 | -0.39559 | -0.76707 | -0.54966 | -0.20248 | 0.43833 |
| substitute_loq | mean | 2000 | 2000 | 0.13916 | 0.08440 | 0.40077 | 0.37593 | 0.37593 | 0.09530 | -0.43077 | -0.12349 | 0.33979 | 0.98181 |
| substitute_loq | median | 2000 | 2000 | 0.04595 | 0.04595 | 0.22787 | 0.22324 | 0.22324 | 0.00512 | -0.22380 | -0.13502 | 0.18276 | 0.57581 |
| substitute_loq | sd | 2000 | 2000 | -0.34216 | -0.15832 | 1.06046 | 1.00400 | 1.00400 | -0.54856 | -1.46028 | -0.95660 | -0.02042 | 2.11909 |
| substitute_loq_half | cv | 2000 | 2000 | -0.17448 | -0.13311 | 0.36554 | 0.32128 | 0.32128 | -0.23992 | -0.59647 | -0.39074 | -0.04347 | 0.61797 |
| substitute_loq_half | mean | 2000 | 2000 | -0.01487 | -0.00902 | 0.39141 | 0.39122 | 0.39122 | -0.05800 | -0.62336 | -0.28892 | 0.20387 | 0.83697 |
| substitute_loq_half | median | 2000 | 2000 | -0.00538 | -0.00538 | 0.30058 | 0.30060 | 0.30060 | 0.00512 | -0.61190 | -0.14244 | 0.18276 | 0.57581 |
| substitute_loq_half | sd | 2000 | 2000 | -0.24131 | -0.11165 | 1.01796 | 0.98919 | 0.98919 | -0.44537 | -1.32327 | -0.85180 | 0.07268 | 2.17773 |
| substitute_zero | cv | 2000 | 2000 | 0.04291 | 0.03274 | 0.34211 | 0.33950 | 0.33950 | -0.02145 | -0.42579 | -0.19067 | 0.18510 | 0.89912 |
| substitute_zero | mean | 2000 | 2000 | -0.16890 | -0.10244 | 0.44250 | 0.40910 | 0.40910 | -0.21213 | -0.81477 | -0.45051 | 0.06365 | 0.69057 |
| substitute_zero | median | 2000 | 2000 | -0.05670 | -0.05670 | 0.40308 | 0.39917 | 0.39917 | 0.00512 | -1.00000 | -0.14244 | 0.18276 | 0.57581 |
| substitute_zero | sd | 2000 | 2000 | -0.12661 | -0.05858 | 0.97978 | 0.97181 | 0.97181 | -0.33052 | -1.16410 | -0.72606 | 0.18068 | 2.24015 |

### n=100, theoretical censoring=40%

| Method | Estimand | N used | Relative N | Bias | Relative bias | RMSE | SD estimates | SD error | Median error | Error Q2.5% | Error Q25% | Error Q75% | Error Q97.5% |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| complete_case | cv | 2000 | 2000 | -0.38996 | -0.29749 | 0.46884 | 0.26034 | 0.26034 | -0.44144 | -0.71237 | -0.54744 | -0.28718 | 0.23506 |
| complete_case | mean | 2000 | 2000 | 0.82632 | 0.50119 | 0.88768 | 0.32437 | 0.32437 | 0.79954 | 0.29125 | 0.60189 | 1.01903 | 1.51038 |
| complete_case | median | 2000 | 2000 | 0.71152 | 0.71152 | 0.73787 | 0.19550 | 0.19550 | 0.69105 | 0.37580 | 0.57792 | 0.82722 | 1.17234 |
| complete_case | sd | 2000 | 2000 | 0.17058 | 0.07893 | 1.00357 | 0.98922 | 0.98922 | -0.01279 | -0.94752 | -0.41005 | 0.49954 | 2.30889 |
| mle_lognormal | cv | 2000 | 2000 | 0.01756 | 0.01340 | 0.21943 | 0.21878 | 0.21878 | -0.00809 | -0.33327 | -0.13729 | 0.14664 | 0.50823 |
| mle_lognormal | mean | 2000 | 2000 | 0.01361 | 0.00826 | 0.20986 | 0.20947 | 0.20947 | 0.00196 | -0.36368 | -0.13025 | 0.14789 | 0.44853 |
| mle_lognormal | median | 2000 | 2000 | 0.00090 | 0.00090 | 0.11447 | 0.11450 | 0.11450 | 0.00041 | -0.22008 | -0.07669 | 0.07778 | 0.22804 |
| mle_lognormal | sd | 2000 | 2000 | 0.06993 | 0.03236 | 0.57374 | 0.56960 | 0.56960 | -0.01727 | -0.81921 | -0.34072 | 0.40562 | 1.33964 |
| reverse_km | mean | 0 | 0 | — | — | — | — | — | — | — | — | — | — |
| reverse_km | q0.025 | 0 | 0 | — | — | — | — | — | — | — | — | — | — |
| reverse_km | q0.1 | 0 | 0 | — | — | — | — | — | — | — | — | — | — |
| reverse_km | q0.5 | 1931 | 1931 | 0.00066 | 0.00066 | 0.11809 | 0.11812 | 0.11812 | -0.00926 | -0.19744 | -0.08802 | 0.07481 | 0.25549 |
| reverse_km | q0.9 | 2000 | 2000 | -0.04918 | -0.01365 | 0.60222 | 0.60036 | 0.60036 | -0.09168 | -1.09456 | -0.46399 | 0.30543 | 1.24732 |
| reverse_km | q0.975 | 2000 | 2000 | 0.22363 | 0.03150 | 2.06476 | 2.05313 | 2.05313 | -0.15617 | -2.77148 | -1.24755 | 1.33033 | 5.23397 |
| ros_log | cv | 2000 | 2000 | -0.07155 | -0.05458 | 0.31684 | 0.30873 | 0.30873 | -0.13015 | -0.46188 | -0.25956 | 0.04503 | 0.66350 |
| ros_log | mean | 2000 | 2000 | 0.00064 | 0.00039 | 0.21396 | 0.21401 | 0.21401 | -0.00928 | -0.37907 | -0.14526 | 0.12786 | 0.46083 |
| ros_log | median | 2000 | 2000 | 0.00481 | 0.00481 | 0.12706 | 0.12700 | 0.12700 | -0.00156 | -0.21711 | -0.08279 | 0.08505 | 0.26521 |
| ros_log | sd | 2000 | 2000 | -0.08732 | -0.04040 | 0.75632 | 0.75145 | 0.75145 | -0.21886 | -0.98088 | -0.54081 | 0.17076 | 1.50439 |
| substitute_loq | cv | 2000 | 2000 | -0.22261 | -0.16982 | 0.36189 | 0.28539 | 0.28539 | -0.27965 | -0.57713 | -0.39274 | -0.11568 | 0.44764 |
| substitute_loq | mean | 2000 | 2000 | 0.14431 | 0.08753 | 0.25515 | 0.21048 | 0.21048 | 0.12893 | -0.21351 | 0.00254 | 0.26938 | 0.59060 |
| substitute_loq | median | 2000 | 2000 | 0.00568 | 0.00568 | 0.12469 | 0.12459 | 0.12459 | -0.00156 | -0.21988 | -0.08292 | 0.08505 | 0.26521 |
| substitute_loq | sd | 2000 | 2000 | -0.17318 | -0.08013 | 0.77683 | 0.75747 | 0.75747 | -0.30558 | -1.07148 | -0.62762 | 0.08978 | 1.43968 |
| substitute_loq_half | cv | 2000 | 2000 | -0.06074 | -0.04634 | 0.29927 | 0.29311 | 0.29311 | -0.12242 | -0.40993 | -0.23657 | 0.04393 | 0.66050 |
| substitute_loq_half | mean | 2000 | 2000 | -0.01150 | -0.00698 | 0.21912 | 0.21887 | 0.21887 | -0.02310 | -0.39411 | -0.16331 | 0.11924 | 0.45758 |
| substitute_loq_half | median | 2000 | 2000 | -0.00499 | -0.00499 | 0.15452 | 0.15447 | 0.15447 | -0.00156 | -0.41277 | -0.08292 | 0.08505 | 0.26521 |
| substitute_loq_half | sd | 2000 | 2000 | -0.08189 | -0.03789 | 0.74957 | 0.74527 | 0.74527 | -0.21260 | -0.95114 | -0.52524 | 0.17053 | 1.49633 |
| substitute_zero | cv | 2000 | 2000 | 0.14736 | 0.11242 | 0.33814 | 0.30442 | 0.30442 | 0.08444 | -0.21687 | -0.03341 | 0.25130 | 0.90757 |
| substitute_zero | mean | 2000 | 2000 | -0.16731 | -0.10148 | 0.28322 | 0.22858 | 0.22858 | -0.17872 | -0.57616 | -0.32500 | -0.03153 | 0.31850 |
| substitute_zero | median | 2000 | 2000 | -0.01566 | -0.01566 | 0.19858 | 0.19801 | 0.19801 | -0.00156 | -0.60682 | -0.08292 | 0.08505 | 0.26521 |
| substitute_zero | sd | 2000 | 2000 | 0.02279 | 0.01055 | 0.73126 | 0.73109 | 0.73109 | -0.10481 | -0.81906 | -0.41047 | 0.26890 | 1.57465 |

### n=500, theoretical censoring=40%

| Method | Estimand | N used | Relative N | Bias | Relative bias | RMSE | SD estimates | SD error | Median error | Error Q2.5% | Error Q25% | Error Q75% | Error Q97.5% |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| complete_case | cv | 2000 | 2000 | -0.32881 | -0.25084 | 0.36872 | 0.16689 | 0.16689 | -0.36440 | -0.54888 | -0.43910 | -0.25842 | 0.08691 |
| complete_case | mean | 2000 | 2000 | 0.80962 | 0.49106 | 0.82192 | 0.14172 | 0.14172 | 0.80465 | 0.54866 | 0.71130 | 0.89929 | 1.09705 |
| complete_case | median | 2000 | 2000 | 0.69076 | 0.69076 | 0.69588 | 0.08429 | 0.08429 | 0.68505 | 0.53562 | 0.63299 | 0.74784 | 0.85837 |
| complete_case | sd | 2000 | 2000 | 0.26638 | 0.12325 | 0.58154 | 0.51707 | 0.51707 | 0.17578 | -0.44738 | -0.07969 | 0.49583 | 1.50591 |
| mle_lognormal | cv | 2000 | 2000 | 0.00013 | 0.00010 | 0.09319 | 0.09322 | 0.09322 | -0.00605 | -0.16645 | -0.06608 | 0.06157 | 0.19150 |
| mle_lognormal | mean | 2000 | 2000 | -0.00030 | -0.00018 | 0.08969 | 0.08971 | 0.08971 | -0.00229 | -0.16931 | -0.06351 | 0.06056 | 0.17788 |
| mle_lognormal | median | 2000 | 2000 | -0.00008 | -0.00008 | 0.04946 | 0.04947 | 0.04947 | 0.00030 | -0.09658 | -0.03325 | 0.03218 | 0.09988 |
| mle_lognormal | sd | 2000 | 2000 | 0.00413 | 0.00191 | 0.23769 | 0.23772 | 0.23772 | -0.00471 | -0.42555 | -0.16163 | 0.16156 | 0.50275 |
| reverse_km | mean | 0 | 0 | — | — | — | — | — | — | — | — | — | — |
| reverse_km | q0.025 | 0 | 0 | — | — | — | — | — | — | — | — | — | — |
| reverse_km | q0.1 | 0 | 0 | — | — | — | — | — | — | — | — | — | — |
| reverse_km | q0.5 | 2000 | 2000 | -0.00110 | -0.00110 | 0.05569 | 0.05569 | 0.05569 | -0.00254 | -0.10947 | -0.03919 | 0.03579 | 0.11191 |
| reverse_km | q0.9 | 2000 | 2000 | -0.02421 | -0.00672 | 0.26838 | 0.26735 | 0.26735 | -0.03718 | -0.52212 | -0.21188 | 0.15615 | 0.50272 |
| reverse_km | q0.975 | 2000 | 2000 | 0.03284 | 0.00463 | 0.85364 | 0.85322 | 0.85322 | -0.03487 | -1.44184 | -0.54436 | 0.58645 | 1.83164 |
| ros_log | cv | 2000 | 2000 | -0.01955 | -0.01492 | 0.18909 | 0.18812 | 0.18812 | -0.05503 | -0.27626 | -0.14439 | 0.06079 | 0.43910 |
| ros_log | mean | 2000 | 2000 | -0.00282 | -0.00171 | 0.09207 | 0.09205 | 0.09205 | -0.00779 | -0.17580 | -0.06832 | 0.05892 | 0.18331 |
| ros_log | median | 2000 | 2000 | 0.00151 | 0.00151 | 0.05571 | 0.05570 | 0.05570 | -0.00035 | -0.10608 | -0.03671 | 0.03814 | 0.11366 |
| ros_log | sd | 2000 | 2000 | -0.02843 | -0.01316 | 0.38627 | 0.38532 | 0.38532 | -0.08962 | -0.57241 | -0.28535 | 0.14453 | 0.93010 |
| substitute_loq | cv | 2000 | 2000 | -0.16505 | -0.12591 | 0.24108 | 0.17578 | 0.17578 | -0.20048 | -0.40270 | -0.28263 | -0.09074 | 0.25960 |
| substitute_loq | mean | 2000 | 2000 | 0.13633 | 0.08269 | 0.16397 | 0.09113 | 0.09113 | 0.13308 | -0.03575 | 0.07194 | 0.19542 | 0.32631 |
| substitute_loq | median | 2000 | 2000 | 0.00151 | 0.00151 | 0.05571 | 0.05570 | 0.05570 | -0.00035 | -0.10608 | -0.03671 | 0.03814 | 0.11366 |
| substitute_loq | sd | 2000 | 2000 | -0.10703 | -0.04952 | 0.40513 | 0.39084 | 0.39084 | -0.16897 | -0.65215 | -0.36724 | 0.06871 | 0.87052 |
| substitute_loq_half | cv | 2000 | 2000 | -0.00326 | -0.00249 | 0.18272 | 0.18274 | 0.18274 | -0.04068 | -0.24331 | -0.12549 | 0.07216 | 0.44734 |
| substitute_loq_half | mean | 2000 | 2000 | -0.01899 | -0.01152 | 0.09638 | 0.09451 | 0.09451 | -0.02226 | -0.19492 | -0.08525 | 0.04359 | 0.17654 |
| substitute_loq_half | median | 2000 | 2000 | 0.00151 | 0.00151 | 0.05571 | 0.05570 | 0.05570 | -0.00035 | -0.10608 | -0.03671 | 0.03814 | 0.11366 |
| substitute_loq_half | sd | 2000 | 2000 | -0.02204 | -0.01020 | 0.38261 | 0.38207 | 0.38207 | -0.08366 | -0.55487 | -0.27655 | 0.14872 | 0.93261 |
| substitute_zero | cv | 2000 | 2000 | 0.20209 | 0.15417 | 0.27804 | 0.19101 | 0.19101 | 0.16333 | -0.04642 | 0.07502 | 0.28466 | 0.67581 |
| substitute_zero | mean | 2000 | 2000 | -0.17431 | -0.10572 | 0.20019 | 0.09849 | 0.09849 | -0.17751 | -0.35755 | -0.24270 | -0.10893 | 0.02445 |
| substitute_zero | median | 2000 | 2000 | 0.00151 | 0.00151 | 0.05571 | 0.05570 | 0.05570 | -0.00035 | -0.10608 | -0.03671 | 0.03814 | 0.11366 |
| substitute_zero | sd | 2000 | 2000 | 0.07610 | 0.03521 | 0.37965 | 0.37204 | 0.37204 | 0.01792 | -0.44086 | -0.17193 | 0.23955 | 1.00483 |

## Sampling variability in the complete samples

The following distributions describe the complete-sample statistics before censoring. Their spread is sampling variability under the generating population; it must not be attributed to a censoring method. `truth_mean - population` and `truth RMSE` summarize those sample statistics against their population targets, not estimator performance after censoring.

| Scenario | Statistic | N | Population target | Mean truth | SD truth | Median truth | Truth Q2.5% | Truth Q25% | Truth Q75% | Truth Q97.5% | Mean truth − population | Truth RMSE vs population |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| n30_LOQ_10 | cv | 2000 | 1.31083 | 1.10994 | 0.30749 | 1.04382 | 0.69887 | 0.90461 | 1.23261 | 1.89521 | -0.20089 | 0.36723 |
| n30_LOQ_10 | mean | 2000 | 1.64872 | 1.65679 | 0.40435 | 1.58130 | 1.05795 | 1.37808 | 1.88944 | 2.63478 | 0.00807 | 0.40433 |
| n30_LOQ_10 | median | 2000 | 1.00000 | 1.03473 | 0.24044 | 1.00573 | 0.64720 | 0.86179 | 1.17190 | 1.59391 | 0.03473 | 0.24287 |
| n30_LOQ_10 | sd | 2000 | 2.16120 | 1.90451 | 0.94860 | 1.66724 | 0.82462 | 1.29534 | 2.21241 | 4.38787 | -0.25669 | 0.98249 |
| n100_LOQ_10 | cv | 2000 | 1.31083 | 1.21561 | 0.27627 | 1.15604 | 0.87983 | 1.03583 | 1.32219 | 1.88117 | -0.09522 | 0.29215 |
| n100_LOQ_10 | mean | 2000 | 1.64872 | 1.63762 | 0.21479 | 1.61706 | 1.28041 | 1.48729 | 1.76227 | 2.13726 | -0.01110 | 0.21502 |
| n100_LOQ_10 | median | 2000 | 1.00000 | 1.00407 | 0.12406 | 0.99752 | 0.78734 | 0.91684 | 1.07882 | 1.27291 | 0.00407 | 0.12409 |
| n100_LOQ_10 | sd | 2000 | 2.16120 | 2.02021 | 0.68002 | 1.87921 | 1.20626 | 1.58735 | 2.28332 | 3.62238 | -0.14099 | 0.69431 |
| n500_LOQ_10 | cv | 2000 | 1.31083 | 1.28336 | 0.19272 | 1.24467 | 1.04809 | 1.16413 | 1.35269 | 1.71416 | -0.02747 | 0.19462 |
| n500_LOQ_10 | mean | 2000 | 1.64872 | 1.64706 | 0.09343 | 1.64211 | 1.48032 | 1.58388 | 1.70288 | 1.84794 | -0.00166 | 0.09342 |
| n500_LOQ_10 | median | 2000 | 1.00000 | 1.00216 | 0.05563 | 0.99955 | 0.90031 | 0.96430 | 1.04051 | 1.11473 | 0.00216 | 0.05566 |
| n500_LOQ_10 | sd | 2000 | 2.16120 | 2.12218 | 0.40237 | 2.05072 | 1.60517 | 1.87031 | 2.29456 | 3.04303 | -0.03902 | 0.40415 |
| n30_LOQ_40 | cv | 2000 | 1.31083 | 1.11969 | 0.32071 | 1.05579 | 0.69742 | 0.90293 | 1.24770 | 1.92020 | -0.19115 | 0.37329 |
| n30_LOQ_40 | mean | 2000 | 1.64872 | 1.65093 | 0.39057 | 1.60664 | 1.03756 | 1.37528 | 1.86758 | 2.48351 | 0.00220 | 0.39048 |
| n30_LOQ_40 | median | 2000 | 1.00000 | 1.03420 | 0.23974 | 1.00512 | 0.63866 | 0.86279 | 1.18276 | 1.57581 | 0.03420 | 0.24210 |
| n30_LOQ_40 | sd | 2000 | 2.16120 | 1.91238 | 0.98956 | 1.71141 | 0.82503 | 1.30396 | 2.22833 | 4.32975 | -0.24881 | 1.02012 |
| n100_LOQ_40 | cv | 2000 | 1.31083 | 1.23161 | 0.29228 | 1.16965 | 0.88041 | 1.05613 | 1.33900 | 1.94508 | -0.07923 | 0.30276 |
| n100_LOQ_40 | mean | 2000 | 1.64872 | 1.65542 | 0.21800 | 1.64260 | 1.27442 | 1.50325 | 1.78569 | 2.12027 | 0.00670 | 0.21804 |
| n100_LOQ_40 | median | 2000 | 1.00000 | 1.00479 | 0.12640 | 0.99844 | 0.77569 | 0.91708 | 1.08505 | 1.26521 | 0.00479 | 0.12646 |
| n100_LOQ_40 | sd | 2000 | 2.16120 | 2.07171 | 0.74609 | 1.94177 | 1.20027 | 1.62621 | 2.32708 | 3.65475 | -0.08948 | 0.75125 |
| n500_LOQ_40 | cv | 2000 | 1.31083 | 1.28921 | 0.18195 | 1.25270 | 1.04995 | 1.16655 | 1.36322 | 1.72942 | -0.02162 | 0.18318 |
| n500_LOQ_40 | mean | 2000 | 1.64872 | 1.64755 | 0.09416 | 1.64440 | 1.47298 | 1.58106 | 1.70953 | 1.84273 | -0.00117 | 0.09414 |
| n500_LOQ_40 | median | 2000 | 1.00000 | 1.00151 | 0.05570 | 0.99965 | 0.89392 | 0.96329 | 1.03814 | 1.11366 | 0.00151 | 0.05571 |
| n500_LOQ_40 | sd | 2000 | 2.16120 | 2.13221 | 0.38271 | 2.07033 | 1.59819 | 1.87700 | 2.30295 | 3.08620 | -0.02898 | 0.38371 |
| n30_LOQ_10 | q0.025 | 2000 | 0.14086 | 0.19367 | 0.06843 | 0.18545 | 0.08324 | 0.14409 | 0.23517 | 0.34549 | 0.05281 | 0.08643 |
| n30_LOQ_10 | q0.1 | 2000 | 0.27761 | 0.31594 | 0.09083 | 0.30318 | 0.16616 | 0.25052 | 0.37249 | 0.52733 | 0.03833 | 0.09857 |
| n30_LOQ_10 | q0.5 | 2000 | 1.00000 | 1.03473 | 0.24044 | 1.00573 | 0.64720 | 0.86179 | 1.17190 | 1.59391 | 0.03473 | 0.24287 |
| n30_LOQ_10 | q0.9 | 2000 | 3.60222 | 3.51591 | 1.12379 | 3.31955 | 1.95852 | 2.73110 | 4.05025 | 6.26547 | -0.08632 | 1.12682 |
| n30_LOQ_10 | q0.975 | 2000 | 7.09907 | 6.34445 | 2.82559 | 5.68255 | 2.89852 | 4.46673 | 7.45092 | 13.46868 | -0.75462 | 2.92394 |
| n100_LOQ_10 | q0.025 | 2000 | 0.14086 | 0.15759 | 0.03679 | 0.15473 | 0.09311 | 0.13138 | 0.18132 | 0.23508 | 0.01673 | 0.04041 |
| n100_LOQ_10 | q0.1 | 2000 | 0.27761 | 0.28709 | 0.04628 | 0.28474 | 0.20875 | 0.25401 | 0.31713 | 0.38961 | 0.00949 | 0.04723 |
| n100_LOQ_10 | q0.5 | 2000 | 1.00000 | 1.00407 | 0.12406 | 0.99752 | 0.78734 | 0.91684 | 1.07882 | 1.27291 | 0.00407 | 0.12409 |
| n100_LOQ_10 | q0.9 | 2000 | 3.60222 | 3.54064 | 0.61077 | 3.46438 | 2.53891 | 3.11843 | 3.89732 | 4.97773 | -0.06158 | 0.61371 |
| n100_LOQ_10 | q0.975 | 2000 | 7.09907 | 6.69729 | 1.72008 | 6.43748 | 4.20927 | 5.50494 | 7.56422 | 10.74500 | -0.40178 | 1.76596 |
| n500_LOQ_10 | q0.025 | 2000 | 0.14086 | 0.14488 | 0.01671 | 0.14454 | 0.11404 | 0.13348 | 0.15627 | 0.17902 | 0.00402 | 0.01718 |
| n500_LOQ_10 | q0.1 | 2000 | 0.27761 | 0.28018 | 0.02145 | 0.27986 | 0.23897 | 0.26535 | 0.29492 | 0.32265 | 0.00258 | 0.02160 |
| n500_LOQ_10 | q0.5 | 2000 | 1.00000 | 1.00216 | 0.05563 | 0.99955 | 0.90031 | 0.96430 | 1.04051 | 1.11473 | 0.00216 | 0.05566 |
| n500_LOQ_10 | q0.9 | 2000 | 3.60222 | 3.59213 | 0.27496 | 3.58222 | 3.09059 | 3.40172 | 3.76191 | 4.21110 | -0.01009 | 0.27507 |
| n500_LOQ_10 | q0.975 | 2000 | 7.09907 | 7.00964 | 0.79937 | 6.95795 | 5.63841 | 6.45538 | 7.48961 | 8.71493 | -0.08944 | 0.80416 |
| n30_LOQ_40 | q0.025 | 2000 | 0.14086 | 0.19469 | 0.07043 | 0.18779 | 0.08284 | 0.14259 | 0.23723 | 0.35858 | 0.05382 | 0.08863 |
| n30_LOQ_40 | q0.1 | 2000 | 0.27761 | 0.31548 | 0.09521 | 0.30451 | 0.16182 | 0.24716 | 0.36983 | 0.54002 | 0.03788 | 0.10244 |
| n30_LOQ_40 | q0.5 | 2000 | 1.00000 | 1.03420 | 0.23974 | 1.00512 | 0.63866 | 0.86279 | 1.18276 | 1.57581 | 0.03420 | 0.24210 |
| n30_LOQ_40 | q0.9 | 2000 | 3.60222 | 3.45579 | 1.04011 | 3.27631 | 1.92384 | 2.70492 | 4.07261 | 5.72474 | -0.14643 | 1.05011 |
| n30_LOQ_40 | q0.975 | 2000 | 7.09907 | 6.34774 | 2.70103 | 5.83386 | 2.95285 | 4.54179 | 7.45754 | 12.86262 | -0.75133 | 2.80292 |
| n100_LOQ_40 | q0.025 | 2000 | 0.14086 | 0.15821 | 0.03691 | 0.15705 | 0.09148 | 0.13185 | 0.18169 | 0.23775 | 0.01735 | 0.04078 |
| n100_LOQ_40 | q0.1 | 2000 | 0.27761 | 0.28914 | 0.04866 | 0.28569 | 0.20556 | 0.25496 | 0.31994 | 0.39640 | 0.01154 | 0.05000 |
| n100_LOQ_40 | q0.5 | 2000 | 1.00000 | 1.00479 | 0.12640 | 0.99844 | 0.77569 | 0.91708 | 1.08505 | 1.26521 | 0.00479 | 0.12646 |
| n100_LOQ_40 | q0.9 | 2000 | 3.60222 | 3.57383 | 0.60274 | 3.53467 | 2.52235 | 3.15753 | 3.93494 | 4.87829 | -0.02839 | 0.60326 |
| n100_LOQ_40 | q0.975 | 2000 | 7.09907 | 6.84624 | 1.78059 | 6.58408 | 4.17065 | 5.59899 | 7.76394 | 11.21585 | -0.25283 | 1.79801 |
| n500_LOQ_40 | q0.025 | 2000 | 0.14086 | 0.14447 | 0.01657 | 0.14362 | 0.11328 | 0.13334 | 0.15520 | 0.17884 | 0.00360 | 0.01695 |
| n500_LOQ_40 | q0.1 | 2000 | 0.27761 | 0.27916 | 0.02102 | 0.27854 | 0.23865 | 0.26484 | 0.29276 | 0.32109 | 0.00156 | 0.02107 |
| n500_LOQ_40 | q0.5 | 2000 | 1.00000 | 1.00151 | 0.05570 | 0.99965 | 0.89392 | 0.96329 | 1.03814 | 1.11366 | 0.00151 | 0.05571 |
| n500_LOQ_40 | q0.9 | 2000 | 3.60222 | 3.58215 | 0.26730 | 3.56865 | 3.08620 | 3.39535 | 3.76270 | 4.10670 | -0.02008 | 0.26799 |
| n500_LOQ_40 | q0.975 | 2000 | 7.09907 | 7.01819 | 0.81674 | 6.96182 | 5.61737 | 6.46054 | 7.52175 | 8.75047 | -0.08088 | 0.82053 |

## Empirical censoring and error association

Censoring summaries are based on one record per scenario/replicate. The empirical proportion is a realized sample characteristic, not an assigned treatment. Correlations below are descriptive associations, not causal effects. Pearson and Spearman correlations use only status-eligible finite errors and are reported separately for the sample-truth and population references.

| Scenario | Theory | Fixed LOQ | Replicates | Mean | Median | SD | Min | Q2.5% | Q25% | Q75% | Q97.5% | Max | Distinct values |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| n30_LOQ_10 | 0.10000 | 0.27761 | 2000 | 0.10023 | 0.10000 | 0.05408 | 0.00000 | 0.00000 | 0.06667 | 0.13333 | 0.23333 | 0.30000 | 10 |
| n100_LOQ_10 | 0.10000 | 0.27761 | 2000 | 0.10119 | 0.10000 | 0.02946 | 0.01000 | 0.05000 | 0.08000 | 0.12000 | 0.16000 | 0.20000 | 20 |
| n500_LOQ_10 | 0.10000 | 0.27761 | 2000 | 0.09987 | 0.10000 | 0.01358 | 0.05800 | 0.07595 | 0.09000 | 0.11000 | 0.12800 | 0.14000 | 40 |
| n30_LOQ_40 | 0.40000 | 0.77620 | 2000 | 0.39688 | 0.40000 | 0.09254 | 0.10000 | 0.23333 | 0.33333 | 0.46667 | 0.56667 | 0.70000 | 19 |
| n100_LOQ_40 | 0.40000 | 0.77620 | 2000 | 0.40146 | 0.40000 | 0.04958 | 0.24000 | 0.31000 | 0.37000 | 0.43000 | 0.50000 | 0.58000 | 35 |
| n500_LOQ_40 | 0.40000 | 0.77620 | 2000 | 0.40020 | 0.40000 | 0.02149 | 0.32000 | 0.35600 | 0.38600 | 0.41400 | 0.44200 | 0.48400 | 64 |

Correlation table columns are scenario, method, statistic, number of usable replicate pairs, Pearson r, and Spearman rho. A missing correlation means fewer than three usable pairs or no variation in one of the variables.

### Censoring vs error: complete-sample reference

#### n=30, theoretical censoring=10%

| Scenario | Method | Statistic | N pairs | Pearson r | Spearman rho |
|---|---|---|---|---|---|
| n30_LOQ_10 | complete_case | cv | 2000 | -0.977 | -0.978 |
| n30_LOQ_10 | complete_case | mean | 2000 | 0.863 | 0.904 |
| n30_LOQ_10 | complete_case | median | 2000 | 0.607 | 0.679 |
| n30_LOQ_10 | complete_case | sd | 2000 | 0.328 | 0.287 |
| n30_LOQ_10 | mle_lognormal | cv | 2000 | 0.483 | 0.521 |
| n30_LOQ_10 | mle_lognormal | mean | 2000 | 0.349 | 0.428 |
| n30_LOQ_10 | mle_lognormal | median | 2000 | -0.253 | -0.269 |
| n30_LOQ_10 | mle_lognormal | sd | 2000 | 0.385 | 0.447 |
| n30_LOQ_10 | reverse_km | mean | 0 | — | — |
| n30_LOQ_10 | reverse_km | q0.025 | 0 | — | — |
| n30_LOQ_10 | reverse_km | q0.1 | 0 | — | — |
| n30_LOQ_10 | reverse_km | q0.5 | 2000 | 0.046 | 0.036 |
| n30_LOQ_10 | reverse_km | q0.9 | 2000 | 0.007 | 0.000 |
| n30_LOQ_10 | reverse_km | q0.975 | 0 | — | — |
| n30_LOQ_10 | ros_log | cv | 2000 | -0.099 | -0.046 |
| n30_LOQ_10 | ros_log | mean | 2000 | 0.110 | 0.052 |
| n30_LOQ_10 | ros_log | median | 2000 | — | — |
| n30_LOQ_10 | ros_log | sd | 2000 | -0.141 | -0.052 |
| n30_LOQ_10 | substitute_loq | cv | 2000 | -0.798 | -0.812 |
| n30_LOQ_10 | substitute_loq | mean | 2000 | 0.817 | 0.821 |
| n30_LOQ_10 | substitute_loq | median | 2000 | — | — |
| n30_LOQ_10 | substitute_loq | sd | 2000 | -0.734 | -0.754 |
| n30_LOQ_10 | substitute_loq_half | cv | 2000 | 0.608 | 0.584 |
| n30_LOQ_10 | substitute_loq_half | mean | 2000 | -0.590 | -0.577 |
| n30_LOQ_10 | substitute_loq_half | median | 2000 | — | — |
| n30_LOQ_10 | substitute_loq_half | sd | 2000 | 0.519 | 0.519 |
| n30_LOQ_10 | substitute_zero | cv | 2000 | 0.891 | 0.909 |
| n30_LOQ_10 | substitute_zero | mean | 2000 | -0.945 | -0.943 |
| n30_LOQ_10 | substitute_zero | median | 2000 | — | — |
| n30_LOQ_10 | substitute_zero | sd | 2000 | 0.848 | 0.855 |

#### n=100, theoretical censoring=10%

| Scenario | Method | Statistic | N pairs | Pearson r | Spearman rho |
|---|---|---|---|---|---|
| n100_LOQ_10 | complete_case | cv | 2000 | -0.961 | -0.972 |
| n100_LOQ_10 | complete_case | mean | 2000 | 0.890 | 0.902 |
| n100_LOQ_10 | complete_case | median | 2000 | 0.572 | 0.592 |
| n100_LOQ_10 | complete_case | sd | 2000 | 0.334 | 0.352 |
| n100_LOQ_10 | mle_lognormal | cv | 2000 | 0.285 | 0.362 |
| n100_LOQ_10 | mle_lognormal | mean | 2000 | 0.307 | 0.349 |
| n100_LOQ_10 | mle_lognormal | median | 2000 | -0.287 | -0.272 |
| n100_LOQ_10 | mle_lognormal | sd | 2000 | 0.262 | 0.342 |
| n100_LOQ_10 | reverse_km | mean | 0 | — | — |
| n100_LOQ_10 | reverse_km | q0.025 | 0 | — | — |
| n100_LOQ_10 | reverse_km | q0.1 | 0 | — | — |
| n100_LOQ_10 | reverse_km | q0.5 | 2000 | 0.036 | 0.038 |
| n100_LOQ_10 | reverse_km | q0.9 | 2000 | 0.008 | -0.004 |
| n100_LOQ_10 | reverse_km | q0.975 | 2000 | -0.031 | -0.040 |
| n100_LOQ_10 | ros_log | cv | 2000 | -0.182 | -0.131 |
| n100_LOQ_10 | ros_log | mean | 2000 | 0.187 | 0.138 |
| n100_LOQ_10 | ros_log | median | 2000 | — | — |
| n100_LOQ_10 | ros_log | sd | 2000 | -0.201 | -0.139 |
| n100_LOQ_10 | substitute_loq | cv | 2000 | -0.806 | -0.801 |
| n100_LOQ_10 | substitute_loq | mean | 2000 | 0.823 | 0.813 |
| n100_LOQ_10 | substitute_loq | median | 2000 | — | — |
| n100_LOQ_10 | substitute_loq | sd | 2000 | -0.700 | -0.690 |
| n100_LOQ_10 | substitute_loq_half | cv | 2000 | 0.608 | 0.584 |
| n100_LOQ_10 | substitute_loq_half | mean | 2000 | -0.592 | -0.570 |
| n100_LOQ_10 | substitute_loq_half | median | 2000 | — | — |
| n100_LOQ_10 | substitute_loq_half | sd | 2000 | 0.504 | 0.483 |
| n100_LOQ_10 | substitute_zero | cv | 2000 | 0.905 | 0.903 |
| n100_LOQ_10 | substitute_zero | mean | 2000 | -0.946 | -0.943 |
| n100_LOQ_10 | substitute_zero | median | 2000 | — | — |
| n100_LOQ_10 | substitute_zero | sd | 2000 | 0.788 | 0.777 |

#### n=500, theoretical censoring=10%

| Scenario | Method | Statistic | N pairs | Pearson r | Spearman rho |
|---|---|---|---|---|---|
| n500_LOQ_10 | complete_case | cv | 2000 | -0.932 | -0.950 |
| n500_LOQ_10 | complete_case | mean | 2000 | 0.909 | 0.907 |
| n500_LOQ_10 | complete_case | median | 2000 | 0.579 | 0.570 |
| n500_LOQ_10 | complete_case | sd | 2000 | 0.286 | 0.337 |
| n500_LOQ_10 | mle_lognormal | cv | 2000 | 0.212 | 0.257 |
| n500_LOQ_10 | mle_lognormal | mean | 2000 | 0.321 | 0.320 |
| n500_LOQ_10 | mle_lognormal | median | 2000 | -0.279 | -0.261 |
| n500_LOQ_10 | mle_lognormal | sd | 2000 | 0.220 | 0.265 |
| n500_LOQ_10 | reverse_km | mean | 0 | — | — |
| n500_LOQ_10 | reverse_km | q0.025 | 0 | — | — |
| n500_LOQ_10 | reverse_km | q0.1 | 0 | — | — |
| n500_LOQ_10 | reverse_km | q0.5 | 2000 | -0.047 | -0.041 |
| n500_LOQ_10 | reverse_km | q0.9 | 2000 | 0.012 | 0.006 |
| n500_LOQ_10 | reverse_km | q0.975 | 2000 | 0.015 | 0.013 |
| n500_LOQ_10 | ros_log | cv | 2000 | -0.172 | -0.154 |
| n500_LOQ_10 | ros_log | mean | 2000 | 0.175 | 0.157 |
| n500_LOQ_10 | ros_log | median | 2000 | — | — |
| n500_LOQ_10 | ros_log | sd | 2000 | -0.169 | -0.150 |
| n500_LOQ_10 | substitute_loq | cv | 2000 | -0.794 | -0.792 |
| n500_LOQ_10 | substitute_loq | mean | 2000 | 0.821 | 0.811 |
| n500_LOQ_10 | substitute_loq | median | 2000 | — | — |
| n500_LOQ_10 | substitute_loq | sd | 2000 | -0.643 | -0.629 |
| n500_LOQ_10 | substitute_loq_half | cv | 2000 | 0.606 | 0.596 |
| n500_LOQ_10 | substitute_loq_half | mean | 2000 | -0.588 | -0.579 |
| n500_LOQ_10 | substitute_loq_half | median | 2000 | — | — |
| n500_LOQ_10 | substitute_loq_half | sd | 2000 | 0.479 | 0.466 |
| n500_LOQ_10 | substitute_zero | cv | 2000 | 0.901 | 0.903 |
| n500_LOQ_10 | substitute_zero | mean | 2000 | -0.945 | -0.944 |
| n500_LOQ_10 | substitute_zero | median | 2000 | — | — |
| n500_LOQ_10 | substitute_zero | sd | 2000 | 0.697 | 0.694 |

#### n=30, theoretical censoring=40%

| Scenario | Method | Statistic | N pairs | Pearson r | Spearman rho |
|---|---|---|---|---|---|
| n30_LOQ_40 | complete_case | cv | 2000 | -0.830 | -0.855 |
| n30_LOQ_40 | complete_case | mean | 2000 | 0.606 | 0.642 |
| n30_LOQ_40 | complete_case | median | 2000 | 0.529 | 0.570 |
| n30_LOQ_40 | complete_case | sd | 2000 | 0.261 | 0.258 |
| n30_LOQ_40 | mle_lognormal | cv | 2000 | 0.433 | 0.434 |
| n30_LOQ_40 | mle_lognormal | mean | 2000 | 0.083 | 0.046 |
| n30_LOQ_40 | mle_lognormal | median | 2000 | -0.081 | -0.091 |
| n30_LOQ_40 | mle_lognormal | sd | 2000 | 0.297 | 0.311 |
| n30_LOQ_40 | reverse_km | mean | 0 | — | — |
| n30_LOQ_40 | reverse_km | q0.025 | 0 | — | — |
| n30_LOQ_40 | reverse_km | q0.1 | 0 | — | — |
| n30_LOQ_40 | reverse_km | q0.5 | 1654 | 0.048 | 0.039 |
| n30_LOQ_40 | reverse_km | q0.9 | 2000 | 0.043 | 0.030 |
| n30_LOQ_40 | reverse_km | q0.975 | 0 | — | — |
| n30_LOQ_40 | ros_log | cv | 2000 | 0.334 | 0.301 |
| n30_LOQ_40 | ros_log | mean | 2000 | -0.288 | -0.279 |
| n30_LOQ_40 | ros_log | median | 2000 | -0.194 | -0.100 |
| n30_LOQ_40 | ros_log | sd | 2000 | 0.259 | 0.294 |
| n30_LOQ_40 | substitute_loq | cv | 2000 | -0.821 | -0.822 |
| n30_LOQ_40 | substitute_loq | mean | 2000 | 0.826 | 0.818 |
| n30_LOQ_40 | substitute_loq | median | 2000 | 0.575 | 0.661 |
| n30_LOQ_40 | substitute_loq | sd | 2000 | -0.493 | -0.481 |
| n30_LOQ_40 | substitute_loq_half | cv | 2000 | 0.182 | 0.170 |
| n30_LOQ_40 | substitute_loq_half | mean | 2000 | -0.159 | -0.148 |
| n30_LOQ_40 | substitute_loq_half | median | 2000 | -0.664 | -0.659 |
| n30_LOQ_40 | substitute_loq_half | sd | 2000 | -0.011 | -0.015 |
| n30_LOQ_40 | substitute_zero | cv | 2000 | 0.865 | 0.878 |
| n30_LOQ_40 | substitute_zero | mean | 2000 | -0.873 | -0.861 |
| n30_LOQ_40 | substitute_zero | median | 2000 | -0.687 | -0.660 |
| n30_LOQ_40 | substitute_zero | sd | 2000 | 0.448 | 0.436 |

#### n=100, theoretical censoring=40%

| Scenario | Method | Statistic | N pairs | Pearson r | Spearman rho |
|---|---|---|---|---|---|
| n100_LOQ_40 | complete_case | cv | 2000 | -0.746 | -0.814 |
| n100_LOQ_40 | complete_case | mean | 2000 | 0.603 | 0.598 |
| n100_LOQ_40 | complete_case | median | 2000 | 0.524 | 0.519 |
| n100_LOQ_40 | complete_case | sd | 2000 | 0.199 | 0.204 |
| n100_LOQ_40 | mle_lognormal | cv | 2000 | 0.266 | 0.352 |
| n100_LOQ_40 | mle_lognormal | mean | 2000 | 0.041 | 0.061 |
| n100_LOQ_40 | mle_lognormal | median | 2000 | -0.096 | -0.113 |
| n100_LOQ_40 | mle_lognormal | sd | 2000 | 0.181 | 0.278 |
| n100_LOQ_40 | reverse_km | mean | 0 | — | — |
| n100_LOQ_40 | reverse_km | q0.025 | 0 | — | — |
| n100_LOQ_40 | reverse_km | q0.1 | 0 | — | — |
| n100_LOQ_40 | reverse_km | q0.5 | 1931 | 0.021 | 0.035 |
| n100_LOQ_40 | reverse_km | q0.9 | 2000 | 0.004 | 0.019 |
| n100_LOQ_40 | reverse_km | q0.975 | 2000 | -0.029 | -0.028 |
| n100_LOQ_40 | ros_log | cv | 2000 | 0.293 | 0.273 |
| n100_LOQ_40 | ros_log | mean | 2000 | -0.256 | -0.251 |
| n100_LOQ_40 | ros_log | median | 2000 | -0.027 | 0.038 |
| n100_LOQ_40 | ros_log | sd | 2000 | 0.251 | 0.266 |
| n100_LOQ_40 | substitute_loq | cv | 2000 | -0.815 | -0.813 |
| n100_LOQ_40 | substitute_loq | mean | 2000 | 0.823 | 0.810 |
| n100_LOQ_40 | substitute_loq | median | 2000 | 0.329 | 0.317 |
| n100_LOQ_40 | substitute_loq | sd | 2000 | -0.400 | -0.397 |
| n100_LOQ_40 | substitute_loq_half | cv | 2000 | 0.242 | 0.230 |
| n100_LOQ_40 | substitute_loq_half | mean | 2000 | -0.182 | -0.179 |
| n100_LOQ_40 | substitute_loq_half | median | 2000 | -0.420 | -0.317 |
| n100_LOQ_40 | substitute_loq_half | sd | 2000 | 0.017 | 0.025 |
| n100_LOQ_40 | substitute_zero | cv | 2000 | 0.874 | 0.876 |
| n100_LOQ_40 | substitute_zero | mean | 2000 | -0.876 | -0.867 |
| n100_LOQ_40 | substitute_zero | median | 2000 | -0.422 | -0.317 |
| n100_LOQ_40 | substitute_zero | sd | 2000 | 0.367 | 0.376 |

#### n=500, theoretical censoring=40%

| Scenario | Method | Statistic | N pairs | Pearson r | Spearman rho |
|---|---|---|---|---|---|
| n500_LOQ_40 | complete_case | cv | 2000 | -0.744 | -0.772 |
| n500_LOQ_40 | complete_case | mean | 2000 | 0.619 | 0.608 |
| n500_LOQ_40 | complete_case | median | 2000 | 0.519 | 0.501 |
| n500_LOQ_40 | complete_case | sd | 2000 | 0.208 | 0.228 |
| n500_LOQ_40 | mle_lognormal | cv | 2000 | 0.141 | 0.183 |
| n500_LOQ_40 | mle_lognormal | mean | 2000 | 0.008 | -0.006 |
| n500_LOQ_40 | mle_lognormal | median | 2000 | -0.064 | -0.055 |
| n500_LOQ_40 | mle_lognormal | sd | 2000 | 0.120 | 0.147 |
| n500_LOQ_40 | reverse_km | mean | 0 | — | — |
| n500_LOQ_40 | reverse_km | q0.025 | 0 | — | — |
| n500_LOQ_40 | reverse_km | q0.1 | 0 | — | — |
| n500_LOQ_40 | reverse_km | q0.5 | 2000 | -0.014 | -0.016 |
| n500_LOQ_40 | reverse_km | q0.9 | 2000 | -0.021 | -0.005 |
| n500_LOQ_40 | reverse_km | q0.975 | 2000 | 0.018 | 0.017 |
| n500_LOQ_40 | ros_log | cv | 2000 | 0.228 | 0.222 |
| n500_LOQ_40 | ros_log | mean | 2000 | -0.217 | -0.213 |
| n500_LOQ_40 | ros_log | median | 2000 | — | — |
| n500_LOQ_40 | ros_log | sd | 2000 | 0.236 | 0.235 |
| n500_LOQ_40 | substitute_loq | cv | 2000 | -0.772 | -0.783 |
| n500_LOQ_40 | substitute_loq | mean | 2000 | 0.826 | 0.811 |
| n500_LOQ_40 | substitute_loq | median | 2000 | — | — |
| n500_LOQ_40 | substitute_loq | sd | 2000 | -0.246 | -0.240 |
| n500_LOQ_40 | substitute_loq_half | cv | 2000 | 0.221 | 0.205 |
| n500_LOQ_40 | substitute_loq_half | mean | 2000 | -0.154 | -0.139 |
| n500_LOQ_40 | substitute_loq_half | median | 2000 | — | — |
| n500_LOQ_40 | substitute_loq_half | sd | 2000 | -0.003 | -0.015 |
| n500_LOQ_40 | substitute_zero | cv | 2000 | 0.844 | 0.837 |
| n500_LOQ_40 | substitute_zero | mean | 2000 | -0.872 | -0.856 |
| n500_LOQ_40 | substitute_zero | median | 2000 | — | — |
| n500_LOQ_40 | substitute_zero | sd | 2000 | 0.209 | 0.201 |

### Censoring vs error: population reference

#### n=30, theoretical censoring=10%

| Scenario | Method | Statistic | N pairs | Pearson r | Spearman rho |
|---|---|---|---|---|---|
| n30_LOQ_10 | complete_case | cv | 2000 | -0.004 | -0.009 |
| n30_LOQ_10 | complete_case | mean | 2000 | 0.011 | 0.011 |
| n30_LOQ_10 | complete_case | median | 2000 | 0.011 | 0.018 |
| n30_LOQ_10 | complete_case | sd | 2000 | 0.007 | -0.002 |
| n30_LOQ_10 | mle_lognormal | cv | 2000 | 0.597 | 0.595 |
| n30_LOQ_10 | mle_lognormal | mean | 2000 | -0.101 | -0.107 |
| n30_LOQ_10 | mle_lognormal | median | 2000 | -0.587 | -0.581 |
| n30_LOQ_10 | mle_lognormal | sd | 2000 | 0.261 | 0.265 |
| n30_LOQ_10 | reverse_km | mean | 0 | — | — |
| n30_LOQ_10 | reverse_km | q0.025 | 0 | — | — |
| n30_LOQ_10 | reverse_km | q0.1 | 0 | — | — |
| n30_LOQ_10 | reverse_km | q0.5 | 2000 | -0.320 | -0.308 |
| n30_LOQ_10 | reverse_km | q0.9 | 2000 | -0.097 | -0.094 |
| n30_LOQ_10 | reverse_km | q0.975 | 0 | — | — |
| n30_LOQ_10 | ros_log | cv | 2000 | 0.148 | 0.168 |
| n30_LOQ_10 | ros_log | mean | 2000 | -0.204 | -0.208 |
| n30_LOQ_10 | ros_log | median | 2000 | -0.321 | -0.308 |
| n30_LOQ_10 | ros_log | sd | 2000 | -0.017 | -0.021 |
| n30_LOQ_10 | substitute_loq | cv | 2000 | 0.130 | 0.148 |
| n30_LOQ_10 | substitute_loq | mean | 2000 | -0.193 | -0.197 |
| n30_LOQ_10 | substitute_loq | median | 2000 | -0.321 | -0.308 |
| n30_LOQ_10 | substitute_loq | sd | 2000 | -0.020 | -0.026 |
| n30_LOQ_10 | substitute_loq_half | cv | 2000 | 0.164 | 0.189 |
| n30_LOQ_10 | substitute_loq_half | mean | 2000 | -0.211 | -0.217 |
| n30_LOQ_10 | substitute_loq_half | median | 2000 | -0.321 | -0.308 |
| n30_LOQ_10 | substitute_loq_half | sd | 2000 | -0.014 | -0.016 |
| n30_LOQ_10 | substitute_zero | cv | 2000 | 0.201 | 0.232 |
| n30_LOQ_10 | substitute_zero | mean | 2000 | -0.229 | -0.235 |
| n30_LOQ_10 | substitute_zero | median | 2000 | -0.321 | -0.308 |
| n30_LOQ_10 | substitute_zero | sd | 2000 | -0.008 | -0.007 |

#### n=100, theoretical censoring=10%

| Scenario | Method | Statistic | N pairs | Pearson r | Spearman rho |
|---|---|---|---|---|---|
| n100_LOQ_10 | complete_case | cv | 2000 | 0.002 | 0.009 |
| n100_LOQ_10 | complete_case | mean | 2000 | 0.038 | 0.035 |
| n100_LOQ_10 | complete_case | median | 2000 | 0.031 | 0.016 |
| n100_LOQ_10 | complete_case | sd | 2000 | 0.017 | 0.018 |
| n100_LOQ_10 | mle_lognormal | cv | 2000 | 0.609 | 0.602 |
| n100_LOQ_10 | mle_lognormal | mean | 2000 | -0.087 | -0.091 |
| n100_LOQ_10 | mle_lognormal | median | 2000 | -0.584 | -0.576 |
| n100_LOQ_10 | mle_lognormal | sd | 2000 | 0.284 | 0.282 |
| n100_LOQ_10 | reverse_km | mean | 0 | — | — |
| n100_LOQ_10 | reverse_km | q0.025 | 0 | — | — |
| n100_LOQ_10 | reverse_km | q0.1 | 0 | — | — |
| n100_LOQ_10 | reverse_km | q0.5 | 2000 | -0.302 | -0.305 |
| n100_LOQ_10 | reverse_km | q0.9 | 2000 | -0.085 | -0.081 |
| n100_LOQ_10 | reverse_km | q0.975 | 2000 | -0.036 | -0.032 |
| n100_LOQ_10 | ros_log | cv | 2000 | 0.095 | 0.125 |
| n100_LOQ_10 | ros_log | mean | 2000 | -0.183 | -0.184 |
| n100_LOQ_10 | ros_log | median | 2000 | -0.303 | -0.306 |
| n100_LOQ_10 | ros_log | sd | 2000 | -0.006 | -0.006 |
| n100_LOQ_10 | substitute_loq | cv | 2000 | 0.085 | 0.113 |
| n100_LOQ_10 | substitute_loq | mean | 2000 | -0.173 | -0.174 |
| n100_LOQ_10 | substitute_loq | median | 2000 | -0.303 | -0.306 |
| n100_LOQ_10 | substitute_loq | sd | 2000 | -0.008 | -0.009 |
| n100_LOQ_10 | substitute_loq_half | cv | 2000 | 0.105 | 0.139 |
| n100_LOQ_10 | substitute_loq_half | mean | 2000 | -0.191 | -0.193 |
| n100_LOQ_10 | substitute_loq_half | median | 2000 | -0.303 | -0.306 |
| n100_LOQ_10 | substitute_loq_half | sd | 2000 | -0.004 | -0.003 |
| n100_LOQ_10 | substitute_zero | cv | 2000 | 0.126 | 0.165 |
| n100_LOQ_10 | substitute_zero | mean | 2000 | -0.210 | -0.211 |
| n100_LOQ_10 | substitute_zero | median | 2000 | -0.303 | -0.306 |
| n100_LOQ_10 | substitute_zero | sd | 2000 | 0.000 | 0.003 |

#### n=500, theoretical censoring=10%

| Scenario | Method | Statistic | N pairs | Pearson r | Spearman rho |
|---|---|---|---|---|---|
| n500_LOQ_10 | complete_case | cv | 2000 | -0.022 | 0.006 |
| n500_LOQ_10 | complete_case | mean | 2000 | 0.039 | 0.048 |
| n500_LOQ_10 | complete_case | median | 2000 | 0.009 | 0.000 |
| n500_LOQ_10 | complete_case | sd | 2000 | -0.008 | 0.016 |
| n500_LOQ_10 | mle_lognormal | cv | 2000 | 0.649 | 0.634 |
| n500_LOQ_10 | mle_lognormal | mean | 2000 | -0.092 | -0.085 |
| n500_LOQ_10 | mle_lognormal | median | 2000 | -0.618 | -0.617 |
| n500_LOQ_10 | mle_lognormal | sd | 2000 | 0.320 | 0.318 |
| n500_LOQ_10 | reverse_km | mean | 0 | — | — |
| n500_LOQ_10 | reverse_km | q0.025 | 0 | — | — |
| n500_LOQ_10 | reverse_km | q0.1 | 0 | — | — |
| n500_LOQ_10 | reverse_km | q0.5 | 2000 | -0.332 | -0.332 |
| n500_LOQ_10 | reverse_km | q0.9 | 2000 | -0.057 | -0.042 |
| n500_LOQ_10 | reverse_km | q0.975 | 2000 | -0.023 | -0.009 |
| n500_LOQ_10 | ros_log | cv | 2000 | 0.040 | 0.088 |
| n500_LOQ_10 | ros_log | mean | 2000 | -0.196 | -0.183 |
| n500_LOQ_10 | ros_log | median | 2000 | -0.329 | -0.329 |
| n500_LOQ_10 | ros_log | sd | 2000 | -0.028 | -0.008 |
| n500_LOQ_10 | substitute_loq | cv | 2000 | 0.033 | 0.079 |
| n500_LOQ_10 | substitute_loq | mean | 2000 | -0.185 | -0.172 |
| n500_LOQ_10 | substitute_loq | median | 2000 | -0.329 | -0.329 |
| n500_LOQ_10 | substitute_loq | sd | 2000 | -0.030 | -0.010 |
| n500_LOQ_10 | substitute_loq_half | cv | 2000 | 0.046 | 0.096 |
| n500_LOQ_10 | substitute_loq_half | mean | 2000 | -0.204 | -0.192 |
| n500_LOQ_10 | substitute_loq_half | median | 2000 | -0.329 | -0.329 |
| n500_LOQ_10 | substitute_loq_half | sd | 2000 | -0.027 | -0.006 |
| n500_LOQ_10 | substitute_zero | cv | 2000 | 0.059 | 0.114 |
| n500_LOQ_10 | substitute_zero | mean | 2000 | -0.224 | -0.211 |
| n500_LOQ_10 | substitute_zero | median | 2000 | -0.329 | -0.329 |
| n500_LOQ_10 | substitute_zero | sd | 2000 | -0.024 | -0.002 |

#### n=30, theoretical censoring=40%

| Scenario | Method | Statistic | N pairs | Pearson r | Spearman rho |
|---|---|---|---|---|---|
| n30_LOQ_40 | complete_case | cv | 2000 | -0.027 | -0.036 |
| n30_LOQ_40 | complete_case | mean | 2000 | 0.026 | -0.005 |
| n30_LOQ_40 | complete_case | median | 2000 | 0.044 | 0.020 |
| n30_LOQ_40 | complete_case | sd | 2000 | 0.008 | -0.020 |
| n30_LOQ_40 | mle_lognormal | cv | 2000 | 0.484 | 0.471 |
| n30_LOQ_40 | mle_lognormal | mean | 2000 | -0.430 | -0.462 |
| n30_LOQ_40 | mle_lognormal | median | 2000 | -0.938 | -0.947 |
| n30_LOQ_40 | mle_lognormal | sd | 2000 | 0.144 | 0.077 |
| n30_LOQ_40 | reverse_km | mean | 0 | — | — |
| n30_LOQ_40 | reverse_km | q0.025 | 0 | — | — |
| n30_LOQ_40 | reverse_km | q0.1 | 0 | — | — |
| n30_LOQ_40 | reverse_km | q0.5 | 1654 | -0.694 | -0.724 |
| n30_LOQ_40 | reverse_km | q0.9 | 2000 | -0.291 | -0.310 |
| n30_LOQ_40 | reverse_km | q0.975 | 0 | — | — |
| n30_LOQ_40 | ros_log | cv | 2000 | 0.256 | 0.268 |
| n30_LOQ_40 | ros_log | mean | 2000 | -0.497 | -0.523 |
| n30_LOQ_40 | ros_log | median | 2000 | -0.807 | -0.816 |
| n30_LOQ_40 | ros_log | sd | 2000 | -0.070 | -0.093 |
| n30_LOQ_40 | substitute_loq | cv | 2000 | 0.078 | 0.080 |
| n30_LOQ_40 | substitute_loq | mean | 2000 | -0.387 | -0.410 |
| n30_LOQ_40 | substitute_loq | median | 2000 | -0.766 | -0.827 |
| n30_LOQ_40 | substitute_loq | sd | 2000 | -0.091 | -0.125 |
| n30_LOQ_40 | substitute_loq_half | cv | 2000 | 0.238 | 0.269 |
| n30_LOQ_40 | substitute_loq_half | mean | 2000 | -0.463 | -0.489 |
| n30_LOQ_40 | substitute_loq_half | median | 2000 | -0.845 | -0.833 |
| n30_LOQ_40 | substitute_loq_half | sd | 2000 | -0.079 | -0.108 |
| n30_LOQ_40 | substitute_zero | cv | 2000 | 0.476 | 0.538 |
| n30_LOQ_40 | substitute_zero | mean | 2000 | -0.531 | -0.556 |
| n30_LOQ_40 | substitute_zero | median | 2000 | -0.844 | -0.833 |
| n30_LOQ_40 | substitute_zero | sd | 2000 | -0.066 | -0.090 |

#### n=100, theoretical censoring=40%

| Scenario | Method | Statistic | N pairs | Pearson r | Spearman rho |
|---|---|---|---|---|---|
| n100_LOQ_40 | complete_case | cv | 2000 | -0.025 | -0.066 |
| n100_LOQ_40 | complete_case | mean | 2000 | -0.000 | -0.002 |
| n100_LOQ_40 | complete_case | median | 2000 | 0.024 | 0.015 |
| n100_LOQ_40 | complete_case | sd | 2000 | -0.008 | -0.047 |
| n100_LOQ_40 | mle_lognormal | cv | 2000 | 0.487 | 0.458 |
| n100_LOQ_40 | mle_lognormal | mean | 2000 | -0.470 | -0.450 |
| n100_LOQ_40 | mle_lognormal | median | 2000 | -0.950 | -0.949 |
| n100_LOQ_40 | mle_lognormal | sd | 2000 | 0.086 | 0.072 |
| n100_LOQ_40 | reverse_km | mean | 0 | — | — |
| n100_LOQ_40 | reverse_km | q0.025 | 0 | — | — |
| n100_LOQ_40 | reverse_km | q0.1 | 0 | — | — |
| n100_LOQ_40 | reverse_km | q0.5 | 1931 | -0.792 | -0.801 |
| n100_LOQ_40 | reverse_km | q0.9 | 2000 | -0.273 | -0.257 |
| n100_LOQ_40 | reverse_km | q0.975 | 2000 | -0.143 | -0.144 |
| n100_LOQ_40 | ros_log | cv | 2000 | 0.148 | 0.154 |
| n100_LOQ_40 | ros_log | mean | 2000 | -0.501 | -0.492 |
| n100_LOQ_40 | ros_log | median | 2000 | -0.818 | -0.814 |
| n100_LOQ_40 | ros_log | sd | 2000 | -0.074 | -0.118 |
| n100_LOQ_40 | substitute_loq | cv | 2000 | 0.042 | 0.020 |
| n100_LOQ_40 | substitute_loq | mean | 2000 | -0.399 | -0.390 |
| n100_LOQ_40 | substitute_loq | median | 2000 | -0.813 | -0.816 |
| n100_LOQ_40 | substitute_loq | sd | 2000 | -0.088 | -0.137 |
| n100_LOQ_40 | substitute_loq_half | cv | 2000 | 0.139 | 0.149 |
| n100_LOQ_40 | substitute_loq_half | mean | 2000 | -0.472 | -0.464 |
| n100_LOQ_40 | substitute_loq_half | median | 2000 | -0.818 | -0.816 |
| n100_LOQ_40 | substitute_loq_half | sd | 2000 | -0.080 | -0.126 |
| n100_LOQ_40 | substitute_zero | cv | 2000 | 0.279 | 0.343 |
| n100_LOQ_40 | substitute_zero | mean | 2000 | -0.536 | -0.528 |
| n100_LOQ_40 | substitute_zero | median | 2000 | -0.764 | -0.816 |
| n100_LOQ_40 | substitute_zero | sd | 2000 | -0.071 | -0.115 |

#### n=500, theoretical censoring=40%

| Scenario | Method | Statistic | N pairs | Pearson r | Spearman rho |
|---|---|---|---|---|---|
| n500_LOQ_40 | complete_case | cv | 2000 | 0.025 | 0.034 |
| n500_LOQ_40 | complete_case | mean | 2000 | 0.030 | 0.031 |
| n500_LOQ_40 | complete_case | median | 2000 | 0.003 | 0.006 |
| n500_LOQ_40 | complete_case | sd | 2000 | 0.030 | 0.038 |
| n500_LOQ_40 | mle_lognormal | cv | 2000 | 0.497 | 0.484 |
| n500_LOQ_40 | mle_lognormal | mean | 2000 | -0.455 | -0.439 |
| n500_LOQ_40 | mle_lognormal | median | 2000 | -0.951 | -0.944 |
| n500_LOQ_40 | mle_lognormal | sd | 2000 | 0.098 | 0.094 |
| n500_LOQ_40 | reverse_km | mean | 0 | — | — |
| n500_LOQ_40 | reverse_km | q0.025 | 0 | — | — |
| n500_LOQ_40 | reverse_km | q0.1 | 0 | — | — |
| n500_LOQ_40 | reverse_km | q0.5 | 2000 | -0.818 | -0.806 |
| n500_LOQ_40 | reverse_km | q0.9 | 2000 | -0.261 | -0.251 |
| n500_LOQ_40 | reverse_km | q0.975 | 2000 | -0.088 | -0.084 |
| n500_LOQ_40 | ros_log | cv | 2000 | 0.140 | 0.168 |
| n500_LOQ_40 | ros_log | mean | 2000 | -0.471 | -0.456 |
| n500_LOQ_40 | ros_log | median | 2000 | -0.817 | -0.805 |
| n500_LOQ_40 | ros_log | sd | 2000 | -0.031 | -0.028 |
| n500_LOQ_40 | substitute_loq | cv | 2000 | 0.070 | 0.087 |
| n500_LOQ_40 | substitute_loq | mean | 2000 | -0.368 | -0.356 |
| n500_LOQ_40 | substitute_loq | median | 2000 | -0.817 | -0.805 |
| n500_LOQ_40 | substitute_loq | sd | 2000 | -0.041 | -0.039 |
| n500_LOQ_40 | substitute_loq_half | cv | 2000 | 0.137 | 0.166 |
| n500_LOQ_40 | substitute_loq_half | mean | 2000 | -0.443 | -0.428 |
| n500_LOQ_40 | substitute_loq_half | median | 2000 | -0.817 | -0.805 |
| n500_LOQ_40 | substitute_loq_half | sd | 2000 | -0.036 | -0.033 |
| n500_LOQ_40 | substitute_zero | cv | 2000 | 0.229 | 0.279 |
| n500_LOQ_40 | substitute_zero | mean | 2000 | -0.509 | -0.493 |
| n500_LOQ_40 | substitute_zero | median | 2000 | -0.817 | -0.805 |
| n500_LOQ_40 | substitute_zero | sd | 2000 | -0.030 | -0.026 |

## Reverse Kaplan–Meier details

All accuracy values below are conditional on a finite estimate with `estimand_status = success`; usable counts and availability identify how many of the 2,000 expected replicates contribute. Availability differs by quantile and scenario, so these conditional accuracy metrics are not unconditional performance summaries and should be read together with the availability table. The mean has no identified nonparametric target under the common LOQ. The low quantiles are theoretically non-identifiable and are excluded even when the implementation returned a finite number. Median non-estimability and missing 97.5th quantiles are retained as availability states, not imputed or automatically treated as computational failures.

| Scenario | Estimand | Expected | Finite | Usable | Finite, status-excluded | Theoretical non-ID | Structural non-ID | Not estimable | Not returned | Failures | Available (%) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| n30_LOQ_10 | mean | 2000 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0.00 |
| n30_LOQ_10 | q0.025 | 2000 | 72 | 0 | 72 | 2000 | 0 | 0 | 0 | 0 | 0.00 |
| n30_LOQ_10 | q0.1 | 2000 | 831 | 0 | 831 | 2000 | 0 | 0 | 0 | 0 | 0.00 |
| n30_LOQ_10 | q0.5 | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| n30_LOQ_10 | q0.9 | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| n30_LOQ_10 | q0.975 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0.00 |
| n100_LOQ_10 | mean | 2000 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0.00 |
| n100_LOQ_10 | q0.025 | 2000 | 4 | 0 | 4 | 2000 | 0 | 0 | 0 | 0 | 0.00 |
| n100_LOQ_10 | q0.1 | 2000 | 901 | 0 | 901 | 2000 | 0 | 0 | 0 | 0 | 0.00 |
| n100_LOQ_10 | q0.5 | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| n100_LOQ_10 | q0.9 | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| n100_LOQ_10 | q0.975 | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| n500_LOQ_10 | mean | 2000 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0.00 |
| n500_LOQ_10 | q0.025 | 2000 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0.00 |
| n500_LOQ_10 | q0.1 | 2000 | 973 | 0 | 973 | 2000 | 0 | 0 | 0 | 0 | 0.00 |
| n500_LOQ_10 | q0.5 | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| n500_LOQ_10 | q0.9 | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| n500_LOQ_10 | q0.975 | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| n30_LOQ_40 | mean | 2000 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0.00 |
| n30_LOQ_40 | q0.025 | 2000 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0.00 |
| n30_LOQ_40 | q0.1 | 2000 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0.00 |
| n30_LOQ_40 | q0.5 | 2000 | 1654 | 1654 | 0 | 0 | 0 | 346 | 0 | 0 | 82.70 |
| n30_LOQ_40 | q0.9 | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| n30_LOQ_40 | q0.975 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 2000 | 0 | 0.00 |
| n100_LOQ_40 | mean | 2000 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0.00 |
| n100_LOQ_40 | q0.025 | 2000 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0.00 |
| n100_LOQ_40 | q0.1 | 2000 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0.00 |
| n100_LOQ_40 | q0.5 | 2000 | 1931 | 1931 | 0 | 0 | 0 | 69 | 0 | 0 | 96.55 |
| n100_LOQ_40 | q0.9 | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| n100_LOQ_40 | q0.975 | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| n500_LOQ_40 | mean | 2000 | 0 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0.00 |
| n500_LOQ_40 | q0.025 | 2000 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0.00 |
| n500_LOQ_40 | q0.1 | 2000 | 0 | 0 | 0 | 2000 | 0 | 0 | 0 | 0 | 0.00 |
| n500_LOQ_40 | q0.5 | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| n500_LOQ_40 | q0.9 | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |
| n500_LOQ_40 | q0.975 | 2000 | 2000 | 2000 | 0 | 0 | 0 | 0 | 0 | 0 | 100.00 |

Across the six scenarios: mean non-identifiable = 12,000; q2.5% theoretically non-identifiable = 12,000; q10% theoretically non-identifiable = 12,000; median not estimable in sample = 415; q97.5% not returned = 4,000. In addition, finite values coexist with theoretical non-identifiability for q2.5% (76 records) and q10% (2,705 records); they were excluded from every primary accuracy metric.

## Figures

Moment-error plots include complete case, the three substitution methods, ROS, and lognormal MLE. Reverse KM is kept separate because it does not estimate those reconstructed-sample moments. All 13 boxplot figures show Tukey outlier points; no observations are removed and the y-axes span the full plotted data range, with scales set independently by panel. Thus the boxplots show complete tails, while central boxes may appear compressed where extremes are large. Censoring scatterplots already display all eligible points and are unchanged. Reverse-KM plots include only identified quantiles with success status; their accuracy is conditional on the usable n shown and must be interpreted with the availability table.

### n=30, theoretical censoring=10%

**Moment errors vs complete-sample truth**  
![Moment errors vs sample truth: n30_LOQ_10](figures/primary_monte_carlo/moment_errors_sample_truth_n30_LOQ_10.png)

**Moment errors vs theoretical population**  
![Moment errors vs population: n30_LOQ_10](figures/primary_monte_carlo/moment_errors_population_n30_LOQ_10.png)

**Empirical censoring vs population error**  
![Censoring and population error: n30_LOQ_10](figures/primary_monte_carlo/censoring_vs_population_error_n30_LOQ_10.png)

### n=100, theoretical censoring=10%

**Moment errors vs complete-sample truth**  
![Moment errors vs sample truth: n100_LOQ_10](figures/primary_monte_carlo/moment_errors_sample_truth_n100_LOQ_10.png)

**Moment errors vs theoretical population**  
![Moment errors vs population: n100_LOQ_10](figures/primary_monte_carlo/moment_errors_population_n100_LOQ_10.png)

**Empirical censoring vs population error**  
![Censoring and population error: n100_LOQ_10](figures/primary_monte_carlo/censoring_vs_population_error_n100_LOQ_10.png)

### n=500, theoretical censoring=10%

**Moment errors vs complete-sample truth**  
![Moment errors vs sample truth: n500_LOQ_10](figures/primary_monte_carlo/moment_errors_sample_truth_n500_LOQ_10.png)

**Moment errors vs theoretical population**  
![Moment errors vs population: n500_LOQ_10](figures/primary_monte_carlo/moment_errors_population_n500_LOQ_10.png)

**Empirical censoring vs population error**  
![Censoring and population error: n500_LOQ_10](figures/primary_monte_carlo/censoring_vs_population_error_n500_LOQ_10.png)

### n=30, theoretical censoring=40%

**Moment errors vs complete-sample truth**  
![Moment errors vs sample truth: n30_LOQ_40](figures/primary_monte_carlo/moment_errors_sample_truth_n30_LOQ_40.png)

**Moment errors vs theoretical population**  
![Moment errors vs population: n30_LOQ_40](figures/primary_monte_carlo/moment_errors_population_n30_LOQ_40.png)

**Empirical censoring vs population error**  
![Censoring and population error: n30_LOQ_40](figures/primary_monte_carlo/censoring_vs_population_error_n30_LOQ_40.png)

### n=100, theoretical censoring=40%

**Moment errors vs complete-sample truth**  
![Moment errors vs sample truth: n100_LOQ_40](figures/primary_monte_carlo/moment_errors_sample_truth_n100_LOQ_40.png)

**Moment errors vs theoretical population**  
![Moment errors vs population: n100_LOQ_40](figures/primary_monte_carlo/moment_errors_population_n100_LOQ_40.png)

**Empirical censoring vs population error**  
![Censoring and population error: n100_LOQ_40](figures/primary_monte_carlo/censoring_vs_population_error_n100_LOQ_40.png)

### n=500, theoretical censoring=40%

**Moment errors vs complete-sample truth**  
![Moment errors vs sample truth: n500_LOQ_40](figures/primary_monte_carlo/moment_errors_sample_truth_n500_LOQ_40.png)

**Moment errors vs theoretical population**  
![Moment errors vs population: n500_LOQ_40](figures/primary_monte_carlo/moment_errors_population_n500_LOQ_40.png)

**Empirical censoring vs population error**  
![Censoring and population error: n500_LOQ_40](figures/primary_monte_carlo/censoring_vs_population_error_n500_LOQ_40.png)

### Reverse-KM identified quantiles

Quantile-error distributions include only rows with a usable, status-success estimate. The n below each box is the contributing replicate count; accuracy is conditional on that availability. Outlier points are shown and no usable values are removed.

![Reverse-KM quantile errors vs population](figures/primary_monte_carlo/reverse_km_quantile_errors_vs_population.png)

## Quantitative reading across the reported strata

The points below describe specific strata and are not an overall method ordering:

- For the complete-sample mean at theoretical 10% censoring, the SD across generated samples was 0.4044 at n=30 and 0.0934 at n=500. This sampling spread is separate from each method's error.
- At n=30 and theoretical 10% censoring, population-reference RMSE for the six moment-output methods ranged from 0.4034 to 0.4727 for the mean; the complete-sample mean itself had SD 0.4044. Read this beside the bias and paired-error tables rather than treating RMSE as censoring-only error.
- At n=500 and theoretical 40% censoring, population-reference mean biases (with RMSE in parentheses) were: complete case 0.8096 (0.8219); zero substitution -0.1743 (0.2002); LOQ/2 -0.0190 (0.0964); LOQ 0.1363 (0.1640); ROS -0.0028 (0.0921); and lognormal MLE -0.0003 (0.0897). These are one scenario/statistic only; MLE moments come from the fitted distribution, not imputation.
- Mean empirical censoring was 10.0233% (SD 5.4081 percentage points) at n=30 and 9.9872% (SD 1.3582 percentage points) at n=500; the realized proportion varies by replicate and its error associations are descriptive.

## Limitations

The audited result files do not retain the latent observations. Therefore it is not possible to independently regenerate the complete-sample `sample_truth` statistics for each replicate. Their stored values agree across the six applicable methods, the theoretical targets are consistent with the configured lognormal distribution, and the stored errors equal estimate minus their designated reference on every record; this verifies internal arithmetic but not the original latent-sample calculation.

The audit classified the Monte Carlo outputs as **suitable with reservations**. Those reservations remain: Reverse-KM numeric availability is not equivalent to estimand identifiability, its accuracy summaries are conditional on usable estimates, and sample-truth provenance cannot be independently revalidated from the retained files. The current combined result content matches the checkpoint payloads, but the historical SHA-256 discrepancy is unresolved and does not establish byte-level provenance. Inference, confidence intervals, UCLs, and bootstrap were not performed. Results describe this single lognormal design and do not establish performance under other generating distributions.

## Conclusions

The stored results are complete and internally reconciled for all six design conditions. The two error references yield distinct questions and are reported separately. Tables show the metric values by scenario, method, and estimand; sampling distributions of the complete-sample statistics are also shown separately, so sample variability is not attributed to a censoring method. The availability tables and figures make clear where a method has fewer usable estimates.

The analysis supports interpretation of these scenario-specific summaries with the stated reservations. Reverse-KM records marked non-identifiable, not estimable, or not returned are excluded; its accuracy is conditional on the remaining usable replicates and is not an unconditional cross-method comparison. No generic cross-method ranking is warranted. No further simulation or inferential claim is made.

