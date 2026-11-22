# Experimental protocol — Part 1: observations below the LOQ

## Objective

The first Monte Carlo evaluates point estimators for left-censored concentration measurements. A censored observation provides the information that its latent value is below the applicable LOQ; it is not an ordinary missing value. Inferential validation remains documented in the existing reports but is excluded from this primary run.

## Initial data-generating model

The first scenario is log(X) ~ Normal(mu = 0, sigma = 1).

Theoretical quantities:
- Mean: exp(mu + sigma^2 / 2) = 1.6487213.
- Median: exp(mu) = 1.
- Standard deviation: sqrt((exp(sigma^2) - 1) * exp(2 * mu + sigma^2)) = 2.1611974.
- Coefficient of variation: sqrt(exp(sigma^2) - 1) = 1.3108325.

This is a plausible positive, right-skewed scenario, not a universal model for environmental concentrations.

## LOQ conditions

Two fixed LOQs are defined from the generating distribution:
- Low LOQ: qlnorm(0.10, meanlog = 0, sdlog = 1), approximately 0.2776.
- High LOQ: qlnorm(0.40, meanlog = 0, sdlog = 1), approximately 0.7761.

The theoretical censoring probabilities are 0.10 and 0.40. LOQs remain fixed across Monte Carlo replicates. They must not be recalibrated from individual samples to force a target observed censoring percentage. Record both theoretical and realized censoring proportions.

## Monte Carlo design

- Sample sizes: n = 30, 100, 500.
- LOQ conditions: 10th and 40th percentiles.
- Initial replicates: 2,000 per sample-size/LOQ combination.
- Approved master seed: 20261122L.
- Total: six conditions and 12,000 generated samples.
- Use the same generated sample for all applicable methods to enable paired comparisons.
- For each replicate, generate the full latent sample, calculate its statistics, then apply the fixed scenario LOQ and run the methods. Compare each estimate with both the complete-sample statistic (paired sample truth) and the population parameter, stored as distinct references.

The first experiment is fixed at 2,000 replicates per condition. Do not extend it to 5,000 or 10,000 replicates.

## Methods retained

No additional method is introduced without scientific rationale and prior discussion.

1. Complete-case analysis.
2. Substitute zero.
3. Substitute LOQ/2.
4. Substitute LOQ.
5. Regression on Order Statistics (ROS), using log-ROS in the lognormal scenario.
6. Maximum likelihood estimation under a left-censored lognormal model.
7. Reverse Kaplan–Meier (reverse KM).

These methods do not necessarily return the same type of output. Parametric likelihood estimation and nonparametric distribution estimation should not be judged solely as data-completion procedures.

## Reverse KM with a common LOQ

A common LOQ does not prevent applying reverse KM. However, when all censored observations share the same threshold, their latent ordering below that threshold is unknown. Reverse KM cannot nonparametrically recover the distributional shape below the LOQ. Quantiles below the censoring boundary are not identified from the observed data alone. Moments requiring integration over the unobserved region need explicit conventions or additional assumptions.

The implementation review must distinguish:
- the reverse-KM distribution estimate and identifiable quantiles;
- the estimator a package uses for moments;
- the procedure used for confidence intervals and one-sided upper confidence limits.

NADA, NADA2 and EnvStats will be inspected at function and source/documentation level. A function combining KM quantiles with a different moment estimator must not be described as an independent reverse-KM estimator.

## Outcomes

### Descriptive accuracy
- Mean, median, standard deviation and coefficient of variation where meaningful.
- 2.5th and 97.5th percentiles, with identifiability qualifications.
- Distributional discrepancy, including tail behaviour.

### Estimation performance
For each estimable parameter: bias, relative bias, empirical standard deviation of estimates and RMSE.

### Inference criteria for separate work (excluded from the primary run)

Inference is excluded from the primary Monte Carlo: no confidence intervals, UCLs, bootstrap refits, or coverage simulation are run. Previously specified inferential outcomes include a two-sided 95% confidence interval for the mean, interval width and empirical coverage, and a one-sided 95% upper confidence limit (UCL) with empirical coverage relative to the true population mean. Existing inference reports remain unchanged and document prior inferential validation work; they are not outcomes of this point-estimator study. Inferential procedures remain method-specific; a common t interval is not imposed on every method.

## Failure and identifiability accounting

For the primary run, record convergence failures, non-finite estimates, non-identifiable statistics, and samples with all observations censored or no censored observations. These outcomes must not be silently removed. Bootstrap failures do not apply because the primary run performs no bootstrap.

## Later model-misspecification scenario

A subsequent scenario will generate data from a Gamma distribution while retaining a lognormal likelihood for MLE. Gamma parameters will be selected and documented explicitly. The LOQs remain fixed at their initial absolute values; censoring probabilities are consequences of the generating distribution, not forced to remain at 10% and 40%.

## Validation sequence

1. Confirm package versions and function interfaces.
2. Inspect documentation and source code for exact estimators.
3. Test a common-LOQ sample with known latent values.
4. Test no-censoring, all-censored, low-censoring and higher-censoring edge cases.
5. Compare package outputs and identify the estimator behind each statistic.
6. Validate a pilot against theoretical quantities and independently calculated likelihoods.
7. After the functional and scientific review, run the approved primary point-estimator Monte Carlo.

## Initial references

- Helsel (2006), *Fabricating data: How substituting values for nondetects can ruin results, and what can be done about it*.
- Hornung and Reed (1990), *Estimation of average concentration in the presence of nondetectable values*.
- Helsel and Cohn (1988), *Estimation of descriptive statistics for multiply censored water quality data*.
- Antweiler and Taylor (2008), *Evaluation of Statistical Treatments of Left-Censored Environmental Data using Coincident Uncensored Data Sets: I. Summary Statistics*.
- Shoari, Dubé and Chenouri (2015), *Estimating the mean and standard deviation of environmental data with below detection limit observations: Considering highly skewed data and model misspecification*.
- Gillespie et al. (2010), work on estimating population distributions with a reverse Kaplan–Meier estimator.

Package documentation: [NADA](https://cran.r-project.org/package=NADA), [NADA2](https://cran.r-project.org/package=NADA2), [EnvStats](https://cran.r-project.org/package=EnvStats).

## Prior validation record

Prior estimator and inference evidence is retained in `reports/estimator_validation.md`, `reports/inference_validation.md`, and `reports/inference_bootstrap_stability.md`; this protocol change does not revise those results. Functional checks cover package interfaces and edge cases, but do not establish estimator validity. The bootstrap stability supplement found material endpoint changes in its targeted B=99 versus B=499 comparison; bootstrap settings are irrelevant to the primary run because inference is disabled.

## Current implementation and validation status

The repository implements the seven methods listed above and a primary-run driver in `R/12_primary_monte_carlo.R`. The driver sends one generated sample to all selected methods, stores latent sample values only for evaluation, derives replicate seeds from the approved master seed, and supports atomic per-scenario checkpoints and signature-checked resume. The exact run uses 2,000 replicates per condition. A small functional smoke check verifies interruption, resume, and reproducibility; it is not evidence of estimator accuracy or inferential validity.

Estimator and inference evidence is summarized in `reports/estimator_validation.md`, `reports/inference_validation.md`, and `reports/inference_bootstrap_stability.md`. The estimator diagnostics have stated reservations and do not validate intervals. Inference decisions remain specific to method, estimand, interval type, and condition; `not validated` and `inconclusive` are distinct outcomes. These reports remain historical evidence and are not modified by the primary-run preparation.

For reverse KM under the common LOQ, the mean is recorded as structurally non-identifiable and only quantiles returned by the fitted curve are considered; structural non-identifiability, sample non-estimability, and a missing returned quantile are kept separate. The runner rejects `inference = TRUE`; the primary run fixes `inference = FALSE` and has no inferential outcomes.

The definitive lognormal Monte Carlo has not been launched. Its approved seed is `20261122L`, its replicate count is fixed at 2,000 per condition, and its dedicated output directory is `results/primary_monte_carlo_lognormal_seed_20261122/`. The directory is reserved and contains no simulation results. After review, run `Rscript simulations/run_primary_monte_carlo.R` from the repository root. Do not execute it during preparation. The later Gamma misspecification scenario remains separate and is not part of this run.
