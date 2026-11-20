# Experimental protocol — Part 1: observations below the LOQ

## Objective

Evaluate by Monte Carlo simulation how alternative treatments of left-censored concentration measurements affect descriptive statistics and inference. A censored observation provides the information that its latent value is below the applicable LOQ; it is not an ordinary missing value.

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
- Use the same generated sample for all applicable methods to enable paired comparisons.

The initial lognormal experiment contains 12,000 generated samples. Review Monte Carlo stability at 2,000, 5,000 and 10,000 replicates.

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

### Inference
- Two-sided 95% confidence interval for the mean.
- Interval width and empirical coverage.
- One-sided 95% upper confidence limit (UCL).
- Empirical coverage of the UCL relative to the true population mean.

Inferential procedures are method-specific. A common t interval is not imposed on every method. Bootstrap procedures must repeat the relevant estimation workflow and preserve censoring information.

## Failure and identifiability accounting

Record convergence failures, non-finite estimates, non-identifiable statistics, failed bootstrap replicates, samples with all observations censored and samples with no censored observations. These outcomes must not be silently removed.

## Later model-misspecification scenario

A subsequent scenario will generate data from a Gamma distribution while retaining a lognormal likelihood for MLE. Gamma parameters will be selected and documented explicitly. The LOQs remain fixed at their initial absolute values; censoring probabilities are consequences of the generating distribution, not forced to remain at 10% and 40%.

## Validation sequence

1. Confirm package versions and function interfaces.
2. Inspect documentation and source code for exact estimators.
3. Test a common-LOQ sample with known latent values.
4. Test no-censoring, all-censored, low-censoring and higher-censoring edge cases.
5. Compare package outputs and identify the estimator behind each statistic.
6. Validate a pilot against theoretical quantities and independently calculated likelihoods.
7. Only then run the full Monte Carlo experiment.

## Initial references

- Helsel (2006), *Fabricating data: How substituting values for nondetects can ruin results, and what can be done about it*.
- Hornung and Reed (1990), *Estimation of average concentration in the presence of nondetectable values*.
- Helsel and Cohn (1988), *Estimation of descriptive statistics for multiply censored water quality data*.
- Antweiler and Taylor (2008), *Evaluation of Statistical Treatments of Left-Censored Environmental Data using Coincident Uncensored Data Sets: I. Summary Statistics*.
- Shoari, Dubé and Chenouri (2015), *Estimating the mean and standard deviation of environmental data with below detection limit observations: Considering highly skewed data and model misspecification*.
- Gillespie et al. (2010), work on estimating population distributions with a reverse Kaplan–Meier estimator.

Package documentation: [NADA](https://cran.r-project.org/package=NADA), [NADA2](https://cran.r-project.org/package=NADA2), [EnvStats](https://cran.r-project.org/package=EnvStats).

## Decisions still requiring validation

- Exact function interfaces and returned object structures for installed package versions.
- How each package estimates moments under a common left-censoring threshold.
- Confidence interval and UCL construction for ROS, MLE and reverse KM.
- Bootstrap design and computational settings.
- Treatment of samples with no detected observations or no censored observations.
