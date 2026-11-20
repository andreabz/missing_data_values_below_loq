# Values below the limit of quantification (LOQ)

Reproducible R study of statistical methods for left-censored environmental and chemical concentration data.

## Scope

The project evaluates the consequences of alternative treatments of observations below a known limit of quantification. The first simulation scenario uses a lognormal data-generating distribution. A later scenario will investigate model misspecification using Gamma-generated data analysed with a lognormal likelihood.

A censored observation provides the information that its latent value is below the applicable LOQ. It is not equivalent to an ordinary missing value.

## Initial experimental design

- Data-generating distribution: Lognormal(meanlog = 0, sdlog = 1)
- Sample sizes: 30, 100, 500
- Fixed LOQs: 10th and 40th percentiles of the generating distribution
- Initial Monte Carlo replicates: 2,000 per sample-size/LOQ scenario
- Initial bootstrap setting for pilot investigations: 2,000 resamples, subject to stability checks
- Methods: complete cases, zero, LOQ/2, LOQ, log-ROS, left-censored lognormal MLE and reverse Kaplan–Meier

LOQs are fixed before sampling and are not recalculated within Monte Carlo replicates. The observed censoring proportion is allowed to vary across samples.

## Repository layout

- `protocol/`: scientific protocol and methodological decisions
- `R/`: reusable R functions and simulation code
- `tests/manual/`: exploratory checks of package APIs and edge cases
- `simulations/`: simulation runners
- `results/`: generated summaries and figures
- `reports/`: Quarto documents

## Getting started

Open `missing_data_values_below_loq.Rproj` in RStudio. Install `data.table` and source:

```r
source("R/00_config.R")
source("R/01_generate_sample.R")
```

Run the generator checks with:

```r
source("tests/manual/test_sample_generator.R")
```

Inspect the installed package APIs with:

```r
source("tests/manual/inspect_censoring_apis.R")
```

The API inspection script reports installed versions and function arguments for NADA, NADA2, EnvStats and survival. It does not yet validate the statistical behaviour of those functions.

## Methodological cautions

1. Reverse Kaplan–Meier does not reveal the ordering or distributional shape of values below a common LOQ.
2. A package function returning a mean or confidence interval does not, by itself, establish that the quantity has the desired interpretation or coverage.
3. ROS, MLE and reverse KM do not necessarily return the same kind of output. Comparisons will distinguish distribution reconstruction, parameter estimation and inference.
4. Inferential procedures must be specified and validated by method.
5. Convergence failures, non-identifiable estimates and failed bootstrap replicates will be recorded, not silently discarded.

## Status

Initial repository scaffold. The protocol and package interfaces are subject to validation before the full Monte Carlo experiment is run.
