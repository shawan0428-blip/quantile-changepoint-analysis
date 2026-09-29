# Quantile Change-Point Analysis

R research code for change-point detection and quantile regression, using simulated time series and annual temperature data from Tuscaloosa, Alabama.

## Scripts

| File | Description |
| --- | --- |
| `simulation_pelt_quantile_regression.R` | Simulates a time series, detects change points, and fits quantile regressions within each segment. |
| `tuscaloosa_pelt_quantile_regression.R` | Applies change-point detection and segment-wise quantile regression to annual temperature data. |
| `Full_NPPELT.R` | Simulates a QAR process, runs a custom PELT implementation, and compares detected changes with EnvCpt. |

## Requirements

Install the R packages once:

```r
install.packages(c("quantreg", "EnvCpt", "changepoint"))
```

- The simulation and Tuscaloosa scripts require a separate definition of `PELT.quantloss()`. Load it first or save it in `PELT_quantloss.R` in the working directory.
- The Tuscaloosa script also requires `Tuscaloosa_annual.txt` in the working directory.
- `Full_NPPELT.R` defines its own detector, `pelt_nonparametric()`, which is a different function from `PELT.quantloss()`.

## How to run

Set the repository folder as your R working directory and run one script at a time. For example:

```r
source("simulation_pelt_quantile_regression.R")
```

These scripts are research implementations under development.
