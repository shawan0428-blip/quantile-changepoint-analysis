# Quantile Change-Point Analysis

An R example that simulates a time series with a change after observation 50, detects change points using PELT, and fits quantile regressions at 0.25, 0.50, and 0.75 within each segment.

The analysis first removes a linear trend and adds it back to the fitted quantiles.

## How to run

Requires R, `quantreg`, and a separate definition of `PELT.quantloss()`.

Save the detector function in `PELT_quantloss.R` alongside the main script. Set that folder as your R working directory, then run:

```r
install.packages("quantreg")  # Only needed once
source("simulation_pelt_quantile_regression.R")
```

## Output

The script prints the detected change points and produces four plots: one change-point detection plot and three quantile regression plots.
