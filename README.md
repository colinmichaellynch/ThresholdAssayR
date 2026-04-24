# ThresholdAssayR

`ThresholdAssayR` provides R functions for designing, simulating, fitting, plotting, and powering binary response/satisfaction threshold assays.

The core model is

```text
P(Y = 1 | x) = logistic{-k(x - theta)}
```

where:

- `theta` is the cue level where the response probability is 0.5.
- `k` controls the direction and steepness of the curve.
- `k < 0` gives an increasing response-threshold curve.
- `k > 0` gives a decreasing satisfaction-threshold curve.

The functions are intended for studies where an investigator presents a set of cue or task-demand levels and observes whether an individual responds, stops, persists, or otherwise changes state.

## Files

```text
thresholdAssayR/
├── R/
│   └── threshold_assay_functions.R
├── examples/
│   └── tutorial_threshold_assay.R
├── DESCRIPTION
├── NAMESPACE
└── README.md
```

## Installation/use from source

For simple use, clone or download this repository and source the functions:

```r
source("R/threshold_assay_functions.R")
```

The fitting and simulation functions use base R. Plotting functions use `ggplot2`:

```r
install.packages("ggplot2")
```

If you convert this repository into a formal R package, you can install it with:

```r
install.packages("devtools")
devtools::install_github("YOUR-USERNAME/thresholdAssayR")
library(thresholdAssayR)
```

## Main workflow

### 1. Generate a sampling sequence

Four sampling strategies are implemented:

- `linear`
- `uniform_random`
- `chebyshev`
- `shifted_chebyshev`

```r
source("R/threshold_assay_functions.R")

sequence <- make_sampling_sequence(
  strategy = "linear",
  n = 16,
  xmin = 0,
  xmax = 50,
  seed = 1
)

sequence
```

### 2. Simulate a pilot dataset

```r
pilot_data <- simulate_threshold_assay(
  strategy = "linear",
  n = 16,
  theta = 25,
  k = -0.25,
  xmin = 0,
  xmax = 50,
  seed = 10
)

head(pilot_data)
```

Here, `theta = 25` means the response probability is 0.5 when the cue is 25. The value `k = -0.25` produces an increasing response-threshold curve.

### 3. Fit the threshold model

```r
fit <- fit_threshold_model(
  data = pilot_data,
  cue_col = "cue",
  response_col = "response",
  alpha = 0.05
)

print(fit)
```

The fitted logistic model is compared with an intercept-only null model using a likelihood-ratio test.

### 4. Extract parameter estimates and uncertainty

```r
parameter_estimates(fit)
```

This returns estimates, standard errors, and confidence intervals for:

- `theta`
- `k`

The estimate of `theta` is obtained from the logistic coefficients using:

```r
theta = -beta0 / beta1
```

and:

```r
k = -beta1
```

Uncertainty in `theta` is calculated using the delta method.

### 5. Run the likelihood-ratio test

```r
likelihood_ratio_test(fit)
```

This tests whether the cue-dependent threshold model fits better than an intercept-only null model.

### 6. Plot the fitted model

```r
plot_threshold_fit(fit)
```

This plots the binary responses, the fitted response curve, and an approximate confidence band.

### 7. Estimate statistical power

```r
power <- estimate_threshold_power(
  alpha = 0.05,
  strategy = "linear",
  n = 16,
  theta = 25,
  k = -0.25,
  nsim = 1000,
  xmin = 0,
  xmax = 50,
  seed = 2026
)

print(power)
```

The power analysis simulates many assays under assumed values of `theta` and `k`, fits the model to each simulated dataset, and estimates power as the fraction of likelihood-ratio tests that reject the intercept-only null model.

The returned object includes:

```r
power$power_summary
power$parameter_summary
power$simulations
```

`parameter_summary` includes expected bias, RMSE, mean standard error, empirical intervals, and coverage for both `theta` and `k`.

### 8. Plot parameter uncertainty from the power analysis

```r
plot_power_parameter_uncertainty(power)
```

The error bars show empirical intervals from the simulated assays, and the x marks show the true assumed parameter values.

## Tutorial

A full tutorial is provided in:

```text
examples/tutorial_threshold_assay.R
```

Run it from the repository root:

```r
source("examples/tutorial_threshold_assay.R")
```

The tutorial:

1. Generates the four sampling sequences.
2. Simulates a fake pilot study.
3. Fits the threshold model.
4. Extracts parameter estimates and the likelihood-ratio test.
5. Plots the fitted model.
6. Runs a power analysis for linear spacing.
7. Plots uncertainty in recovered parameter estimates.
8. Compares power across a small grid of sample sizes.

## Function reference

### Design functions

```r
make_sampling_sequence()
plot_sampling_sequences()
```

### Simulation functions

```r
threshold_probability()
simulate_threshold_data()
simulate_threshold_assay()
```

### Model-fitting functions

```r
fit_threshold_model()
parameter_estimates()
likelihood_ratio_test()
predict_threshold_fit()
plot_threshold_fit()
```

### Power-analysis functions

```r
estimate_threshold_power()
plot_power_parameter_uncertainty()
```

## Notes on interpretation

A significant likelihood-ratio test indicates that the response probability changes with cue level relative to an intercept-only model. It does not, by itself, guarantee precise recovery of both parameters. For that reason, the power-analysis function also reports uncertainty and empirical recovery error for `theta` and `k`.

When the true curve is very shallow, when sample size is small, or when cue levels do not span the transition region, estimates of `theta` and `k` can be unstable. Pilot data should therefore be used both to test for a cue-dependent response and to plan follow-up sampling designs with adequate power and acceptable parameter uncertainty.

## Reproducibility

Most functions accept a `seed` argument. Set this argument to reproduce sampling sequences, simulated datasets, and power analyses.

## License

Add your preferred license here before uploading to GitHub.
