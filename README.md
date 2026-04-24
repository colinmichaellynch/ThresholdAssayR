# ThresholdAssayR

`ThresholdAssayR` provides simple R functions for designing, simulating, fitting, plotting, and powering binary response/satisfaction threshold assays.

The core model is

```text
\[
P(Y = 1 \mid x) = \operatorname{logistic}\{-k(x-\theta)\},
\]
```
where:

- `theta` is the cue level where the response probability is 0.5.
- `k` controls the direction and steepness of the curve.
- `k < 0` gives an increasing response-threshold curve.
- `k > 0` gives a decreasing satisfaction-threshold curve.

The functions are intended for studies where an investigator presents a set of cue or task-demand levels and observes whether an individual responds, stops, persists, or otherwise changes state.

## Files

Save both files in the same directory:

```text
threshold_assay_functions.R
tutorial_threshold_assay.R
```

Then run the tutorial from that directory:

```r
source("tutorial_threshold_assay.R")
```

For your own analysis, load the functions with:

```r
source("threshold_assay_functions.R")
```

## Dependencies

The fitting, simulation, and power-analysis functions use base R. Plotting functions use `ggplot2`:

```r
install.packages("ggplot2")
```

## Required data structure

Your dataset should have one row per cue presentation or trial. At minimum, it should include these columns:

| column | description |
|---|---|
| `presentation` | Trial or presentation order, usually `1, 2, ..., n`. This column is useful for documenting the experiment, although the model-fitting function only uses `cue` and `response`. |
| `cue` | Numeric cue or task-demand value presented on that trial. |
| `response` | Binary outcome coded as `0` or `1`. |

Example:

```r
pilot_data <- data.frame(
  presentation = 1:6,
  cue = c(0, 10, 20, 30, 40, 50),
  response = c(1, 1, 1, 0, 0, 0)
)
```

A CSV file should have the same structure:

```text
presentation,cue,response
1,0,1
2,10,1
3,20,1
4,30,0
5,40,0
6,50,0
```

You can load a real dataset with:

```r
pilot_data <- read.csv("my_real_threshold_assay_data.csv")
```

Then fit the model with:

```r
fit <- fit_threshold_model(
  data = pilot_data,
  cue_col = "cue",
  response_col = "response",
  alpha = 0.05
)
```

## Main workflow

### 1. Generate a sampling sequence

Four sampling strategies are implemented:

- `linear`
- `uniform_random`
- `chebyshev`
- `shifted_chebyshev`

```r
source("threshold_assay_functions.R")

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

The simulated pilot dataset is deliberately structured like a real dataset, with `presentation`, `cue`, and `response` columns.

```r
pilot_simulation <- simulate_threshold_assay(
  strategy = "linear",
  n = 16,
  theta = 25,
  k = 0.075,
  xmin = 0,
  xmax = 50,
  seed = 10
)

pilot_data <- pilot_simulation[, c("presentation", "cue", "response")]
head(pilot_data)
```

Here, `theta = 25` means the response probability is 0.5 when the cue is 25. The value `k = 0.075` produces a decreasing satisfaction-threshold curve.

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

The estimate of `theta` is obtained from the logistic coefficients using

```r
theta = -beta0 / beta1
```

and

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
  k = 0.075,
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
tutorial_threshold_assay.R
```

Run it from the directory containing both files:

```r
source("tutorial_threshold_assay.R")
```

The tutorial:

1. Describes the required real-data structure.
2. Generates the four sampling sequences.
3. Simulates a fake pilot study with the same columns as a real dataset.
4. Fits the threshold model.
5. Extracts parameter estimates and the likelihood-ratio test.
6. Plots the fitted model.
7. Runs a power analysis for linear spacing.
8. Plots uncertainty in recovered parameter estimates.
9. Compares power across a small grid of sample sizes.

## Function reference

### Design functions

```r
make_sampling_sequence()
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
