# Tutorial: fitting and powering a threshold assay
#
# This tutorial shows how to use the threshold-assay functions on a simulated
# pilot study. The pilot study is structured the same way a real dataset should
# be structured: one row per presentation/trial, with columns for presentation
# order, cue level, and binary response.

# Save this tutorial file and threshold_assay_functions.R in the same directory.
source("threshold_assay_functions.R")

# ggplot2 is only needed for plotting.
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

set.seed(101)

# -------------------------------------------------------------------------
# 1. Expected data structure
# -------------------------------------------------------------------------

# Real datasets should have at least these columns:
#
#   presentation : the trial or presentation order, usually 1, 2, ..., n
#   cue          : the cue or task-demand value presented on that trial
#   response     : binary response coded as 0/1
#
# The fitting function uses only cue and response, but presentation is useful
# for documenting the actual experimental order.

# -------------------------------------------------------------------------
# 2. Define the assumed pilot-study parameters
# -------------------------------------------------------------------------

theta_pilot <- 25     # threshold: P(response) = 0.5 at cue = 25
k_pilot <- 0.075      # positive k gives a decreasing satisfaction-threshold curve
n_pilot <- 16
xmin <- 0
xmax <- 50

# -------------------------------------------------------------------------
# 3. Generate sampling sequences for the four candidate designs
# -------------------------------------------------------------------------

linear_sequence <- make_sampling_sequence(
  strategy = "linear",
  n = n_pilot,
  xmin = xmin,
  xmax = xmax,
  seed = 1
)

uniform_sequence <- make_sampling_sequence(
  strategy = "uniform_random",
  n = n_pilot,
  xmin = xmin,
  xmax = xmax,
  seed = 2
)

chebyshev_sequence <- make_sampling_sequence(
  strategy = "chebyshev",
  n = n_pilot,
  xmin = xmin,
  xmax = xmax,
  seed = 3
)

shifted_chebyshev_sequence <- make_sampling_sequence(
  strategy = "shifted_chebyshev",
  n = n_pilot,
  xmin = xmin,
  xmax = xmax,
  seed = 4
)

print(linear_sequence)
print(uniform_sequence)
print(chebyshev_sequence)
print(shifted_chebyshev_sequence)

# Optional plot of the four sampling sequences.
if (requireNamespace("ggplot2", quietly = TRUE)) {
  all_sequences <- rbind(
    data.frame(linear_sequence, strategy = "linear"),
    data.frame(uniform_sequence, strategy = "uniform_random"),
    data.frame(chebyshev_sequence, strategy = "chebyshev"),
    data.frame(shifted_chebyshev_sequence, strategy = "shifted_chebyshev")
  )

  ggplot2::ggplot(all_sequences, ggplot2::aes(x = presentation, y = cue)) +
    ggplot2::geom_point() +
    ggplot2::geom_line(alpha = 0.5) +
    ggplot2::facet_wrap(~ strategy) +
    ggplot2::labs(
      x = "Presentation index",
      y = "Cue level",
      title = paste("Sampling sequences for n =", n_pilot)
    ) +
    ggplot2::theme_classic()
}

# -------------------------------------------------------------------------
# 4. Simulate a pilot dataset using linear spacing
# -------------------------------------------------------------------------

# This simulated pilot_data has the same required columns as a real dataset:
# presentation, cue, and response.
pilot_simulation <- simulate_threshold_assay(
  strategy = "linear",
  n = n_pilot,
  theta = theta_pilot,
  k = k_pilot,
  xmin = xmin,
  xmax = xmax,
  seed = 10
)

pilot_data <- pilot_simulation[, c("presentation", "cue", "response")]
print(pilot_data)

# If you had a real dataset in a CSV file, it should look like pilot_data and
# could be loaded like this:
#
# pilot_data <- read.csv("my_real_threshold_assay_data.csv")
#
# The required columns are:
#
#   presentation,cue,response
#   1,0,1
#   2,3.33,1
#   3,6.67,1
#   ...
#
# response must be coded as 0/1.

# -------------------------------------------------------------------------
# 5. Fit the threshold model to the pilot data
# -------------------------------------------------------------------------

pilot_fit <- fit_threshold_model(
  data = pilot_data,
  cue_col = "cue",
  response_col = "response",
  alpha = 0.05
)

print(pilot_fit)

# Extract parameter estimates.
parameter_estimates(pilot_fit)

# Extract likelihood-ratio test comparing cue-dependent model to null model.
likelihood_ratio_test(pilot_fit)

# Plot the fitted response curve.
plot_threshold_fit(pilot_fit)

# -------------------------------------------------------------------------
# 6. Run a power analysis for a planned linear-spacing assay
# -------------------------------------------------------------------------

# For a fast tutorial, use nsim = 200. For a final analysis, increase nsim to
# 1000 or more.
linear_power <- estimate_threshold_power(
  alpha = 0.05,
  strategy = "linear",
  n = n_pilot,
  theta = theta_pilot,
  k = k_pilot,
  nsim = 200,
  xmin = xmin,
  xmax = xmax,
  seed = 2026
)

print(linear_power)

# The power summary gives the estimated probability that the likelihood-ratio
# test rejects the intercept-only null model.
linear_power$power_summary

# The parameter summary gives expected error, uncertainty, and coverage for
# theta and k.
linear_power$parameter_summary

# Plot empirical uncertainty in the recovered parameter estimates.
plot_power_parameter_uncertainty(linear_power)

# -------------------------------------------------------------------------
# 7. Optional: compare power across multiple sample sizes for linear spacing
# -------------------------------------------------------------------------

n_grid <- c(8, 12, 16, 24, 32)

power_by_n <- lapply(n_grid, function(n_now) {
  estimate_threshold_power(
    alpha = 0.05,
    strategy = "linear",
    n = n_now,
    theta = theta_pilot,
    k = k_pilot,
    nsim = 200,
    xmin = xmin,
    xmax = xmax,
    seed = 5000 + n_now
  )$power_summary
})

power_by_n <- do.call(rbind, power_by_n)
print(power_by_n)

if (requireNamespace("ggplot2", quietly = TRUE)) {
  ggplot2::ggplot(power_by_n, ggplot2::aes(x = n, y = power)) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = power_conf_low, ymax = power_conf_high),
      width = 0.5
    ) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(
      x = "Sample size",
      y = "Estimated power",
      title = "Power for linear-spacing threshold assays"
    ) +
    ggplot2::theme_classic()
}

# -------------------------------------------------------------------------
# 8. Optional: save key tutorial outputs
# -------------------------------------------------------------------------

# Uncomment these lines if you want to save tutorial output tables.
#
# if (!dir.exists("tutorial_outputs")) dir.create("tutorial_outputs")
#
# write.csv(pilot_data, "tutorial_outputs/pilot_data.csv", row.names = FALSE)
# write.csv(parameter_estimates(pilot_fit), "tutorial_outputs/pilot_parameter_estimates.csv", row.names = FALSE)
# write.csv(likelihood_ratio_test(pilot_fit), "tutorial_outputs/pilot_lrt.csv", row.names = FALSE)
# write.csv(linear_power$power_summary, "tutorial_outputs/linear_power_summary.csv", row.names = FALSE)
# write.csv(linear_power$parameter_summary, "tutorial_outputs/linear_power_parameter_summary.csv", row.names = FALSE)
# write.csv(power_by_n, "tutorial_outputs/linear_power_by_sample_size.csv", row.names = FALSE)
