# Threshold assay design and analysis functions
#
# This file contains functions for designing, simulating, fitting, plotting,
# and conducting power analyses for binary response/satisfaction threshold
# assays. The model is
#
#   P(Y = 1 | x) = logistic[-k * (x - theta)]
#
# where theta is the cue level at which P(Y = 1) = 0.5 and k controls the
# direction and steepness of the cue-response curve. Negative k gives an
# increasing response-threshold curve; positive k gives a decreasing
# satisfaction-threshold curve.
#
# Expected data structure
#   The fitting functions expect one row per cue presentation/trial. At minimum,
#   the dataset should contain:
#
#     presentation : trial/presentation order, usually 1, 2, ..., n
#     cue          : numeric cue or task-demand value
#     response     : binary outcome coded 0/1
#
#   The presentation column is not required for fitting, but it is useful for
#   documenting the trial order and matching the real experimental data format.

#' Logistic threshold probability
#'
#' Computes P(Y = 1 | x) under the threshold model.
#'
#' @param x Numeric vector of cue values.
#' @param theta Threshold location.
#' @param k Direction/steepness parameter.
#'
#' @return Numeric vector of probabilities.
threshold_probability <- function(x, theta, k) {
  stats::plogis(-k * (x - theta))
}

#' Make a sampling sequence
#'
#' Generates cue levels for one of four sampling strategies.
#'
#' @param strategy One of "linear", "uniform_random", "chebyshev",
#'   or "shifted_chebyshev". Case and spaces are ignored.
#' @param n Number of observations/cue presentations.
#' @param xmin Minimum cue value.
#' @param xmax Maximum cue value.
#' @param randomize_order Logical. If TRUE, cue presentation order is randomized
#'   after cue levels are generated.
#' @param seed Optional random seed.
#'
#' @return A data.frame with presentation index and cue value.
make_sampling_sequence <- function(strategy,
                                   n,
                                   xmin = 0,
                                   xmax = 50,
                                   randomize_order = TRUE,
                                   seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  if (length(n) != 1 || n < 1 || n != as.integer(n)) {
    stop("n must be a positive integer.")
  }
  if (xmax <= xmin) {
    stop("xmax must be greater than xmin.")
  }

  strategy_clean <- tolower(gsub("[ -]", "_", strategy))

  if (strategy_clean == "linear") {
    cue <- seq(xmin, xmax, length.out = n)

  } else if (strategy_clean %in% c("uniform_random", "uniform")) {
    cue <- sort(stats::runif(n, xmin, xmax))

  } else if (strategy_clean == "chebyshev") {
    i <- seq_len(n)
    nodes <- cos((2 * i - 1) * pi / (2 * n))
    cue <- xmin + (xmax - xmin) * sort((nodes + 1) / 2)

  } else if (strategy_clean %in% c("shifted_chebyshev", "shifted")) {
    i <- seq_len(n)
    nodes <- cos((2 * i - 1) * pi / (2 * n))
    scaled_nodes <- (nodes + 1) / 2

    # Pull ordinary Chebyshev nodes toward the center. This keeps deterministic
    # boundary coverage while reducing the concentration of points at the
    # extremes.
    shifted_nodes <- 0.5 + 0.70 * (scaled_nodes - 0.5)
    cue <- xmin + (xmax - xmin) * sort(shifted_nodes)

  } else {
    stop(
      "Unknown strategy. Use one of: linear, uniform_random, chebyshev, shifted_chebyshev."
    )
  }

  if (randomize_order && n > 1) {
    cue <- sample(cue, size = n, replace = FALSE)
  }

  data.frame(
    presentation = seq_len(n),
    cue = as.numeric(cue)
  )
}

#' Simulate threshold-assay data
#'
#' Generates binary observations from the logistic threshold model.
#'
#' @param cue Numeric vector of cue values.
#' @param theta Threshold location.
#' @param k Direction/steepness parameter.
#' @param presentation Optional presentation/trial index. If NULL, this is
#'   set to 1:length(cue).
#' @param seed Optional random seed.
#'
#' @return A data.frame with presentation, cue, response, true probability,
#'   theta, and k.
simulate_threshold_data <- function(cue,
                                    theta,
                                    k,
                                    presentation = NULL,
                                    seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  if (is.null(presentation)) {
    presentation <- seq_along(cue)
  }

  p <- threshold_probability(cue, theta = theta, k = k)
  response <- stats::rbinom(length(cue), size = 1, prob = p)

  data.frame(
    presentation = as.integer(presentation),
    cue = as.numeric(cue),
    response = as.integer(response),
    probability = as.numeric(p),
    theta_true = theta,
    k_true = k
  )
}

#' Simulate an assay from a sampling strategy
#'
#' Convenience wrapper that creates a sampling sequence and then simulates
#' binary responses.
#'
#' @param strategy Sampling strategy.
#' @param n Number of observations.
#' @param theta Threshold location.
#' @param k Direction/steepness parameter.
#' @param xmin Minimum cue value.
#' @param xmax Maximum cue value.
#' @param randomize_order Logical; randomize the presentation order.
#' @param seed Optional random seed.
#'
#' @return A data.frame containing presentation, cue, response, probability,
#'   and true parameter values.
simulate_threshold_assay <- function(strategy,
                                     n,
                                     theta,
                                     k,
                                     xmin = 0,
                                     xmax = 50,
                                     randomize_order = TRUE,
                                     seed = NULL) {
  sequence <- make_sampling_sequence(
    strategy = strategy,
    n = n,
    xmin = xmin,
    xmax = xmax,
    randomize_order = randomize_order,
    seed = seed
  )

  # Use a second seed stream only if a seed was provided. This avoids making
  # the same random numbers drive both the design randomization and response
  # generation.
  response_seed <- if (is.null(seed)) NULL else seed + 1
  simulate_threshold_data(
    cue = sequence$cue,
    theta = theta,
    k = k,
    presentation = sequence$presentation,
    seed = response_seed
  )
}

#' Fit the logistic threshold model
#'
#' Fits the cue-dependent logistic threshold model and compares it with an
#' intercept-only null model using a likelihood-ratio test.
#'
#' @param data A data.frame containing one row per trial/presentation.
#' @param cue_col Name of cue column.
#' @param response_col Name of binary response column.
#' @param alpha Significance level used for confidence intervals.
#' @param conf_level Confidence level for Wald intervals.
#'
#' @return A list containing fitted models, parameter estimates, LRT results,
#'   and original data used for fitting.
fit_threshold_model <- function(data,
                                cue_col = "cue",
                                response_col = "response",
                                alpha = 0.05,
                                conf_level = 1 - alpha) {
  if (!cue_col %in% names(data)) stop("cue_col not found in data.")
  if (!response_col %in% names(data)) stop("response_col not found in data.")

  d <- data.frame(
    cue = as.numeric(data[[cue_col]]),
    response = as.integer(data[[response_col]])
  )

  if (!all(d$response %in% c(0, 1))) {
    stop("response column must contain only 0/1 values.")
  }

  full_model <- stats::glm(
    response ~ cue,
    family = stats::binomial(link = "logit"),
    data = d
  )

  null_model <- stats::glm(
    response ~ 1,
    family = stats::binomial(link = "logit"),
    data = d
  )

  beta_hat <- stats::coef(full_model)
  vc <- stats::vcov(full_model)

  beta0 <- unname(beta_hat[1])
  beta1 <- unname(beta_hat[2])

  theta_hat <- NA_real_
  theta_se <- NA_real_
  theta_ci <- c(NA_real_, NA_real_)

  k_hat <- -beta1
  k_se <- sqrt(vc[2, 2])

  zcrit <- stats::qnorm(1 - (1 - conf_level) / 2)

  if (is.finite(beta1) && abs(beta1) > .Machine$double.eps) {
    theta_hat <- -beta0 / beta1

    # Delta-method gradient for theta = -beta0 / beta1.
    grad_theta <- c(-1 / beta1, beta0 / beta1^2)
    theta_var <- as.numeric(t(grad_theta) %*% vc %*% grad_theta)
    theta_se <- sqrt(max(theta_var, 0))

    theta_ci <- theta_hat + c(-1, 1) * zcrit * theta_se
  }

  k_ci <- k_hat + c(-1, 1) * zcrit * k_se

  estimates <- data.frame(
    parameter = c("theta", "k"),
    estimate = c(theta_hat, k_hat),
    std_error = c(theta_se, k_se),
    conf_low = c(theta_ci[1], k_ci[1]),
    conf_high = c(theta_ci[2], k_ci[2]),
    row.names = NULL
  )

  lrt_stat <- as.numeric(2 * (stats::logLik(full_model) - stats::logLik(null_model)))
  lrt_df <- attr(stats::logLik(full_model), "df") - attr(stats::logLik(null_model), "df")
  lrt_p <- stats::pchisq(lrt_stat, df = lrt_df, lower.tail = FALSE)

  lrt <- data.frame(
    test = "full logistic threshold model vs intercept-only null model",
    statistic = lrt_stat,
    df = lrt_df,
    p_value = lrt_p,
    alpha = alpha,
    reject_null = lrt_p < alpha
  )

  out <- list(
    data = d,
    full_model = full_model,
    null_model = null_model,
    estimates = estimates,
    lrt = lrt,
    alpha = alpha,
    conf_level = conf_level
  )

  class(out) <- "threshold_fit"
  out
}

#' Extract threshold parameter estimates
#'
#' @param fit Object returned by fit_threshold_model().
#'
#' @return A data.frame of parameter estimates, standard errors, and intervals.
parameter_estimates <- function(fit) {
  if (!inherits(fit, "threshold_fit")) {
    stop("fit must be an object returned by fit_threshold_model().")
  }
  fit$estimates
}

#' Extract likelihood-ratio test results
#'
#' @param fit Object returned by fit_threshold_model().
#'
#' @return A data.frame containing the LRT result.
likelihood_ratio_test <- function(fit) {
  if (!inherits(fit, "threshold_fit")) {
    stop("fit must be an object returned by fit_threshold_model().")
  }
  fit$lrt
}

#' Predict from a threshold fit
#'
#' @param fit Object returned by fit_threshold_model().
#' @param new_cue Numeric vector of cue values.
#' @param conf_level Confidence level for prediction bands on the response scale.
#'
#' @return A data.frame with cue, fitted probability, and confidence limits.
predict_threshold_fit <- function(fit,
                                  new_cue = NULL,
                                  conf_level = fit$conf_level) {
  if (!inherits(fit, "threshold_fit")) {
    stop("fit must be an object returned by fit_threshold_model().")
  }

  if (is.null(new_cue)) {
    cue_range <- range(fit$data$cue, finite = TRUE)
    new_cue <- seq(cue_range[1], cue_range[2], length.out = 200)
  }

  newdata <- data.frame(cue = as.numeric(new_cue))
  pred <- stats::predict(
    fit$full_model,
    newdata = newdata,
    type = "link",
    se.fit = TRUE
  )

  zcrit <- stats::qnorm(1 - (1 - conf_level) / 2)

  eta <- pred$fit
  eta_low <- eta - zcrit * pred$se.fit
  eta_high <- eta + zcrit * pred$se.fit

  data.frame(
    cue = newdata$cue,
    probability = stats::plogis(eta),
    conf_low = stats::plogis(eta_low),
    conf_high = stats::plogis(eta_high)
  )
}

#' Plot the fitted threshold model
#'
#' @param fit Object returned by fit_threshold_model().
#' @param show_ci Logical; show confidence band for fitted probability.
#' @param jitter_height Vertical jitter applied to binary responses.
#'
#' @return A ggplot object if ggplot2 is installed.
plot_threshold_fit <- function(fit,
                               show_ci = TRUE,
                               jitter_height = 0.035) {
  if (!inherits(fit, "threshold_fit")) {
    stop("fit must be an object returned by fit_threshold_model().")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot_threshold_fit() requires ggplot2. Install it with install.packages('ggplot2').")
  }

  ggplot2 <- asNamespace("ggplot2")
  pred <- predict_threshold_fit(fit)

  p <- ggplot2$ggplot() +
    ggplot2$geom_point(
      data = fit$data,
      ggplot2$aes(x = cue, y = response),
      position = ggplot2$position_jitter(width = 0, height = jitter_height),
      alpha = 0.75
    )

  if (show_ci) {
    p <- p +
      ggplot2$geom_ribbon(
        data = pred,
        ggplot2$aes(x = cue, ymin = conf_low, ymax = conf_high),
        alpha = 0.20
      )
  }

  p +
    ggplot2$geom_line(
      data = pred,
      ggplot2$aes(x = cue, y = probability),
      linewidth = 1
    ) +
    ggplot2$coord_cartesian(ylim = c(0, 1)) +
    ggplot2$labs(
      x = "Cue level",
      y = "P(response)",
      title = "Fitted logistic threshold model"
    ) +
    ggplot2$theme_classic()
}

#' Estimate statistical power for a threshold assay
#'
#' Simulates many assays under assumed theta and k values, fits the threshold
#' model to each simulated assay, and estimates power as the probability of
#' rejecting the intercept-only null model with a likelihood-ratio test.
#'
#' @param alpha Significance level for the LRT.
#' @param strategy Sampling strategy.
#' @param n Sample size.
#' @param theta Assumed threshold location.
#' @param k Assumed direction/steepness parameter.
#' @param nsim Number of simulated assays.
#' @param xmin Minimum cue value.
#' @param xmax Maximum cue value.
#' @param randomize_order Logical; randomize presentation order.
#' @param seed Optional random seed.
#' @param conf_level Confidence level for parameter intervals.
#'
#' @return A list with settings, simulation-level results, and summary tables.
estimate_threshold_power <- function(alpha,
                                     strategy,
                                     n,
                                     theta,
                                     k,
                                     nsim = 1000,
                                     xmin = 0,
                                     xmax = 50,
                                     randomize_order = TRUE,
                                     seed = NULL,
                                     conf_level = 1 - alpha) {
  if (!is.null(seed)) set.seed(seed)

  sim_results <- vector("list", nsim)

  for (i in seq_len(nsim)) {
    # For uniform-random designs, this generates a new random design in each
    # simulation. For deterministic designs, it preserves the cue set but can
    # randomize the presentation order.
    sequence <- make_sampling_sequence(
      strategy = strategy,
      n = n,
      xmin = xmin,
      xmax = xmax,
      randomize_order = randomize_order
    )

    dat <- simulate_threshold_data(
      cue = sequence$cue,
      theta = theta,
      k = k,
      presentation = sequence$presentation
    )

    fit <- try(
      fit_threshold_model(dat, alpha = alpha, conf_level = conf_level),
      silent = TRUE
    )

    if (inherits(fit, "try-error")) {
      sim_results[[i]] <- data.frame(
        simulation = i,
        detected = NA,
        lrt_p = NA_real_,
        theta_hat = NA_real_,
        theta_se = NA_real_,
        theta_low = NA_real_,
        theta_high = NA_real_,
        theta_covered = NA,
        k_hat = NA_real_,
        k_se = NA_real_,
        k_low = NA_real_,
        k_high = NA_real_,
        k_covered = NA
      )
    } else {
      est <- parameter_estimates(fit)
      lrt <- likelihood_ratio_test(fit)

      theta_row <- est[est$parameter == "theta", ]
      k_row <- est[est$parameter == "k", ]

      sim_results[[i]] <- data.frame(
        simulation = i,
        detected = lrt$reject_null,
        lrt_p = lrt$p_value,
        theta_hat = theta_row$estimate,
        theta_se = theta_row$std_error,
        theta_low = theta_row$conf_low,
        theta_high = theta_row$conf_high,
        theta_covered = theta >= theta_row$conf_low && theta <= theta_row$conf_high,
        k_hat = k_row$estimate,
        k_se = k_row$std_error,
        k_low = k_row$conf_low,
        k_high = k_row$conf_high,
        k_covered = k >= k_row$conf_low && k <= k_row$conf_high
      )
    }
  }

  simulations <- do.call(rbind, sim_results)
  valid <- !is.na(simulations$detected)

  power_hat <- mean(simulations$detected[valid])
  power_ci <- stats::binom.test(
    sum(simulations$detected[valid]),
    sum(valid),
    conf.level = conf_level
  )$conf.int

  parameter_summary <- data.frame(
    parameter = c("theta", "k"),
    true_value = c(theta, k),
    mean_estimate = c(
      mean(simulations$theta_hat, na.rm = TRUE),
      mean(simulations$k_hat, na.rm = TRUE)
    ),
    median_estimate = c(
      stats::median(simulations$theta_hat, na.rm = TRUE),
      stats::median(simulations$k_hat, na.rm = TRUE)
    ),
    bias = c(
      mean(simulations$theta_hat - theta, na.rm = TRUE),
      mean(simulations$k_hat - k, na.rm = TRUE)
    ),
    rmse = c(
      sqrt(mean((simulations$theta_hat - theta)^2, na.rm = TRUE)),
      sqrt(mean((simulations$k_hat - k)^2, na.rm = TRUE))
    ),
    mean_std_error = c(
      mean(simulations$theta_se, na.rm = TRUE),
      mean(simulations$k_se, na.rm = TRUE)
    ),
    empirical_low = c(
      stats::quantile(simulations$theta_hat, probs = (1 - conf_level) / 2, na.rm = TRUE),
      stats::quantile(simulations$k_hat, probs = (1 - conf_level) / 2, na.rm = TRUE)
    ),
    empirical_high = c(
      stats::quantile(simulations$theta_hat, probs = 1 - (1 - conf_level) / 2, na.rm = TRUE),
      stats::quantile(simulations$k_hat, probs = 1 - (1 - conf_level) / 2, na.rm = TRUE)
    ),
    coverage = c(
      mean(simulations$theta_covered, na.rm = TRUE),
      mean(simulations$k_covered, na.rm = TRUE)
    ),
    row.names = NULL
  )

  power_summary <- data.frame(
    alpha = alpha,
    conf_level = conf_level,
    strategy = tolower(gsub("[ -]", "_", strategy)),
    n = n,
    theta = theta,
    k = k,
    nsim_requested = nsim,
    nsim_successful = sum(valid),
    power = power_hat,
    power_conf_low = power_ci[1],
    power_conf_high = power_ci[2]
  )

  out <- list(
    settings = list(
      alpha = alpha,
      conf_level = conf_level,
      strategy = strategy,
      n = n,
      theta = theta,
      k = k,
      nsim = nsim,
      xmin = xmin,
      xmax = xmax,
      randomize_order = randomize_order,
      seed = seed
    ),
    power_summary = power_summary,
    parameter_summary = parameter_summary,
    simulations = simulations
  )

  class(out) <- "threshold_power"
  out
}

#' Plot parameter-estimate uncertainty from a power analysis
#'
#' Shows the mean estimate for theta and k with empirical error bars from the
#' simulated assays. True parameter values are shown as points.
#'
#' @param power Object returned by estimate_threshold_power().
#'
#' @return A ggplot object if ggplot2 is installed.
plot_power_parameter_uncertainty <- function(power) {
  if (!inherits(power, "threshold_power")) {
    stop("power must be an object returned by estimate_threshold_power().")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot_power_parameter_uncertainty() requires ggplot2.")
  }

  ggplot2 <- asNamespace("ggplot2")

  df <- power$parameter_summary

  ggplot2$ggplot(df, ggplot2$aes(x = parameter, y = mean_estimate)) +
    ggplot2$geom_errorbar(
      ggplot2$aes(ymin = empirical_low, ymax = empirical_high),
      width = 0.15
    ) +
    ggplot2$geom_point(size = 2) +
    ggplot2$geom_point(
      ggplot2$aes(y = true_value),
      shape = 4,
      size = 3,
      stroke = 1
    ) +
    ggplot2$labs(
      x = "Parameter",
      y = "Estimate",
      title = "Parameter recovery across simulated assays",
      subtitle = "Error bars show empirical intervals; x marks true parameter values"
    ) +
    ggplot2$theme_classic()
}

#' Print method for fitted threshold models
#'
#' @param x Object returned by fit_threshold_model().
#' @param ... Ignored.
print.threshold_fit <- function(x, ...) {
  cat("Threshold model fit\n")
  cat("===================\n\n")
  print(x$estimates)
  cat("\nLikelihood-ratio test:\n")
  print(x$lrt)
  invisible(x)
}

#' Print method for threshold power analyses
#'
#' @param x Object returned by estimate_threshold_power().
#' @param ... Ignored.
print.threshold_power <- function(x, ...) {
  cat("Threshold assay power analysis\n")
  cat("==============================\n\n")
  print(x$power_summary)
  cat("\nParameter recovery summary:\n")
  print(x$parameter_summary)
  invisible(x)
}
