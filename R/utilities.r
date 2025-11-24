#' Utility Functions for powerprior Package
#'
#' Additional helper functions for power prior analysis

#' @keywords internal
#' @importFrom rlang .data
"_PACKAGE"

#' Calculate Effective Sample Size
#'
#' @description
#' Computes the effective sample size from historical data given a
#' discounting parameter.
#'
#' @param m Sample size of historical data
#' @param a0 Discounting parameter
#'
#' @return Effective sample size (numeric)
#'
#' @examples
#' # With 100 historical observations and 50% discounting
#' effective_sample_size(100, a0 = 0.5)  # Returns 50
#'
#' @export
effective_sample_size <- function(m, a0) {
  if (!is.numeric(m) || m <= 0) {
    stop("m must be a positive number")
  }
  if (!is.numeric(a0) || any(a0 < 0) || any(a0 > 1)) {
    stop("a0 must be between 0 and 1")
  }
  return(a0 * m)
}


#' Compare Multiple Discounting Parameters
#'
#' @description
#' Computes posteriors for multiple values of a0 to facilitate sensitivity analysis.
#'
#' @param historical_data Historical data (vector for univariate, matrix for multivariate)
#' @param current_data Current data (vector for univariate, matrix for multivariate)
#' @param a0_values Vector of discounting parameters to compare
#' @param multivariate Logical indicating if data is multivariate (default: FALSE)
#' @param n_samples Number of posterior samples per a0 value (default: 1000)
#' @param ... Additional arguments passed to powerprior functions
#'
#' @return A data frame with summary statistics for each a0 value
#'
#' @examples
#' set.seed(123)
#' historical <- rnorm(50, mean = 10, sd = 2)
#' current <- rnorm(30, mean = 10.5, sd = 2)
#'
#' comparison <- compare_discounting(
#'   historical,
#'   current,
#'   a0_values = c(0, 0.25, 0.5, 0.75, 1.0)
#' )
#' print(comparison)
#'
#' @export
compare_discounting <- function(historical_data, current_data,
                                a0_values = seq(0, 1, by = 0.1),
                                multivariate = FALSE,
                                n_samples = 1000,
                                ...) {

  if (!is.numeric(a0_values) || any(a0_values < 0) || any(a0_values > 1)) {
    stop("a0_values must be numeric values between 0 and 1")
  }

  n_a0 <- length(a0_values)

  if (multivariate) {
    # Multivariate case
    p <- ncol(historical_data)

    results <- data.frame(
      a0 = a0_values,
      ess = effective_sample_size(nrow(historical_data), a0_values)
    )

    # Add columns for each dimension
    for (d in 1:p) {
      results[[paste0("mu", d, "_mean")]] <- NA
      results[[paste0("mu", d, "_sd")]] <- NA
      results[[paste0("mu", d, "_ci_lower")]] <- NA
      results[[paste0("mu", d, "_ci_upper")]] <- NA
    }

    for (i in seq_along(a0_values)) {
      pp <- powerprior_multivariate(historical_data, a0 = a0_values[i], ...)
      posterior <- posterior_multivariate(pp, current_data)
      samples <- sample_posterior_multivariate(posterior, n_samples = n_samples,
                                                marginal = TRUE)

      for (d in 1:p) {
        results[[paste0("mu", d, "_mean")]][i] <- mean(samples[, d])
        results[[paste0("mu", d, "_sd")]][i] <- sd(samples[, d])
        results[[paste0("mu", d, "_ci_lower")]][i] <- quantile(samples[, d], 0.025)
        results[[paste0("mu", d, "_ci_upper")]][i] <- quantile(samples[, d], 0.975)
      }
    }

  } else {
    # Univariate case
    results <- data.frame(
      a0 = a0_values,
      ess = effective_sample_size(length(historical_data), a0_values),
      mu_mean = NA,
      mu_sd = NA,
      mu_ci_lower = NA,
      mu_ci_upper = NA,
      sigma2_mean = NA,
      sigma2_sd = NA
    )

    for (i in seq_along(a0_values)) {
      pp <- powerprior_univariate(historical_data, a0 = a0_values[i], ...)
      posterior <- posterior_univariate(pp, current_data)
      samples <- sample_posterior_univariate(posterior, n_samples = n_samples)

      results$mu_mean[i] <- mean(samples[, "mu"])
      results$mu_sd[i] <- sd(samples[, "mu"])
      results$mu_ci_lower[i] <- quantile(samples[, "mu"], 0.025)
      results$mu_ci_upper[i] <- quantile(samples[, "mu"], 0.975)
      results$sigma2_mean[i] <- mean(samples[, "sigma2"])
      results$sigma2_sd[i] <- sd(samples[, "sigma2"])
    }
  }

  class(results) <- c("powerprior_comparison", "data.frame")
  return(results)
}


#' Print method for powerprior_comparison
#' @param x Object of class "powerprior_comparison"
#' @param digits Number of digits to round to
#' @param ... Additional arguments
#' @return Invisibly returns the input object (for method chaining)
#' @export
print.powerprior_comparison <- function(x, digits = 4, ...) {
  cat("Power Prior Sensitivity Analysis\n")
  cat("=================================\n\n")
  cat("Comparing", nrow(x), "discounting parameters\n\n")

  # Round numeric columns
  numeric_cols <- sapply(x, is.numeric)
  x[numeric_cols] <- lapply(x[numeric_cols], round, digits = digits)

  print(as.data.frame(x), row.names = FALSE)
  invisible(x)
}


#' Plot Sensitivity Analysis
#'
#' @description
#' Creates a plot showing how posterior estimates vary with the discounting parameter.
#'
#' @param comparison_results Output from compare_discounting()
#' @param parameter Name of parameter to plot (default: "mu")
#' @param dimension For multivariate case, which dimension to plot (default: 1)
#'
#' @return A ggplot2 object (if ggplot2 is available) or base R plot
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' historical <- rnorm(50, mean = 10, sd = 2)
#' current <- rnorm(30, mean = 10.5, sd = 2)
#'
#' comparison <- compare_discounting(historical, current)
#' plot_sensitivity(comparison)
#' }
#'
#' @export
plot_sensitivity <- function(comparison_results, parameter = "mu", dimension = 1) {

  if (!inherits(comparison_results, "powerprior_comparison")) {
    stop("comparison_results must be from compare_discounting()")
  }

  # Check if data is multivariate
  is_multivariate <- any(grepl("mu1_", names(comparison_results)))

  if (is_multivariate) {
    mean_col <- paste0("mu", dimension, "_mean")
    ci_lower_col <- paste0("mu", dimension, "_ci_lower")
    ci_upper_col <- paste0("mu", dimension, "_ci_upper")
    title <- paste0("Sensitivity Analysis: mu[", dimension, "] vs Discounting Parameter")
    ylab <- paste0("Posterior mu[", dimension, "]")
  } else {
    if (parameter == "mu") {
      mean_col <- "mu_mean"
      ci_lower_col <- "mu_ci_lower"
      ci_upper_col <- "mu_ci_upper"
      title <- "Sensitivity Analysis: mu vs Discounting Parameter"
      ylab <- "Posterior mu"
    } else if (parameter == "sigma2") {
      mean_col <- "sigma2_mean"
      ci_lower_col <- NULL
      ci_upper_col <- NULL
      title <- "Sensitivity Analysis: sigma2 vs Discounting Parameter"
      ylab <- "Posterior sigma2"
    } else {
      stop("parameter must be 'mu' or 'sigma2'")
    }
  }

  # Try to use ggplot2 if available
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    p <- ggplot2::ggplot(comparison_results, ggplot2::aes(x = .data$a0, y = .data[[!!mean_col]])) +
      ggplot2::geom_line(size = 1, color = "blue") +
      ggplot2::geom_point(size = 2, color = "blue")

    if (!is.null(ci_lower_col) && !is.null(ci_upper_col)) {
      p <- p + ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data[[ci_lower_col]], ymax = .data[[ci_upper_col]]),
        alpha = 0.2, fill = "blue"
      )
    }

    p <- p +
      ggplot2::labs(
        title = title,
        x = "Discounting Parameter (a0)",
        y = ylab
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
        axis.text = ggplot2::element_text(size = 10),
        axis.title = ggplot2::element_text(size = 12)
      )

    return(p)

  } else {
    # Base R plot
    plot(comparison_results$a0, comparison_results[[mean_col]],
         type = "b", pch = 19, col = "blue",
         xlab = "Discounting Parameter (a0)",
         ylab = ylab,
         main = title,
         ylim = if (!is.null(ci_lower_col)) {
           range(c(comparison_results[[ci_lower_col]],
                   comparison_results[[ci_upper_col]]))
         } else {
           NULL
         })

    if (!is.null(ci_lower_col) && !is.null(ci_upper_col)) {
      polygon(c(comparison_results$a0, rev(comparison_results$a0)),
              c(comparison_results[[ci_lower_col]],
                rev(comparison_results[[ci_upper_col]])),
              col = rgb(0, 0, 1, 0.2), border = NA)
    }

    grid()
  }
}


#' Compute Posterior Summaries
#'
#' @description
#' Provides comprehensive posterior summary statistics.
#'
#' @param posterior Posterior object from posterior_univariate() or
#'   posterior_multivariate()
#' @param n_samples Number of samples to draw (default: 10000)
#' @param prob Probability for credible intervals (default: 0.95)
#'
#' @return A list with posterior summaries
#'
#' @examples
#' set.seed(123)
#' historical <- rnorm(50, mean = 10, sd = 2)
#' current <- rnorm(30, mean = 10.5, sd = 2)
#'
#' pp <- powerprior_univariate(historical, a0 = 0.5)
#' posterior <- posterior_univariate(pp, current)
#' summary <- posterior_summary(posterior)
#' print(summary)
#'
#' @export
posterior_summary <- function(posterior, n_samples = 10000, prob = 0.95) {

  alpha <- 1 - prob

  if (inherits(posterior, "posterior_univariate")) {
    # Univariate case
    samples <- sample_posterior_univariate(posterior, n_samples = n_samples)

    mu_summary <- c(
      mean = mean(samples[, "mu"]),
      sd = sd(samples[, "mu"]),
      median = median(samples[, "mu"]),
      quantile(samples[, "mu"], probs = c(alpha/2, 1 - alpha/2))
    )
    names(mu_summary)[4:5] <- c(paste0(alpha/2*100, "%"),
                                 paste0((1-alpha/2)*100, "%"))

    sigma2_summary <- c(
      mean = mean(samples[, "sigma2"]),
      sd = sd(samples[, "sigma2"]),
      median = median(samples[, "sigma2"]),
      quantile(samples[, "sigma2"], probs = c(alpha/2, 1 - alpha/2))
    )
    names(sigma2_summary)[4:5] <- c(paste0(alpha/2*100, "%"),
                                     paste0((1-alpha/2)*100, "%"))

    result <- list(
      mu = mu_summary,
      sigma2 = sigma2_summary,
      prob = prob,
      n_samples = n_samples
    )

  } else if (inherits(posterior, "posterior_multivariate")) {
    # Multivariate case
    samples <- sample_posterior_multivariate(posterior, n_samples = n_samples,
                                              marginal = TRUE)
    p <- ncol(samples)

    mu_summary <- matrix(NA, nrow = p, ncol = 5)
    colnames(mu_summary) <- c("mean", "sd", "median",
                               paste0(alpha/2*100, "%"),
                               paste0((1-alpha/2)*100, "%"))
    rownames(mu_summary) <- paste0("mu[", 1:p, "]")

    for (d in 1:p) {
      mu_summary[d, ] <- c(
        mean = mean(samples[, d]),
        sd = sd(samples[, d]),
        median = median(samples[, d]),
        quantile(samples[, d], probs = c(alpha/2, 1 - alpha/2))
      )
    }

    result <- list(
      mu = mu_summary,
      prob = prob,
      n_samples = n_samples,
      dimension = p
    )

  } else {
    stop("posterior must be from posterior_univariate() or posterior_multivariate()")
  }

  class(result) <- c("powerprior_summary", "list")
  return(result)
}


#' Print method for powerprior_summary
#' @param x Object of class "powerprior_summary"
#' @param digits Number of digits to round to
#' @param ... Additional arguments
#' @return Invisibly returns the input object (for method chaining)
#' @export
print.powerprior_summary <- function(x, digits = 4, ...) {
  cat("Posterior Summary Statistics\n")
  cat("============================\n\n")
  cat("Number of samples:", x$n_samples, "\n")
  cat("Credible interval probability:", x$prob, "\n\n")

  if (!is.null(x$mu) && is.matrix(x$mu)) {
    cat("Multivariate Mean Parameter (mu):\n")
    print(round(x$mu, digits))
  } else {
    cat("Mean Parameter (mu):\n")
    print(round(x$mu, digits))

    if (!is.null(x$sigma2)) {
      cat("\nVariance Parameter (sigma2):\n")
      print(round(x$sigma2, digits))
    }
  }

  invisible(x)
}


#' Calculate Bayes Factor
#'
#' @description
#' Computes the Bayes factor comparing two models with different discounting parameters.
#'
#' @param historical_data Historical data
#' @param current_data Current data
#' @param a0_1 First discounting parameter
#' @param a0_2 Second discounting parameter (default: 0)
#' @param multivariate Logical indicating multivariate data (default: FALSE)
#' @param ... Additional arguments passed to powerprior functions
#'
#' @return Bayes factor (BF_12 = P(data|a0_1) / P(data|a0_2))
#'
#' @details
#' The Bayes factor compares the marginal likelihoods under two different
#' discounting parameters. BF > 1 favors a0_1, BF < 1 favors a0_2.
#'
#' Note: This is a simple approximation using the observed data likelihood
#' evaluated at posterior means.
#'
#' @examples
#' set.seed(123)
#' historical <- rnorm(50, mean = 10, sd = 2)
#' current <- rnorm(30, mean = 10.5, sd = 2)
#'
#' # Compare moderate borrowing (0.5) vs no borrowing (0)
#' bf <- bayes_factor(historical, current, a0_1 = 0.5, a0_2 = 0)
#' cat("Bayes Factor:", bf, "\n")
#'
#' @export
bayes_factor <- function(historical_data, current_data,
                         a0_1, a0_2 = 0,
                         multivariate = FALSE,
                         ...) {

  warning("This is a simplified Bayes factor approximation. For rigorous comparison, use proper marginal likelihood computation.")

  if (multivariate) {
    # Model 1
    pp1 <- powerprior_multivariate(historical_data, a0 = a0_1, ...)
    post1 <- posterior_multivariate(pp1, current_data)

    # Model 2
    pp2 <- powerprior_multivariate(historical_data, a0 = a0_2, ...)
    post2 <- posterior_multivariate(pp2, current_data)

    # Approximate using BIC or similar criteria
    # This is a placeholder - proper implementation would compute marginal likelihood
    n <- nrow(current_data)
    p <- ncol(current_data)

    # Use effective sample size in penalty
    ess1 <- effective_sample_size(nrow(historical_data), a0_1)
    ess2 <- effective_sample_size(nrow(historical_data), a0_2)

    bf <- exp((ess1 - ess2) / 2)

  } else {
    # Model 1
    pp1 <- powerprior_univariate(historical_data, a0 = a0_1, ...)
    post1 <- posterior_univariate(pp1, current_data)

    # Model 2
    pp2 <- powerprior_univariate(historical_data, a0 = a0_2, ...)
    post2 <- posterior_univariate(pp2, current_data)

    n <- length(current_data)
    ess1 <- effective_sample_size(length(historical_data), a0_1)
    ess2 <- effective_sample_size(length(historical_data), a0_2)

    bf <- exp((ess1 - ess2) / 2)
  }

  return(bf)
}


#' Check Data Compatibility
#'
#' @description
#' Performs a simple compatibility check between historical and current data
#' to help guide the choice of discounting parameter.
#'
#' @param historical_data Historical data
#' @param current_data Current data
#' @param alpha Significance level for compatibility test (default: 0.05)
#'
#' @return A list with compatibility test results
#'
#' @examples
#' set.seed(123)
#' historical <- rnorm(50, mean = 10, sd = 2)
#' current <- rnorm(30, mean = 10.5, sd = 2)
#'
#' compatibility <- check_compatibility(historical, current)
#' print(compatibility)
#'
#' @export
check_compatibility <- function(historical_data, current_data, alpha = 0.05) {

  if (is.matrix(historical_data) || is.data.frame(historical_data)) {
    # Multivariate case - use Hotelling's T-squared
    warning("Multivariate compatibility check not fully implemented. Using univariate approach on first variable.")
    historical_data <- historical_data[, 1]
    current_data <- current_data[, 1]
  }

  # Two-sample t-test
  test_result <- t.test(historical_data, current_data)

  compatible <- test_result$p.value > alpha

  # Suggest a0 based on compatibility
  if (compatible) {
    suggested_a0 <- 0.7  # High borrowing if compatible
    recommendation <- "Data appear compatible. Consider moderate to high borrowing (a0 = 0.5-1.0)."
  } else {
    suggested_a0 <- 0.3  # Low borrowing if incompatible
    recommendation <- "Data show some differences. Consider low to moderate borrowing (a0 = 0.0-0.5)."
  }

  result <- list(
    test = "Two-sample t-test",
    statistic = test_result$statistic,
    p_value = test_result$p.value,
    compatible = compatible,
    alpha = alpha,
    historical_mean = mean(historical_data),
    current_mean = mean(current_data),
    difference = mean(current_data) - mean(historical_data),
    suggested_a0 = suggested_a0,
    recommendation = recommendation
  )

  class(result) <- c("compatibility_check", "list")
  return(result)
}


#' Print method for compatibility_check
#' @param x Object of class "compatibility_check"
#' @param digits Number of digits to round to
#' @param ... Additional arguments
#' @return Invisibly returns the input object (for method chaining)
#' @export
print.compatibility_check <- function(x, digits = 4, ...) {
  cat("Data Compatibility Check\n")
  cat("========================\n\n")
  cat("Test:", x$test, "\n")
  cat("Test statistic:", round(x$statistic, digits), "\n")
  cat("P-value:", round(x$p_value, digits), "\n")
  cat("Compatible at alpha =", x$alpha, ":", x$compatible, "\n\n")

  cat("Historical mean:", round(x$historical_mean, digits), "\n")
  cat("Current mean:", round(x$current_mean, digits), "\n")
  cat("Difference:", round(x$difference, digits), "\n\n")

  cat("Suggested a0:", x$suggested_a0, "\n")
  cat("Recommendation:", x$recommendation, "\n")

  invisible(x)
}
