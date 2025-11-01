#' Compute Power Prior for Univariate Normal Data
#'
#' @description
#' Computes the power prior for univariate normal data using a conjugate
#' Normal-Inverse-Chi-squared (NIX) representation.
#'
#' @param historical_data Numeric vector of historical observations. Must contain
#'   at least 2 observations. Missing values (NAs) are automatically removed.
#'   The function assumes data are independent and identically distributed from
#'   a normal distribution with unknown mean and variance.
#'
#' @param a0 Discounting parameter in `[0, 1]`. Controls the degree of borrowing from
#'   historical data. Specific values have intuitive interpretations:
#'   \itemize{
#'     \item a0 = 0: No borrowing from historical data; power prior equals the initial prior
#'     \item a0 = 0.5: Moderate borrowing; historical data contribute with half weight
#'     \item a0 = 1: Full borrowing; historical data weighted equally with current data
#'   }
#'   The parameter acts as a multiplicative discount on the historical likelihood.
#'
#' @param mu0 Prior mean for the normal distribution of the mean parameter.
#'   Only used when specifying an informative initial prior. If NULL, a vague
#'   (non-informative) prior is used instead. Represents the prior belief about
#'   the center of the data distribution.
#'
#' @param kappa0 Prior precision parameter (inverse of scaled variance) for the mean.
#'   Only used for informative priors. If NULL, vague prior is applied. Higher values
#'   indicate greater prior confidence in mu0. Interpreted as the "prior sample size"
#'   contributing to the mean estimate. For example, kappa0 = 1 is roughly equivalent
#'   to one prior observation.
#'
#' @param nu0 Prior degrees of freedom for the inverse chi-squared distribution of
#'   variance. Only used for informative priors. If NULL, vague prior is applied.
#'   Higher values indicate greater prior confidence in sigma2_0. Typically set to
#'   small positive values (e.g., 1-5) for weakly informative priors.
#'
#' @param sigma2_0 Prior scale parameter for the inverse chi-squared distribution.
#'   Only used for informative priors. If NULL, vague prior is applied.
#'   Represents the prior belief about the scale (spread) of the data. Should be
#'   on the variance scale (not standard deviation).
#'
#' @return A list of class "powerprior_univariate" containing:
#'   \item{mu_n}{Updated posterior mean parameter from power prior}
#'   \item{kappa_n}{Updated posterior precision parameter (effective sample size)}
#'   \item{nu_n}{Updated posterior degrees of freedom}
#'   \item{sigma2_n}{Updated posterior scale parameter (variance scale)}
#'   \item{a0}{Discounting parameter used}
#'   \item{m}{Sample size of historical data}
#'   \item{xbar}{Sample mean of historical data}
#'   \item{S}{Sum of squared deviations of historical data}
#'   \item{vague_prior}{Logical indicating if vague prior was used}
#'   \item{mu0}{Prior mean (if informative prior used)}
#'   \item{kappa0}{Prior precision (if informative prior used)}
#'   \item{nu0}{Prior degrees of freedom (if informative prior used)}
#'   \item{sigma2_0}{Prior scale parameter (if informative prior used)}
#'
#' @details
#' ## Background on Power Priors
#'
#' Power priors provide a principled framework for incorporating historical information
#' in Bayesian analysis while controlling the degree of borrowing through a discounting
#' parameter. The power prior is defined as:
#'
#' \deqn{P(\theta | x, a_0) = L(\theta | x)^{a_0} P(\theta)}
#'
#' where
#'
#' \itemize{
#'   \item \eqn{L(\theta | x)} is the likelihood function evaluated on historical data \eqn{x}
#'   \item \eqn{a_0 \in [0, 1]} is the discounting parameter
#'   \item \eqn{P(\theta)} is the initial prior distribution
#' }
#'
#' This approach is especially valuable in clinical trial design, where historical trial
#' data can improve statistical efficiency while maintaining appropriate skepticism about
#' whether historical and current populations are similar.
#'
#' ## Conjugacy and Computational Efficiency
#'
#' Typically, power priors result in non-closed-form posterior distributions requiring
#' computationally expensive Markov Chain Monte Carlo (MCMC) sampling, especially
#' problematic for simulation studies requiring thousands of posterior samples.
#'
#' This implementation exploits a key mathematical result: when applying power priors to
#' normal data with conjugate initial priors (Normal-Inverse-Chi-squared), the resulting
#' power prior and posterior distributions remain in closed form as NIX distributions.
#' This allows:
#'
#' \itemize{
#'   \item Direct computation without MCMC approximation
#'   \item Closed-form parameter updates
#'   \item Exact posterior sampling from standard distributions
#'   \item Efficient sensitivity analyses across different a0 values
#' }
#'
#' ## Informative vs. Vague Priors
#'
#' **Informative Initial Prior (all of mu0, kappa0, nu0, sigma2_0 provided):**
#'
#' Uses a Normal-Inverse-Chi-squared (NIX) conjugate prior with structure:
#'
#' \deqn{\mu | \sigma^2 \sim N(\mu_0, \sigma^2 / \kappa_0)}
#' \deqn{\sigma^2 \sim \text{Inv-}\chi^2(\nu_0, \sigma^2_0)}
#'
#' The power prior parameters are updated:
#'
#' \deqn{\mu_n = \frac{a_0 m \bar{x} + \kappa_0 \mu_0}{a_0 m + \kappa_0}}
#'
#' \deqn{\kappa_n = a_0 m + \kappa_0}
#'
#' \deqn{\nu_n = a_0 m + \nu_0}
#'
#' \deqn{\sigma^2_n = \frac{1}{\nu_n} \left[ a_0 S + \nu_0 \sigma^2_0 +
#'   \frac{a_0 m \kappa_0 (\bar{x} - \mu_0)^2}{a_0 m + \kappa_0} \right]}
#'
#' where \eqn{m} is the sample size, \eqn{\bar{x}} is the sample mean, and
#' \eqn{S = \sum_{i=1}^m (x_i - \bar{x})^2} is the sum of squared deviations.
#'
#' **Vague (Non-informative) Initial Prior (all of mu0, kappa0, nu0, sigma2_0 are NULL):**
#'
#' Uses a non-informative prior \eqn{P(\mu, \sigma^2) \propto \sigma^{-2}} that is
#' locally uniform in log-space. This places minimal prior constraints on the parameters.
#' The power prior parameters simplify to:
#'
#' \deqn{\mu_n = \bar{x}}
#'
#' \deqn{\kappa_n = a_0 m}
#'
#' \deqn{\nu_n = a_0 m - 1}
#'
#' \deqn{\sigma^2_n = \frac{S}{m}}
#'
#' The vague prior approach is recommended when there is no strong prior information,
#' or when you want the analysis to be primarily driven by the (discounted) historical data.
#'
#' ## Parameter Interpretation
#'
#' **Effective Sample Size (kappa_n):** The updated precision parameter can be interpreted
#' as an effective sample size. Higher values indicate more concentrated posterior distributions
#' for the mean. The formula \eqn{\kappa_n = a_0 m + \kappa_0} shows that the historical
#' sample size \eqn{m} is discounted by \eqn{a_0} before combining with the prior's
#' contribution \eqn{\kappa_0}.
#'
#' **Posterior Mean (mu_n):** A weighted average of the historical sample mean \eqn{\bar{x}},
#' prior mean \eqn{\mu_0}, and their relative precision:
#' \eqn{\mu_n = \frac{a_0 m \bar{x} + \kappa_0 \mu_0}{a_0 m + \kappa_0}}
#' This naturally interpolates between the data and prior, with weights determined by
#' their precision.
#'
#' **Degrees of Freedom (nu_n):** Controls the tail behavior of posterior distributions
#' derived from the NIX. Higher values produce lighter tails, indicating greater confidence.
#'
#' **Scale Parameter (sigma2_n):** Estimates the variability in the data. The term involving
#' \eqn{(\bar{x} - \mu_0)^2} captures disagreement between the historical mean and prior mean,
#' which increases uncertainty in variance estimation when they conflict.
#'
#' @references
#' Huang, Y., Yamaguchi, Y., Homma, G., Maruo, K., & Takeda, K. (2024).
#' "Conjugate Representation of Power Priors for Efficient Bayesian Analysis
#' of Normal Data." Statistical Science (unpublished).
#'
#' Ibrahim, J. G., & Chen, M. H. (2000). "Power prior distributions for regression models."
#' Statistical Science, 15(1), 46-60.
#'
#' Gelman, A., Carlin, J. B., Stern, H. S., et al. (2013).
#' Bayesian Data Analysis (3rd ed.). CRC Press.
#'
#' @examples
#' # Generate historical data
#' historical <- rnorm(50, mean = 10, sd = 2)
#'
#' # Compute power prior with informative initial prior
#' pp_inform <- powerprior_univariate(
#'   historical_data = historical,
#'   a0 = 0.5,
#'   mu0 = 10,
#'   kappa0 = 1,
#'   nu0 = 3,
#'   sigma2_0 = 4
#' )
#' print(pp_inform)
#'
#' # Compute power prior with vague prior (no prior specification)
#' pp_vague <- powerprior_univariate(
#'   historical_data = historical,
#'   a0 = 0.5
#' )
#' print(pp_vague)
#'
#' # Compare different discounting levels
#' pp_weak <- powerprior_univariate(historical_data = historical, a0 = 0.1)
#' pp_strong <- powerprior_univariate(historical_data = historical, a0 = 0.9)
#' cat("Strong discounting (a0=0.1) - kappa_n:", pp_weak$kappa_n, "\n")
#' cat("Weak discounting (a0=0.9) - kappa_n:", pp_strong$kappa_n, "\n")
#'
#' @export
powerprior_univariate <- function(historical_data, a0,
                                  mu0 = NULL, kappa0 = NULL,
                                  nu0 = NULL, sigma2_0 = NULL) {

  # Validate inputs
  if (!is.numeric(historical_data) || length(historical_data) < 2) {
    stop("historical_data must be a numeric vector with at least 2 observations")
  }

  if (!is.numeric(a0) || length(a0) != 1 || a0 < 0 || a0 > 1) {
    stop("a0 must be a single numeric value between 0 and 1")
  }

  # Remove NAs
  historical_data <- na.omit(historical_data)

  # Calculate sufficient statistics
  m <- length(historical_data)
  xbar <- mean(historical_data)
  S <- sum((historical_data - xbar)^2)

  # Check if using informative or vague prior
  use_informative <- !is.null(mu0) && !is.null(kappa0) &&
    !is.null(nu0) && !is.null(sigma2_0)

  if (use_informative) {
    # Theorem 1: Informative NIX initial prior
    # Equations (1)-(4)
    mu_n <- (a0 * m * xbar + kappa0 * mu0) / (a0 * m + kappa0)
    kappa_n <- a0 * m + kappa0
    nu_n <- a0 * m + nu0
    sigma2_n <- (1 / nu_n) * (a0 * S + nu0 * sigma2_0 +
                                (a0 * m * kappa0 * (xbar - mu0)^2) / (a0 * m + kappa0))

  } else {
    # Theorem 3: Vague initial prior
    # Equations (9)-(12)
    mu_n <- xbar
    kappa_n <- a0 * m
    nu_n <- a0 * m - 1
    sigma2_n <- S / m
  }

  result <- list(
    mu_n = mu_n,
    kappa_n = kappa_n,
    nu_n = nu_n,
    sigma2_n = sigma2_n,
    a0 = a0,
    m = m,
    xbar = xbar,
    S = S,
    vague_prior = !use_informative
  )

  if (use_informative) {
    result$mu0 <- mu0
    result$kappa0 <- kappa0
    result$nu0 <- nu0
    result$sigma2_0 <- sigma2_0
  }

  class(result) <- c("powerprior_univariate", "list")
  return(result)
}


#' Posterior Obtained by Updating Power Prior with Current Data (Univariate)
#'
#' @description
#' Updates a power prior with current trial data to obtain the posterior
#' distribution.
#'
#' @param powerprior Object of class "powerprior_univariate" from powerprior_univariate()
#' @param current_data Numeric vector of current trial observations. Must contain
#'   at least 2 observations. Missing values (NAs) are automatically removed.
#'
#' @return A list of class "posterior_univariate" containing:
#'   \item{mu_star}{Posterior mean parameter}
#'   \item{kappa_star}{Posterior precision parameter}
#'   \item{nu_star}{Posterior degrees of freedom}
#'   \item{sigma2_star}{Posterior scale parameter}
#'   \item{n}{Sample size of current data}
#'   \item{ybar}{Sample mean of current data}
#'   \item{Sy}{Sum of squared deviations of current data}
#'   \item{powerprior}{Original power prior object}
#'
#' @details
#' ## Posterior Updating
#'
#' Given a power prior distribution P(\eqn{\mu}, \eqn{\sigma^2} | x, \eqn{a_0}) and new current data y,
#' the posterior distribution is computed by combining the power prior with
#' the likelihood of the current data using Bayes' theorem.
#'
#' The conjugate structure ensures the posterior remains a NIX distribution.
#' For both informative and vague initial priors, the updating follows standard
#' conjugate rules, leveraging the fact that both the power prior and likelihood
#' are in the NIX family.
#'
#' This eliminates the computational burden of MCMC and allows direct posterior
#' inference and sampling (see \code{\link{sample_posterior_univariate}}).
#'
#' @examples
#' # Generate data
#' historical <- rnorm(50, mean = 10, sd = 2)
#' current <- rnorm(30, mean = 10.5, sd = 2)
#'
#' # Compute power prior and posterior
#' pp <- powerprior_univariate(historical, a0 = 0.5)
#' posterior <- posterior_univariate(pp, current)
#' print(posterior)
#'
#' # With informative prior
#' pp_inform <- powerprior_univariate(
#'   historical, a0 = 0.5,
#'   mu0 = 10, kappa0 = 1, nu0 = 3, sigma2_0 = 4
#' )
#' posterior_inform <- posterior_univariate(pp_inform, current)
#' print(posterior_inform)
#'
#' @export
posterior_univariate <- function(powerprior, current_data) {

  if (!inherits(powerprior, "powerprior_univariate")) {
    stop("powerprior must be an object of class 'powerprior_univariate'")
  }

  if (!is.numeric(current_data) || length(current_data) < 2) {
    stop("current_data must be a numeric vector with at least 2 observations")
  }

  # Remove NAs
  current_data <- na.omit(current_data)

  # Calculate sufficient statistics for current data
  n <- length(current_data)
  ybar <- mean(current_data)
  Sy <- sum((current_data - ybar)^2)

  # Extract power prior parameters
  mu_n <- powerprior$mu_n
  kappa_n <- powerprior$kappa_n
  nu_n <- powerprior$nu_n
  sigma2_n <- powerprior$sigma2_n
  a0 <- powerprior$a0
  m <- powerprior$m

  if (powerprior$vague_prior) {
    # Theorem 4: Posterior with vague initial prior
    # Equations (13)-(16)
    mu_star <- (a0 * m * powerprior$xbar + n * ybar) / (a0 * m + n)
    kappa_star <- a0 * m + n
    nu_star <- a0 * m + n - 1
    sigma2_star <- (1 / (a0 * m + n)) *
      (a0 * powerprior$S + Sy +
         (a0 * m * n * (powerprior$xbar - ybar)^2) / (a0 * m + n))

  } else {
    # Theorem 2: Posterior with informative initial prior
    # Equations (5)-(8)
    xbar <- powerprior$xbar
    S <- powerprior$S
    mu0 <- powerprior$mu0
    kappa0 <- powerprior$kappa0
    nu0 <- powerprior$nu0
    sigma2_0 <- powerprior$sigma2_0

    mu_star <- (a0 * m * xbar + kappa0 * mu0 + n * ybar) / (a0 * m + kappa0 + n)
    kappa_star <- a0 * m + kappa0 + n
    nu_star <- a0 * m + nu0 + n
    sigma2_star <- (1 / nu_star) *
      (a0 * S + nu0 * sigma2_0 +
         (a0 * m * kappa0 * (xbar - mu0)^2) / (a0 * m + kappa0) +
         Sy +
         (n * (a0 * m + kappa0) * (mu_n - ybar)^2) / (a0 * m + kappa0 + n))
  }

  result <- list(
    mu_star = mu_star,
    kappa_star = kappa_star,
    nu_star = nu_star,
    sigma2_star = sigma2_star,
    n = n,
    ybar = ybar,
    Sy = Sy,
    powerprior = powerprior
  )

  class(result) <- c("posterior_univariate", "list")
  return(result)
}


#' Sample from Posterior Distribution (Univariate)
#'
#' @description
#' Generates samples from the posterior distribution using the conjugate
#' representation. Can sample the joint distribution (mu, sigma2) or just
#' the marginal distribution of mu.
#'
#' @param posterior Object of class "posterior_univariate" from posterior_univariate()
#' @param n_samples Number of samples to generate (default: 1000)
#' @param marginal Logical. If TRUE, samples only mu from t-distribution.
#'   If FALSE, samples joint (mu, sigma2) from NIX distribution (default: FALSE)
#'
#' @return If marginal=FALSE, a matrix with columns "mu" and "sigma2".
#'   If marginal=TRUE, a vector of mu samples.
#'
#' @examples
#' # Generate data and compute posterior
#' historical <- rnorm(50, mean = 10, sd = 2)
#' current <- rnorm(30, mean = 10.5, sd = 2)
#' pp <- powerprior_univariate(historical, a0 = 0.5)
#' posterior <- posterior_univariate(pp, current)
#'
#' # Sample from joint distribution
#' samples_joint <- sample_posterior_univariate(posterior, n_samples = 1000)
#'
#' # Sample from marginal distribution of mu
#' samples_marginal <- sample_posterior_univariate(posterior, n_samples = 1000,
#'                                                  marginal = TRUE)
#'
#' @export
sample_posterior_univariate <- function(posterior, n_samples = 1000,
                                        marginal = FALSE) {

  if (!inherits(posterior, "posterior_univariate")) {
    stop("posterior must be an object of class 'posterior_univariate'")
  }

  if (!is.numeric(n_samples) || n_samples < 1) {
    stop("n_samples must be a positive integer")
  }

  mu_star <- posterior$mu_star
  kappa_star <- posterior$kappa_star
  nu_star <- posterior$nu_star
  sigma2_star <- posterior$sigma2_star

  if (marginal) {
    # Sample from marginal t-distribution
    # mu | y, x, a0 ~ t_{nu_star}(mu_star, sigma2_star / kappa_star)
    samples <- mu_star + sqrt(sigma2_star / kappa_star) *
      rt(n_samples, df = nu_star)

    return(samples)

  } else {
    # Sample from joint NIX distribution
    # Algorithm from Section 3.5

    # Step 1: Draw sigma2 from Inv-Chi-squared distribution
    # sigma2 ~ Inv-chi-squared(nu_star, sigma2_star)
    # This is equivalent to: sigma2 = nu_star * sigma2_star / chi-squared(nu_star)
    chi_sq_samples <- rchisq(n_samples, df = nu_star)
    sigma2_samples <- (nu_star * sigma2_star) / chi_sq_samples

    # Step 2: Draw mu | sigma2 from Normal distribution
    # mu | sigma2 ~ N(mu_star, sigma2 / kappa_star)
    mu_samples <- rnorm(n_samples, mean = mu_star,
                        sd = sqrt(sigma2_samples / kappa_star))

    samples <- cbind(mu = mu_samples, sigma2 = sigma2_samples)

    return(samples)
  }
}


#' Print method for powerprior_univariate
#' @param x Object of class "powerprior_univariate"
#' @param ... Additional arguments
#' @export
print.powerprior_univariate <- function(x, ...) {
  cat("Univariate Power Prior (NIX Distribution)\n")
  cat("==========================================\n\n")
  cat("Historical data:\n")
  cat("  Sample size (m):", x$m, "\n")
  cat("  Sample mean:", round(x$xbar, 4), "\n")
  cat("  Sum of squares:", round(x$S, 4), "\n\n")
  cat("Discounting parameter (a0):", x$a0, "\n\n")
  cat("Prior type:", ifelse(x$vague_prior, "Vague (non-informative)", "Informative NIX"), "\n\n")
  cat("Power prior parameters:\n")
  cat("  mu_n:", round(x$mu_n, 4), "\n")
  cat("  kappa_n:", round(x$kappa_n, 4), "\n")
  cat("  nu_n:", round(x$nu_n, 4), "\n")
  cat("  sigma2_n:", round(x$sigma2_n, 4), "\n")
  invisible(x)
}


#' Print method for posterior_univariate
#' @param x Object of class "posterior_univariate"
#' @param ... Additional arguments
#' @export
print.posterior_univariate <- function(x, ...) {
  cat("Univariate Posterior Distribution (NIX)\n")
  cat("=======================================\n\n")
  cat("Current data:\n")
  cat("  Sample size (n):", x$n, "\n")
  cat("  Sample mean:", round(x$ybar, 4), "\n\n")
  cat("Posterior parameters:\n")
  cat("  mu*:", round(x$mu_star, 4), "\n")
  cat("  kappa*:", round(x$kappa_star, 4), "\n")
  cat("  nu*:", round(x$nu_star, 4), "\n")
  cat("  sigma2*:", round(x$sigma2_star, 4), "\n\n")
  cat("Posterior mean of mu:", round(x$mu_star, 4), "\n")
  cat("Posterior variance of mu:", round(x$sigma2_star / x$kappa_star, 4), "\n")
  invisible(x)
}
