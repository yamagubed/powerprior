#' Compute Power Prior for Multivariate Normal Data
#'
#' @description
#' Computes the power prior for multivariate normal data using a conjugate
#' Normal-Inverse-Wishart (NIW) representation.
#'
#' @param historical_data Matrix or data frame where rows are observations and
#'   columns are variables. Must have at least 2 rows and 2 columns. Missing values
#'   are not automatically handled; rows with any NA values should be removed before
#'   calling this function.
#'
#' @param a0 Discounting parameter in `[0, 1]`. Controls the degree of borrowing from
#'   historical data in the multivariate setting. Same interpretation as univariate case:
#'   \itemize{
#'     \item a0 = 0: No borrowing from historical data
#'     \item a0 = 0.5: Moderate borrowing with half weight
#'     \item a0 = 1: Full borrowing of historical information
#'   }
#'
#' @param mu0 Prior mean vector of length p (number of variables).
#'   Only used when specifying an informative initial prior. If NULL, a vague
#'   (non-informative) prior is used. Represents the prior belief about the center
#'   of the multivariate distribution across all variables.
#'
#' @param kappa0 Prior precision parameter (scalar). Controls the concentration of
#'   the prior around mu0. Only used for informative priors. If NULL, vague prior
#'   is applied. Interpreted as the "prior effective sample size" for the mean
#'   estimate. Higher values indicate greater prior confidence.
#'
#' @param nu0 Prior degrees of freedom for the inverse Wishart distribution.
#'   Only used for informative priors. If NULL, vague prior is applied.
#'   Must satisfy nu0 >= p. Higher values indicate greater prior confidence
#'   in the covariance structure. Typical values range from p to 2p for informative priors.
#'
#' @param Lambda0 Prior scale matrix (p x p, symmetric positive definite).
#'   Only used for informative priors. If NULL, vague prior is applied.
#'   Represents the prior belief about the covariance structure. Larger values
#'   correspond to greater prior uncertainty about variable spread and relationships.
#'   If Lambda0 is provided, must be symmetric and positive definite.
#'
#' @return A list of class "powerprior_multivariate" containing:
#'   \item{mu_n}{Updated mean vector (p-dimensional)}
#'   \item{kappa_n}{Updated precision parameter (effective sample size)}
#'   \item{nu_n}{Updated degrees of freedom}
#'   \item{Lambda_n}{Updated scale matrix (p x p)}
#'   \item{a0}{Discounting parameter used}
#'   \item{m}{Sample size of historical data}
#'   \item{p}{Dimension (number of variables) of the data}
#'   \item{xbar}{Sample mean vector of historical data}
#'   \item{Sx}{Sum of squares and crossproducts matrix}
#'   \item{vague_prior}{Logical indicating if vague prior was used}
#'   \item{mu0}{Prior mean vector (if informative prior used)}
#'   \item{kappa0}{Prior precision parameter (if informative prior used)}
#'   \item{nu0}{Prior degrees of freedom (if informative prior used)}
#'   \item{Lambda0}{Prior scale matrix (if informative prior used)}
#'
#' @details
#' ## Background on Multivariate Power Priors
#'
#' The power prior framework extends naturally to multivariate normal data when
#' the mean vector \eqn{\mu} and covariance matrix \eqn{\Sigma} are jointly unknown. This is essential
#' for modern applications involving multiple correlated endpoints, such as clinical
#' trials measuring multiple health outcomes simultaneously.
#'
#' The power prior for multivariate data is defined analogously to the univariate case:
#'
#' \deqn{P(\boldsymbol{\mu}, \boldsymbol{\Sigma} | \mathbf{X}, a_0) =
#'   L(\boldsymbol{\mu}, \boldsymbol{\Sigma} | \mathbf{X})^{a_0} P(\boldsymbol{\mu}, \boldsymbol{\Sigma})}
#'
#' where:
#' - \eqn{\mathbf{X}} is the m x p historical data matrix
#' - \eqn{\boldsymbol{\mu}} is the p-dimensional mean vector
#' - \eqn{\boldsymbol{\Sigma}} is the p x p covariance matrix
#' - \eqn{a_0 \in [0, 1]} is the discounting parameter
#'
#' ## Conjugacy via Normal-Inverse-Wishart Distribution
#'
#' The key advantage of this implementation is that power priors applied to multivariate
#' normal data with Normal-Inverse-Wishart (NIW) conjugate initial priors remain in
#' closed form as NIW distributions. This preserves:
#'
#' \itemize{
#'   \item Exact posterior computation without approximation
#'   \item Closed-form parameter updates and marginalization
#'   \item Efficient sampling from standard distributions
#'   \item Computational tractability for high-dimensional problems (within reason)
#'   \item Natural joint modeling of correlations via the covariance structure
#' }
#'
#' For practitioners, this means you can incorporate historical information on
#' multiple correlated endpoints while maintaining full Bayesian rigor and computational efficiency.
#'
#' ## Informative vs. Vague Priors
#'
#' **Informative Initial Prior (all of mu0, kappa0, nu0, Lambda0 provided):**
#'
#' Uses a Normal-Inverse-Wishart (NIW) conjugate prior with hierarchical structure:
#'
#' \deqn{\boldsymbol{\mu} | \boldsymbol{\Sigma} \sim N_p(\boldsymbol{\mu}_0, \boldsymbol{\Sigma} / \kappa_0)}
#' \deqn{\boldsymbol{\Sigma} \sim \text{Inv-Wishart}(\nu_0, \boldsymbol{\Lambda}_0)}
#'
#' The power prior parameters are updated:
#'
#' \deqn{\boldsymbol{\mu}_n = \frac{a_0 m \overline{\mathbf{x}} + \kappa_0 \boldsymbol{\mu}_0}{a_0 m + \kappa_0}}
#'
#' \deqn{\kappa_n = a_0 m + \kappa_0}
#'
#' \deqn{\nu_n = a_0 m + \nu_0}
#'
#' \deqn{\boldsymbol{\Lambda}_n = a_0 \mathbf{S}_x + \boldsymbol{\Lambda}_0 +
#'   \frac{\kappa_0 a_0 m}{a_0 m + \kappa_0} (\overline{\mathbf{x}} - \boldsymbol{\mu}_0)
#'   (\overline{\mathbf{x}} - \boldsymbol{\mu}_0)^T}
#'
#' where \eqn{m} is sample size, \eqn{\overline{\mathbf{x}}} is the sample mean vector,
#' and \eqn{\mathbf{S}_x = \sum_{i=1}^m (\mathbf{x}_i - \overline{\mathbf{x}})
#' (\mathbf{x}_i - \overline{\mathbf{x}})^T} is the sum of squares and crossproducts matrix.
#'
#' **Vague (Non-informative) Initial Prior (all of mu0, kappa0, nu0, Lambda0 are NULL):**
#'
#' Uses a non-informative prior \eqn{P(\boldsymbol{\mu}, \boldsymbol{\Sigma}) \propto
#' |\boldsymbol{\Sigma}|^{-(p+1)/2}} that places minimal constraints on parameters.
#' The power prior parameters simplify to:
#'
#' \deqn{\boldsymbol{\mu}_n = \overline{\mathbf{x}}}
#'
#' \deqn{\kappa_n = a_0 m}
#'
#' \deqn{\nu_n = a_0 m - 1}
#'
#' \deqn{\boldsymbol{\Lambda}_n = a_0 \mathbf{S}_x}
#'
#' The vague prior is recommended when there is minimal prior information, or when
#' you want the analysis driven primarily by the historical data.
#'
#' ## Parameter Interpretation in Multivariate Setting
#'
#' **Effective Sample Size (\eqn{\kappa_n}):** Quantifies how much "effective historical data"
#' has been incorporated. The formula \eqn{\kappa_n = a_0 m + \kappa_0} shows the
#' discounted historical sample size combined with prior precision. This controls the
#' concentration of the posterior distribution for the mean vector.
#'
#' **Mean Vector (\eqn{\mu_n}):** The updated mean is a precision-weighted average:
#' \eqn{\boldsymbol{\mu}_n = \frac{a_0 m \overline{\mathbf{x}} + \kappa_0 \boldsymbol{\mu}_0}{a_0 m + \kappa_0}}
#' This naturally balances the historical sample mean and prior mean, with weights
#' proportional to their respective precisions.
#'
#' **Degrees of Freedom (\eqn{\nu_n}):** Controls tail behavior and the concentration of the
#' Wishart distribution. Higher values indicate greater confidence in the covariance
#' estimate. The minimum value needed is p (number of variables) for the Inverse-Wishart
#' to be well-defined; \eqn{\nu_n \geq p} is always maintained.
#'
#' **Scale Matrix (\eqn{\Lambda_n}):** The p x p scale matrix that captures both the dispersion
#' of individual variables and their correlations. The term
#' \eqn{\frac{\kappa_0 a_0 m}{a_0 m + \kappa_0} (\overline{\mathbf{x}} - \boldsymbol{\mu}_0)
#' (\overline{\mathbf{x}} - \boldsymbol{\mu}_0)^T} adds uncertainty when the historical
#' mean conflicts with the prior mean, naturally reflecting disagreement between data sources.
#'
#' ## Practical Considerations
#'
#' **Dimension:** This implementation works efficiently for moderate-dimensional problems
#' (typically p \eqn{\leq} 10). For higher dimensions, consider data reduction techniques or
#' structural assumptions on the covariance matrix.
#'
#' **Prior Specification:** When specifying Lambda0, ensure it is symmetric positive
#' definite. A simple approach is to use a multiple of the identity matrix (e.g.,
#' Lambda0 = diag(p)) for a weakly informative prior.
#'
#' **Discounting:** The same a0 parameter is used for all variables and their correlations.
#' If you suspect differential reliability of historical information across variables,
#' consider multiple analyses with different a0 values and sensitivity analyses.
#'
#' @references
#' Huang, Y., Yamaguchi, Y., Homma, G., Maruo, K., & Takeda, K. (2024).
#' "Conjugate Representation of Power Priors for Efficient Bayesian Analysis
#' of Normal Data." (unpublished).
#'
#' Ibrahim, J. G., & Chen, M. H. (2000). "Power prior distributions for regression models."
#' Statistical Science, 15(1), 46-60.
#'
#' Gelman, A., Carlin, J. B., Stern, H. S., et al. (2013).
#' Bayesian Data Analysis (3rd ed.). CRC Press.
#'
#' @examples
#' # Generate multivariate historical data with correlation
#' library(MASS)
#' Sigma_true <- matrix(c(4, 1, 1, 2), 2, 2)
#' historical <- mvrnorm(50, mu = c(10, 5), Sigma = Sigma_true)
#'
#' # Compute power prior with vague prior
#' pp <- powerprior_multivariate(historical, a0 = 0.5)
#' print(pp)
#'
#' # Compute power prior with informative prior
#' pp_inform <- powerprior_multivariate(
#'   historical,
#'   a0 = 0.5,
#'   mu0 = c(10, 5),
#'   kappa0 = 1,
#'   nu0 = 5,
#'   Lambda0 = diag(2)
#' )
#' print(pp_inform)
#'
#' @export
powerprior_multivariate <- function(historical_data, a0,
                                    mu0 = NULL, kappa0 = NULL,
                                    nu0 = NULL, Lambda0 = NULL) {

  # Convert to matrix if data frame
  if (is.data.frame(historical_data)) {
    historical_data <- as.matrix(historical_data)
  }

  # Validate inputs
  if (!is.matrix(historical_data) || !is.numeric(historical_data)) {
    stop("historical_data must be a numeric matrix or data frame")
  }

  if (nrow(historical_data) < 2) {
    stop("historical_data must have at least 2 observations")
  }

  if (ncol(historical_data) < 2) {
    stop("historical_data must have at least 2 variables (columns)")
  }

  if (!is.numeric(a0) || length(a0) != 1 || a0 < 0 || a0 > 1) {
    stop("a0 must be a single numeric value between 0 and 1")
  }

  # Get dimensions
  m <- nrow(historical_data)
  p <- ncol(historical_data)

  # Calculate sufficient statistics
  xbar <- colMeans(historical_data)
  X_centered <- sweep(historical_data, 2, xbar)
  Sx <- t(X_centered) %*% X_centered

  # Check if using informative or vague prior
  use_informative <- !is.null(mu0) && !is.null(kappa0) &&
    !is.null(nu0) && !is.null(Lambda0)

  if (use_informative) {
    # Validate informative prior parameters
    if (length(mu0) != p) {
      stop("mu0 must have length equal to number of columns in historical_data")
    }
    if (!is.matrix(Lambda0) || nrow(Lambda0) != p || ncol(Lambda0) != p) {
      stop("Lambda0 must be a p x p matrix")
    }
    if (!isSymmetric(Lambda0)) {
      stop("Lambda0 must be symmetric")
    }

    # Theorem 5: Informative NIW initial prior
    # Equations (17)-(20)
    mu_n <- (a0 * m * xbar + kappa0 * mu0) / (a0 * m + kappa0)
    kappa_n <- a0 * m + kappa0
    nu_n <- a0 * m + nu0

    diff_vec <- xbar - mu0
    Lambda_n <- a0 * Sx + Lambda0 +
      (kappa0 * a0 * m / (a0 * m + kappa0)) * (diff_vec %*% t(diff_vec))

  } else {
    # Corollary 1: Vague initial prior
    mu_n <- xbar
    kappa_n <- a0 * m
    nu_n <- a0 * m - 1
    Lambda_n <- a0 * Sx
  }

  result <- list(
    mu_n = mu_n,
    kappa_n = kappa_n,
    nu_n = nu_n,
    Lambda_n = Lambda_n,
    a0 = a0,
    m = m,
    p = p,
    xbar = xbar,
    Sx = Sx,
    vague_prior = !use_informative
  )

  if (use_informative) {
    result$mu0 <- mu0
    result$kappa0 <- kappa0
    result$nu0 <- nu0
    result$Lambda0 <- Lambda0
  }

  class(result) <- c("powerprior_multivariate", "list")
  return(result)
}


#' Posterior Obtained by Updating Power Prior with Current Data (Multivariate)
#'
#' @description
#' Updates a multivariate power prior with current trial data to obtain
#' the posterior distribution.
#'
#' @param powerprior Object of class "powerprior_multivariate" from
#'   powerprior_multivariate()
#' @param current_data Matrix or data frame of current observations. Must have the
#'   same number of columns as the historical data used to create the power prior.
#'   Rows represent observations, columns represent variables.
#'
#' @return A list of class "posterior_multivariate" containing:
#'   \item{mu_star}{Posterior mean vector}
#'   \item{kappa_star}{Posterior precision parameter}
#'   \item{nu_star}{Posterior degrees of freedom}
#'   \item{Lambda_star}{Posterior scale matrix}
#'   \item{n}{Sample size of current data}
#'   \item{p}{Dimension of data}
#'   \item{ybar}{Sample mean vector of current data}
#'   \item{Sy}{Sum of squares and crossproducts matrix for current data}
#'   \item{powerprior}{Original power prior object}
#'
#' @details
#' ## Posterior Updating in the Multivariate Setting
#'
#' Given a power prior P(\eqn{\mu}, \eqn{\Sigma} | X, \eqn{a_0}) and new current data Y, the posterior
#' distribution combines both through Bayes' theorem. The conjugate NIW structure
#' ensures the posterior remains a NIW distribution with closed-form parameters.
#'
#' The updating mechanism mirrors the univariate case but extends to handle the
#' covariance matrix structure. The scale matrix \eqn{\Lambda^*} incorporates:
#'
#' 1. Discounted historical sum of squares and crossproducts (\eqn{a_0} * S_x)
#' 2. Prior scale information (\eqn{\Lambda_0}, if using informative prior)
#' 3. Current data sum of squares and crossproducts (S_y)
#' 4. Disagreement terms that increase uncertainty when historical and current means differ
#'
#' The posterior mean vector is a precision-weighted average of the historical mean,
#' prior mean (if provided), and current mean, allowing for flexible incorporation
#' of multiple information sources with different precisions.
#'
#' @examples
#' library(MASS)
#' Sigma_true <- matrix(c(4, 1, 1, 2), 2, 2)
#' historical <- mvrnorm(50, mu = c(10, 5), Sigma = Sigma_true)
#' current <- mvrnorm(30, mu = c(10.5, 5.2), Sigma = Sigma_true)
#'
#' # With vague prior
#' pp <- powerprior_multivariate(historical, a0 = 0.5)
#' posterior <- posterior_multivariate(pp, current)
#' print(posterior)
#'
#' # With informative prior
#' pp_inform <- powerprior_multivariate(
#'   historical, a0 = 0.5,
#'   mu0 = c(10, 5), kappa0 = 1, nu0 = 5, Lambda0 = diag(2)
#' )
#' posterior_inform <- posterior_multivariate(pp_inform, current)
#' print(posterior_inform)
#'
#' @export
posterior_multivariate <- function(powerprior, current_data) {

  if (!inherits(powerprior, "powerprior_multivariate")) {
    stop("powerprior must be an object of class 'powerprior_multivariate'")
  }

  # Convert to matrix if data frame
  if (is.data.frame(current_data)) {
    current_data <- as.matrix(current_data)
  }

  if (!is.matrix(current_data) || !is.numeric(current_data)) {
    stop("current_data must be a numeric matrix or data frame")
  }

  if (ncol(current_data) != powerprior$p) {
    stop("current_data must have the same number of columns as historical data")
  }

  if (nrow(current_data) < 2) {
    stop("current_data must have at least 2 observations")
  }

  # Get dimensions
  n <- nrow(current_data)
  p <- powerprior$p

  # Calculate sufficient statistics for current data
  ybar <- colMeans(current_data)
  Y_centered <- sweep(current_data, 2, ybar)
  Sy <- t(Y_centered) %*% Y_centered

  # Extract power prior parameters
  mu_n <- powerprior$mu_n
  kappa_n <- powerprior$kappa_n
  nu_n <- powerprior$nu_n
  Lambda_n <- powerprior$Lambda_n
  a0 <- powerprior$a0
  m <- powerprior$m

  if (powerprior$vague_prior) {
    # Corollary 2: Posterior with vague initial prior
    # Equations (25)-(28)
    xbar <- powerprior$xbar
    Sx <- powerprior$Sx

    mu_star <- (a0 * m * xbar + n * ybar) / (a0 * m + n)
    kappa_star <- a0 * m + n
    nu_star <- a0 * m + n - 1

    diff_vec <- ybar - xbar
    Lambda_star <- a0 * Sx + Sy +
      (a0 * m * n / (a0 * m + n)) * (diff_vec %*% t(diff_vec))

  } else {
    # Theorem 6: Posterior with informative initial prior
    # Equations (21)-(24)
    xbar <- powerprior$xbar
    Sx <- powerprior$Sx
    mu0 <- powerprior$mu0
    kappa0 <- powerprior$kappa0
    nu0 <- powerprior$nu0
    Lambda0 <- powerprior$Lambda0

    mu_star <- (a0 * m * xbar + kappa0 * mu0 + n * ybar) / (a0 * m + kappa0 + n)
    kappa_star <- a0 * m + kappa0 + n
    nu_star <- a0 * m + nu0 + n

    diff_vec1 <- xbar - mu0
    diff_vec2 <- ybar - mu_n
    Lambda_star <- a0 * Sx + Lambda0 +
      (kappa0 * a0 * m / (a0 * m + kappa0)) * (diff_vec1 %*% t(diff_vec1)) +
      Sy +
      ((a0 * m + kappa0) * n / (a0 * m + kappa0 + n)) *
      (diff_vec2 %*% t(diff_vec2))
  }

  result <- list(
    mu_star = mu_star,
    kappa_star = kappa_star,
    nu_star = nu_star,
    Lambda_star = Lambda_star,
    n = n,
    p = p,
    ybar = ybar,
    Sy = Sy,
    powerprior = powerprior
  )

  class(result) <- c("posterior_multivariate", "list")
  return(result)
}


#' Sample from Posterior Distribution (Multivariate)
#'
#' @description
#' Generates samples from the multivariate posterior distribution using exact
#' closed-form expressions from the Normal-Inverse-Wishart conjugate family.
#'
#' @param posterior Object of class "posterior_multivariate" from
#'   posterior_multivariate()
#' @param n_samples Number of samples to generate (default: 1000).
#'   For large n_samples or high dimensions, computation time increases.
#' @param marginal Logical. If TRUE, samples only \eqn{\mu} from the multivariate
#'   t-distribution. If FALSE (default), samples the joint (\eqn{\mu}, \eqn{\Sigma}) from the
#'   NIW distribution, which is more computationally intensive but provides
#'   uncertainty in the covariance structure.
#'
#' @return If marginal=FALSE, a list with components:
#'   \item{mu}{n_samples x p matrix of mean samples}
#'   \item{Sigma}{p x p x n_samples array of covariance samples}
#'   If marginal=TRUE, an n_samples x p matrix of mean samples.
#'
#' @details
#' ## Sampling Algorithms
#'
#' **Joint Sampling (marginal=FALSE):**
#'
#' Implements the standard hierarchical sampling algorithm for the NIW distribution:
#'
#' \enumerate{
#'   \item Draw \eqn{\Sigma} ~ Inverse-Wishart(\eqn{\nu^*}, \eqn{\Lambda^*})
#'   \item Draw \eqn{\mu} | \eqn{\Sigma} ~ N_p(\eqn{\mu^*}, \eqn{\Sigma} / \eqn{\kappa^*})
#' }
#'
#' This produces samples from the joint distribution P(\eqn{\mu}, \eqn{\Sigma} | Y, X, \eqn{a_0}) and captures
#' both uncertainty in the mean and uncertainty in the covariance structure, including
#' their dependence.
#'
#' **Marginal Sampling (marginal=TRUE):**
#'
#' Uses the marginal t-distribution of the mean:
#'
#' \deqn{\mu | Y, X, a_0 \sim t_{\nu^*-p+1}(\mu^*, \Lambda^* / (\kappa^*(\nu^*-p+1)))}
#'
#' This is computationally faster and useful when you primarily care about inference
#' on the mean vector, marginalizing over uncertainty in the covariance.
#'
#' @examples
#' library(MASS)
#' Sigma_true <- matrix(c(4, 1, 1, 2), 2, 2)
#' historical <- mvrnorm(50, mu = c(10, 5), Sigma = Sigma_true)
#' current <- mvrnorm(30, mu = c(10.5, 5.2), Sigma = Sigma_true)
#'
#' pp <- powerprior_multivariate(historical, a0 = 0.5)
#' posterior <- posterior_multivariate(pp, current)
#'
#' # Sample from joint distribution (covariance included)
#' samples_joint <- sample_posterior_multivariate(posterior, n_samples = 100)
#' cat("Mean samples shape:", dim(samples_joint$mu), "\n")
#' cat("Covariance samples shape:", dim(samples_joint$Sigma), "\n")
#'
#' # Sample only means (faster)
#' samples_marginal <- sample_posterior_multivariate(posterior, n_samples = 100,
#'                                                    marginal = TRUE)
#'
#' @export
sample_posterior_multivariate <- function(posterior, n_samples = 1000,
                                          marginal = FALSE) {

  if (!inherits(posterior, "posterior_multivariate")) {
    stop("posterior must be an object of class 'posterior_multivariate'")
  }

  if (!is.numeric(n_samples) || n_samples < 1) {
    stop("n_samples must be a positive integer")
  }

  mu_star <- posterior$mu_star
  kappa_star <- posterior$kappa_star
  nu_star <- posterior$nu_star
  Lambda_star <- posterior$Lambda_star
  p <- posterior$p

  if (marginal) {
    # Sample from marginal multivariate t-distribution
    # mu | Y, X, a0 ~ t_{nu_star - p + 1}(mu_star, Lambda_star / (kappa_star * (nu_star - p + 1)))

    df <- nu_star - p + 1
    if (df <= 0) {
      stop("Degrees of freedom for marginal distribution is non-positive. Need more data.")
    }

    scale_matrix <- Lambda_star / (kappa_star * df)

    # Use LaplacesDemon::rmvt for multivariate t
    samples <- LaplacesDemon::rmvt(n_samples, mu = mu_star, S = scale_matrix, df = df)

    return(samples)

  } else {
    # Sample from joint NIW distribution

    # Step 1: Draw Sigma from Inverse-Wishart distribution
    Sigma_samples <- array(0, dim = c(p, p, n_samples))

    for (i in 1:n_samples) {
      # Draw from Inverse-Wishart(nu_star, Lambda_star)
      Sigma_samples[,,i] <- LaplacesDemon::rinvwishart(nu_star, Lambda_star)
    }

    # Step 2: Draw mu | Sigma from multivariate normal
    mu_samples <- matrix(0, nrow = n_samples, ncol = p)

    for (i in 1:n_samples) {
      Sigma_i <- Sigma_samples[,,i]
      cov_matrix <- Sigma_i / kappa_star
      mu_samples[i,] <- MASS::mvrnorm(1, mu = mu_star, Sigma = cov_matrix)
    }

    return(list(mu = mu_samples, Sigma = Sigma_samples))
  }
}


#' Print method for powerprior_multivariate
#' @param x Object of class "powerprior_multivariate"
#' @param ... Additional arguments
#' @export
print.powerprior_multivariate <- function(x, ...) {
  cat("Multivariate Power Prior (NIW Distribution)\n")
  cat("============================================\n\n")
  cat("Historical data:\n")
  cat("  Sample size (m):", x$m, "\n")
  cat("  Dimension (p):", x$p, "\n")
  cat("  Sample mean:", round(x$xbar, 4), "\n\n")
  cat("Discounting parameter (a0):", x$a0, "\n\n")
  cat("Prior type:", ifelse(x$vague_prior, "Vague (non-informative)", "Informative NIW"), "\n\n")
  cat("Power prior parameters:\n")
  cat("  mu_n:", round(x$mu_n, 4), "\n")
  cat("  kappa_n:", round(x$kappa_n, 4), "\n")
  cat("  nu_n:", round(x$nu_n, 4), "\n")
  cat("  Lambda_n:\n")
  print(round(x$Lambda_n, 4))
  invisible(x)
}


#' Print method for posterior_multivariate
#' @param x Object of class "posterior_multivariate"
#' @param ... Additional arguments
#' @export
print.posterior_multivariate <- function(x, ...) {
  cat("Multivariate Posterior Distribution (NIW)\n")
  cat("==========================================\n\n")
  cat("Current data:\n")
  cat("  Sample size (n):", x$n, "\n")
  cat("  Dimension (p):", x$p, "\n")
  cat("  Sample mean:", round(x$ybar, 4), "\n\n")
  cat("Posterior parameters:\n")
  cat("  mu*:", round(x$mu_star, 4), "\n")
  cat("  kappa*:", round(x$kappa_star, 4), "\n")
  cat("  nu*:", round(x$nu_star, 4), "\n")
  cat("  Lambda*:\n")
  print(round(x$Lambda_star, 4))
  cat("\nPosterior mean of mu:", round(x$mu_star, 4), "\n")
  invisible(x)
}
