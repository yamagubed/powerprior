# Unit Tests for powerprior Package
# Using testthat framework

library(testthat)
library(powerprior)
library(MASS)

# ==============================================================================
# Tests for Univariate Functions
# ==============================================================================

test_that("powerprior_univariate works with vague prior", {
  set.seed(123)
  historical <- rnorm(50, mean = 10, sd = 2)

  pp <- powerprior_univariate(historical, a0 = 0.5)

  expect_s3_class(pp, "powerprior_univariate")
  expect_true(pp$vague_prior)
  expect_equal(pp$a0, 0.5)
  expect_equal(pp$m, 50)
  expect_true(is.numeric(pp$mu_n))
  expect_true(is.numeric(pp$kappa_n))
  expect_true(is.numeric(pp$nu_n))
  expect_true(is.numeric(pp$sigma2_n))
})

test_that("powerprior_univariate works with informative prior", {
  set.seed(123)
  historical <- rnorm(50, mean = 10, sd = 2)

  pp <- powerprior_univariate(
    historical,
    a0 = 0.5,
    mu0 = 9.5,
    kappa0 = 2,
    nu0 = 5,
    sigma2_0 = 4
  )

  expect_s3_class(pp, "powerprior_univariate")
  expect_false(pp$vague_prior)
  expect_equal(pp$mu0, 9.5)
  expect_equal(pp$kappa0, 2)
})

test_that("powerprior_univariate validates inputs", {
  expect_error(powerprior_univariate(c(1), a0 = 0.5),
               "at least 2 observations")
  expect_error(powerprior_univariate(c(1, 2, 3), a0 = 1.5),
               "between 0 and 1")
  expect_error(powerprior_univariate(c(1, 2, 3), a0 = -0.1),
               "between 0 and 1")
})

test_that("posterior_univariate works correctly", {
  set.seed(123)
  historical <- rnorm(50, mean = 10, sd = 2)
  current <- rnorm(30, mean = 10.5, sd = 2)

  pp <- powerprior_univariate(historical, a0 = 0.5)
  posterior <- posterior_univariate(pp, current)

  expect_s3_class(posterior, "posterior_univariate")
  expect_equal(posterior$n, 30)
  expect_true(is.numeric(posterior$mu_star))
  expect_true(is.numeric(posterior$kappa_star))
  expect_true(is.numeric(posterior$nu_star))
  expect_true(is.numeric(posterior$sigma2_star))
})

test_that("sample_posterior_univariate generates valid samples", {
  set.seed(123)
  historical <- rnorm(50, mean = 10, sd = 2)
  current <- rnorm(30, mean = 10.5, sd = 2)

  pp <- powerprior_univariate(historical, a0 = 0.5)
  posterior <- posterior_univariate(pp, current)

  # Joint sampling
  samples_joint <- sample_posterior_univariate(posterior, n_samples = 100)
  expect_true(is.matrix(samples_joint))
  expect_equal(nrow(samples_joint), 100)
  expect_equal(ncol(samples_joint), 2)
  expect_true(all(samples_joint[, "sigma2"] > 0))

  # Marginal sampling
  samples_marginal <- sample_posterior_univariate(posterior, n_samples = 100,
                                                  marginal = TRUE)
  expect_true(is.numeric(samples_marginal))
  expect_equal(length(samples_marginal), 100)
})

test_that("extreme a0 values work correctly", {
  set.seed(123)
  historical <- rnorm(50, mean = 10, sd = 2)
  current <- rnorm(30, mean = 10.5, sd = 2)

  # a0 = 0: no borrowing
  pp0 <- powerprior_univariate(historical, a0 = 0)
  post0 <- posterior_univariate(pp0, current)
  expect_true(abs(post0$mu_star - mean(current)) < 0.5)

  # a0 = 1: full borrowing
  pp1 <- powerprior_univariate(historical, a0 = 1)
  post1 <- posterior_univariate(pp1, current)
  combined_mean <- mean(c(historical, current))
  expect_true(abs(post1$mu_star - combined_mean) < 0.5)
})


# ==============================================================================
# Tests for Multivariate Functions
# ==============================================================================

test_that("powerprior_multivariate works with vague prior", {
  set.seed(456)
  Sigma_true <- matrix(c(4, 1.5, 1.5, 3), nrow = 2)
  historical <- mvrnorm(60, mu = c(10, 5), Sigma = Sigma_true)

  pp <- powerprior_multivariate(historical, a0 = 0.6)

  expect_s3_class(pp, "powerprior_multivariate")
  expect_true(pp$vague_prior)
  expect_equal(pp$a0, 0.6)
  expect_equal(pp$m, 60)
  expect_equal(pp$p, 2)
  expect_true(is.numeric(pp$mu_n))
  expect_true(is.matrix(pp$Lambda_n))
})

test_that("powerprior_multivariate works with informative prior", {
  set.seed(456)
  Sigma_true <- matrix(c(4, 1.5, 1.5, 3), nrow = 2)
  historical <- mvrnorm(60, mu = c(10, 5), Sigma = Sigma_true)

  pp <- powerprior_multivariate(
    historical,
    a0 = 0.6,
    mu0 = c(9.5, 4.8),
    kappa0 = 2,
    nu0 = 6,
    Lambda0 = diag(2) * 2
  )

  expect_s3_class(pp, "powerprior_multivariate")
  expect_false(pp$vague_prior)
  expect_equal(pp$mu0, c(9.5, 4.8))
  expect_equal(pp$kappa0, 2)
})

test_that("powerprior_multivariate validates inputs", {
  expect_error(powerprior_multivariate(matrix(1:4, nrow = 1), a0 = 0.5),
               "at least 2 observations")
  expect_error(powerprior_multivariate(matrix(1:10, nrow = 5), a0 = 1.5),
               "between 0 and 1")
})

test_that("posterior_multivariate works correctly", {
  set.seed(456)
  Sigma_true <- matrix(c(4, 1.5, 1.5, 3), nrow = 2)
  historical <- mvrnorm(60, mu = c(10, 5), Sigma = Sigma_true)
  current <- mvrnorm(40, mu = c(10.3, 5.2), Sigma = Sigma_true)

  pp <- powerprior_multivariate(historical, a0 = 0.6)
  posterior <- posterior_multivariate(pp, current)

  expect_s3_class(posterior, "posterior_multivariate")
  expect_equal(posterior$n, 40)
  expect_equal(posterior$p, 2)
  expect_true(is.numeric(posterior$mu_star))
  expect_true(is.matrix(posterior$Lambda_star))
})

test_that("sample_posterior_multivariate generates valid samples", {
  set.seed(456)
  Sigma_true <- matrix(c(4, 1.5, 1.5, 3), nrow = 2)
  historical <- mvrnorm(60, mu = c(10, 5), Sigma = Sigma_true)
  current <- mvrnorm(40, mu = c(10.3, 5.2), Sigma = Sigma_true)

  pp <- powerprior_multivariate(historical, a0 = 0.6)
  posterior <- posterior_multivariate(pp, current)

  # Marginal sampling
  samples_marginal <- sample_posterior_multivariate(posterior, n_samples = 50,
                                                    marginal = TRUE)
  expect_true(is.matrix(samples_marginal))
  expect_equal(nrow(samples_marginal), 50)
  expect_equal(ncol(samples_marginal), 2)

  # Joint sampling
  samples_joint <- sample_posterior_multivariate(posterior, n_samples = 20,
                                                 marginal = FALSE)
  expect_true(is.list(samples_joint))
  expect_true("mu" %in% names(samples_joint))
  expect_true("Sigma" %in% names(samples_joint))
  expect_equal(nrow(samples_joint$mu), 20)
  expect_equal(dim(samples_joint$Sigma), c(2, 2, 20))
})

test_that("multivariate dimension consistency", {
  set.seed(456)
  Sigma_true <- matrix(c(4, 1.5, 1.5, 3), nrow = 2)
  historical <- mvrnorm(60, mu = c(10, 5), Sigma = Sigma_true)
  current <- mvrnorm(40, mu = c(10.3, 5.2, 0), Sigma = diag(3))

  pp <- powerprior_multivariate(historical, a0 = 0.6)

  expect_error(posterior_multivariate(pp, current),
               "same number of columns")
})


# ==============================================================================
# Integration Tests
# ==============================================================================

test_that("full workflow univariate produces reasonable results", {
  set.seed(789)
  true_mu <- 10
  true_sigma <- 2

  historical <- rnorm(100, mean = true_mu, sd = true_sigma)
  current <- rnorm(50, mean = true_mu, sd = true_sigma)

  pp <- powerprior_univariate(historical, a0 = 0.5)
  posterior <- posterior_univariate(pp, current)
  samples <- sample_posterior_univariate(posterior, n_samples = 1000,
                                         marginal = TRUE)

  # Check that posterior mean is close to true value
  expect_true(abs(mean(samples) - true_mu) < 1)

  # Check that 95% CI contains true value
  ci <- quantile(samples, c(0.025, 0.975))
  expect_true(true_mu >= ci[1] && true_mu <= ci[2])
})

test_that("full workflow multivariate produces reasonable results", {
  set.seed(789)
  true_mu <- c(10, 5)
  Sigma_true <- matrix(c(4, 1.5, 1.5, 3), nrow = 2)

  historical <- mvrnorm(100, mu = true_mu, Sigma = Sigma_true)
  current <- mvrnorm(50, mu = true_mu, Sigma = Sigma_true)

  pp <- powerprior_multivariate(historical, a0 = 0.5)
  posterior <- posterior_multivariate(pp, current)
  samples <- sample_posterior_multivariate(posterior, n_samples = 500,
                                           marginal = TRUE)

  # Check that posterior means are close to true values
  post_means <- colMeans(samples)
  expect_true(all(abs(post_means - true_mu) < 1))
})
