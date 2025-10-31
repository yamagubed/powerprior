# powerprior: Conjugate Power Priors for Bayesian Analysis of Normal Data

## Overview

The `powerprior` package implements conjugate representations of power priors for efficient Bayesian analysis of normal data. This package provides closed-form computation of power priors and posterior distributions, eliminating the need for MCMC sampling.

**Key Features:**
- Univariate normal data using Normal-Inverse-Chi-squared (NIX) distributions
- Multivariate normal data using Normal-Inverse-Wishart (NIW) distributions
- Closed-form posterior updating (no MCMC required)
- Support for both informative and vague (non-informative) initial priors
- Efficient posterior sampling algorithms
- Ideal for simulation studies and clinical trial design

## Theoretical Background

Power priors provide a flexible framework for incorporating historical data with discounting:

$$P(\theta|x, a_0) = L(\theta|x)^{a_0}P(\theta)$$

where:
- $x$ is historical data
- $\theta$ are parameters of interest
- $a_0 \in [0,1]$ is the discounting parameter
- $a_0 = 0$: no borrowing from historical data
- $a_0 = 1$: full borrowing (treat as equivalent to current data)

### Dependencies

The package requires:
- `stats` (base R)
- `MASS`
- `LaplacesDemon`

Install dependencies:
```r
install.packages(c("MASS", "LaplacesDemon"))
```

## Quick Start

### Univariate Example

```r
library(powerprior)

# Generate historical and current data
set.seed(123)
historical <- rnorm(50, mean = 10, sd = 2)
current <- rnorm(30, mean = 10.5, sd = 2)

# Compute power prior with moderate discounting
pp <- powerprior_univariate(historical, a0 = 0.5)
print(pp)

# Update with current data to get posterior
posterior <- posterior_univariate(pp, current)
print(posterior)

# Sample from posterior
samples <- sample_posterior_univariate(posterior, n_samples = 1000)

# Summarize
cat("Posterior mean of mu:", mean(samples[, "mu"]), "\n")
cat("95% CI for mu:", quantile(samples[, "mu"], c(0.025, 0.975)), "\n")
```

### Multivariate Example

```r
library(powerprior)
library(MASS)

# Generate bivariate data
Sigma_true <- matrix(c(4, 1.5, 1.5, 3), nrow = 2)
historical <- mvrnorm(60, mu = c(10, 5), Sigma = Sigma_true)
current <- mvrnorm(40, mu = c(10.3, 5.2), Sigma = Sigma_true)

# Compute power prior
pp <- powerprior_multivariate(historical, a0 = 0.6)

# Update with current data
posterior <- posterior_multivariate(pp, current)

# Sample from marginal distribution of mu
samples <- sample_posterior_multivariate(posterior, n_samples = 1000, marginal = TRUE)

cat("Posterior mean:", colMeans(samples), "\n")
```

## Main Functions

### Univariate Analysis

| Function | Description |
|----------|-------------|
| `powerprior_univariate()` | Compute power prior for univariate normal data |
| `posterior_univariate()` | Update power prior with current data |
| `sample_posterior_univariate()` | Sample from posterior distribution |

### Multivariate Analysis

| Function | Description |
|----------|-------------|
| `powerprior_multivariate()` | Compute power prior for multivariate normal data |
| `posterior_multivariate()` | Update power prior with current data |
| `sample_posterior_multivariate()` | Sample from posterior distribution |

## Detailed Examples

### Example 1: Sensitivity Analysis

Explore how different discounting parameters affect posterior inference:

```r
a0_values <- seq(0, 1, by = 0.1)
results <- data.frame(a0 = a0_values, mean = NA, sd = NA)

for (i in seq_along(a0_values)) {
  pp <- powerprior_univariate(historical, a0 = a0_values[i])
  post <- posterior_univariate(pp, current)
  samples <- sample_posterior_univariate(post, n_samples = 2000, marginal = TRUE)
  
  results$mean[i] <- mean(samples)
  results$sd[i] <- sd(samples)
}

print(results)
```

### Example 2: Informative Prior

Use an informative NIX initial prior:

```r
pp_inform <- powerprior_univariate(
  historical_data = historical,
  a0 = 0.5,
  mu0 = 9.5,      # Prior mean
  kappa0 = 2,     # Prior precision
  nu0 = 5,        # Prior degrees of freedom
  sigma2_0 = 4    # Prior scale parameter
)

posterior_inform <- posterior_univariate(pp_inform, current)
```

### Example 3: Clinical Trial Simulation

Assess operating characteristics through simulation:

```r
n_simulations <- 1000
coverage_count <- 0
true_mean <- 10.5

for (sim in 1:n_simulations) {
  hist <- rnorm(50, mean = 10, sd = 2)
  curr <- rnorm(30, mean = true_mean, sd = 2)
  
  pp <- powerprior_univariate(hist, a0 = 0.5)
  post <- posterior_univariate(pp, curr)
  samples <- sample_posterior_univariate(post, n_samples = 1000, marginal = TRUE)
  
  ci <- quantile(samples, c(0.025, 0.975))
  if (true_mean >= ci[1] && true_mean <= ci[2]) {
    coverage_count <- coverage_count + 1
  }
}

cat("Coverage probability:", coverage_count / n_simulations, "\n")
```

## Computational Advantages

The conjugate representation provides significant computational benefits:

- **No MCMC required**: Closed-form posterior computation
- **Fast sampling**: Direct sampling from known distributions
- **Simulation studies**: Ideal for repeated posterior sampling
- **Clinical trial design**: Efficient operating characteristic evaluation

Traditional power prior analysis with MCMC might take hours for 1000 simulations; this package completes them in seconds.
