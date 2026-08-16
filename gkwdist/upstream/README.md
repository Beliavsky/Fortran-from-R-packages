
<!-- README.md is generated from README.Rmd. Please edit that file -->

# gkwdist: Generalized Kumaraswamy Distribution Family <img src="man/figures/gkwdist.png" align="right" height="138" alt="gkwdist logo" />

[![CRAN
status](https://www.r-pkg.org/badges/version/gkwdist)](https://CRAN.R-project.org/package=gkwdist)
[![R-CMD-check](https://github.com/evandeilton/gkwdist/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/evandeilton/gkwdist/actions/workflows/R-CMD-check.yaml)
[![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/gkwdist)](https://cran.r-project.org/package=gkwdist)
[![License:MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

`gkwdist` implements the **Generalized Kumaraswamy (GKw)** distribution
family and its seven nested sub-models for bounded continuous data on
$(0,1)$. All functions are implemented in **C++** via RcppArmadillo for
maximum computational efficiency.

**Key Features:** - Seven flexible distributions for proportions, rates,
and bounded data - Standard R distribution API: `d*`, `p*`, `q*`, `r*` -
Analytical log-likelihood, gradient, and Hessian functions - **All
functions implemented in C++** for optimal performance - 10-50× faster
than equivalent R implementations - Numerically stable for extreme
parameter values

------------------------------------------------------------------------

## Installation

``` r
# Install from CRAN
install.packages("gkwdist")

# Or Install from GitHub
# install.packages("devtools")
devtools::install_github("evandeilton/gkwdist")
```

------------------------------------------------------------------------

## The Distribution Family

### Complete Function Table

| Distribution | Code | Parameters | Functions |
|:---|:--:|:---|:---|
| **Generalized Kumaraswamy** | `gkw` | $\alpha, \beta, \gamma, \delta, \lambda$ | `dgkw`, `pgkw`, `qgkw`, `rgkw`, `llgkw`, `grgkw`, `hsgkw` |
| **Beta-Kumaraswamy** | `bkw` | $\alpha, \beta, \gamma, \delta$ | `dbkw`, `pbkw`, `qbkw`, `rbkw`, `llbkw`, `grbkw`, `hsbkw` |
| **Kumaraswamy-Kumaraswamy** | `kkw` | $\alpha, \beta, \delta, \lambda$ | `dkkw`, `pkkw`, `qkkw`, `rkkw`, `llkkw`, `grkkw`, `hskkw` |
| **Exponentiated Kumaraswamy** | `ekw` | $\alpha, \beta, \lambda$ | `dekw`, `pekw`, `qekw`, `rekw`, `llekw`, `grekw`, `hsekw` |
| **McDonald (Beta Power)** | `mc` | $\gamma, \delta, \lambda$ | `dmc`, `pmc`, `qmc`, `rmc`, `llmc`, `grmc`, `hsmc` |
| **Kumaraswamy** | `kw` | $\alpha, \beta$ | `dkw`, `pkw`, `qkw`, `rkw`, `llkw`, `grkw`, `hskw` |
| **Beta** | `beta_` | $\gamma, \delta$ | `dbeta_`, `pbeta_`, `qbeta_`, `rbeta_`, `llbeta`, `grbeta`, `hsbeta` |

### Function Types

**Distribution Functions (C++ implementation with R interface):**

- `d*()` — Probability density function (PDF)
- `p*()` — Cumulative distribution function (CDF)
- `q*()` — Quantile function (inverse CDF)
- `r*()` — Random number generation

**Analytical Functions for Maximum Likelihood (C++ implementation):**

All analytical functions use the signature: `function(par, data)` where
`par` is a numeric vector of parameters.

- `ll*(par, data)` — Negative log-likelihood:

$$-\ell(\boldsymbol{\theta}; \mathbf{x}) = -\sum_{i=1}^n \log f(x_i; \boldsymbol{\theta})$$

- `gr*(par, data)` — Negative gradient (negative score vector):

$$-\nabla_{\boldsymbol{\theta}} \ell(\boldsymbol{\theta}; \mathbf{x})$$

- `hs*(par, data)` — Negative Hessian matrix:

$$-\nabla^2_{\boldsymbol{\theta}} \ell(\boldsymbol{\theta}; \mathbf{x})$$

**Note:** These functions return **negative** values to facilitate
direct use with optimization routines like `optim()`, which perform
minimization by default.

------------------------------------------------------------------------

## Mathematical Specification

**Notation:** All parameters are strictly positive
($\alpha, \beta, \gamma, \delta, \lambda > 0$); support $x \in (0, 1)$.
The beta function is $B(a,b) = \Gamma(a)\Gamma(b)/\Gamma(a+b)$, and the
regularized incomplete beta function is $I_z(a,b) = B_z(a,b)/B(a,b)$
where $B_z(a,b) = \int_0^z t^{a-1}(1-t)^{b-1}dt$.

### 1. Generalized Kumaraswamy (GKw)

**Parameters:** $\alpha, \beta, \gamma, \delta, \lambda > 0$

**PDF:**

$$f_{\text{GKw}}(x; \alpha, \beta, \gamma, \delta, \lambda) = \frac{\lambda \alpha \beta}{B(\gamma, \delta+1)} x^{\alpha-1} (1-x^\alpha)^{\beta-1} \left[1-(1-x^\alpha)^\beta\right]^{\gamma\lambda-1} \left(1-\left[1-(1-x^\alpha)^\beta\right]^\lambda\right)^{\delta}$$

**CDF:**

$$F_{\text{GKw}}(x; \alpha, \beta, \gamma, \delta, \lambda) = I_{\left[1-(1-x^\alpha)^\beta\right]^\lambda}(\gamma, \delta+1)$$

**Quantile:** Numerical inversion of the CDF via root-finding
algorithms.

**Moments:** Analytical expressions not available in closed form.
Numerical integration or simulation methods required.

------------------------------------------------------------------------

### 2. Beta-Kumaraswamy (BKw)

**Relationship:** Special case of GKw with $\lambda = 1$

**PDF:**

$$f_{\text{BKw}}(x; \alpha, \beta, \gamma, \delta) = \frac{\alpha \beta}{B(\gamma, \delta+1)} x^{\alpha-1} (1-x^\alpha)^{\beta\delta+\beta-1} \left[1-(1-x^\alpha)^\beta\right]^{\gamma-1}$$

**CDF:**

$$F_{\text{BKw}}(x; \alpha, \beta, \gamma, \delta) = I_{1-(1-x^\alpha)^\beta}(\gamma, \delta+1)$$

**Quantile:** Numerical inversion via root-finding. For the inverse:

$$u = I_y(\gamma, \delta+1) \quad \text{where} \quad y = 1-(1-x^\alpha)^\beta$$

Solving for $x$:

$$x = \left(1-\left(1-I_u^{-1}(\gamma, \delta+1)\right)^{1/\beta}\right)^{1/\alpha}$$

------------------------------------------------------------------------

### 3. Kumaraswamy-Kumaraswamy (KKw)

**Relationship:** Special case of GKw with $\gamma = 1$

**PDF:**

$$f_{\text{KKw}}(x; \alpha, \beta, \delta, \lambda) = \lambda \alpha \beta (\delta+1) \, x^{\alpha-1} (1-x^\alpha)^{\beta-1} \left[1-(1-x^\alpha)^\beta\right]^{\lambda-1} \left(1-\left[1-(1-x^\alpha)^\beta\right]^\lambda\right)^{\delta}$$

**CDF:**

$$F_{\text{KKw}}(x; \alpha, \beta, \delta, \lambda) = 1 - \left(1-\left[1-(1-x^\alpha)^\beta\right]^\lambda\right)^{\delta+1}$$

**Quantile (closed-form):**

$$Q_{\text{KKw}}(p; \alpha, \beta, \delta, \lambda) = \left(1 - \left(1 - \left(1-(1-p)^{1/(\delta+1)}\right)^{1/\lambda}\right)^{1/\beta}\right)^{1/\alpha}$$

**Moments:** Analytical expressions not available in closed form.

------------------------------------------------------------------------

### 4. Exponentiated Kumaraswamy (EKw)

**Relationship:** Special case of GKw with $\gamma = 1$ and $\delta = 0$

**PDF:**

$$f_{\text{EKw}}(x; \alpha, \beta, \lambda) = \lambda \alpha \beta \, x^{\alpha-1} (1-x^\alpha)^{\beta-1} \left[1-(1-x^\alpha)^\beta\right]^{\lambda-1}$$

**CDF:**

$$F_{\text{EKw}}(x; \alpha, \beta, \lambda) = \left[1-(1-x^\alpha)^\beta\right]^\lambda$$

**Quantile (closed-form):**

$$Q_{\text{EKw}}(p; \alpha, \beta, \lambda) = \left(1-\left(1-p^{1/\lambda}\right)^{1/\beta}\right)^{1/\alpha}$$

**Moments:**

$$\mathbb{E}(X^r) = \lambda \sum_{k=0}^{\infty} \frac{(-1)^k \binom{\lambda-1}{k}}{k+1} \cdot \beta \, B\left(1 + \frac{r}{\alpha}, (k+1)\beta\right)$$

where the binomial coefficient is generalized:
$\binom{\lambda-1}{k} = \frac{(\lambda-1)(\lambda-2)\cdots(\lambda-k)}{k!}$.

------------------------------------------------------------------------

### 5. McDonald (Beta Power)

**Relationship:** Special case of GKw with $\alpha = \beta = 1$

**PDF:**

$$f_{\text{MC}}(x; \gamma, \delta, \lambda) = \frac{\lambda}{B(\gamma, \delta+1)} x^{\gamma\lambda-1} (1-x^\lambda)^{\delta}$$

**CDF:**

$$F_{\text{MC}}(x; \gamma, \delta, \lambda) = I_{x^\lambda}(\gamma, \delta+1)$$

**Quantile:**

$$Q_{\text{MC}}(p; \gamma, \delta, \lambda) = \left(I_p^{-1}(\gamma, \delta+1)\right)^{1/\lambda}$$

where $I_p^{-1}(\gamma, \delta+1)$ is the inverse regularized incomplete
beta function (quantile function of the Beta distribution with
parameters $\gamma$ and $\delta+1$).

**Moments:**

$$\mathbb{E}(X^r) = \frac{B(\gamma + r/\lambda, \delta+1)}{B(\gamma, \delta+1)}$$

which is valid for $r/\lambda > -\gamma$.

------------------------------------------------------------------------

### 6. Kumaraswamy (Kw)

**Relationship:** Special case of GKw with $\gamma = 1$, $\delta = 0$,
$\lambda = 1$

**PDF:**

$$f_{\text{Kw}}(x; \alpha, \beta) = \alpha \beta \, x^{\alpha-1} (1-x^\alpha)^{\beta-1}$$

**CDF:**

$$F_{\text{Kw}}(x; \alpha, \beta) = 1 - (1-x^\alpha)^\beta$$

**Quantile (closed-form):**

$$Q_{\text{Kw}}(p; \alpha, \beta) = \left(1-(1-p)^{1/\beta}\right)^{1/\alpha}$$

**Moments:**

$$\mathbb{E}(X^r) = \beta \, B\left(1 + \frac{r}{\alpha}, \beta\right) = \frac{\beta \, \Gamma(1+r/\alpha) \, \Gamma(\beta)}{\Gamma(1+r/\alpha+\beta)}$$

which is valid for $r/\alpha > -1$.

**Special Cases:**

$$\mathbb{E}(X) = \frac{\beta \, \Gamma(1+1/\alpha) \, \Gamma(\beta)}{\Gamma(1+1/\alpha+\beta)}$$

$$\mathbb{E}(X^2) = \frac{\beta \, \Gamma(1+2/\alpha) \, \Gamma(\beta)}{\Gamma(1+2/\alpha+\beta)}$$

$$\text{Var}(X) = \mathbb{E}(X^2) - \left(\mathbb{E}(X)\right)^2$$

------------------------------------------------------------------------

### 7. Beta

**Relationship:** Special case of GKw with
$\alpha = \beta = \lambda = 1$

**Note on Parameterization:** This package uses a parameterization where
the Beta distribution has parameters $(\gamma, \delta+1)$, corresponding
to the standard $\text{Beta}(\text{shape}_1, \text{shape}_2)$ with
$\text{shape}_1 = \gamma$ and $\text{shape}_2 = \delta + 1$.

**PDF:**

$$f_{\text{Beta}}(x; \gamma, \delta) = \frac{1}{B(\gamma, \delta+1)} x^{\gamma-1} (1-x)^{\delta}$$

**CDF:**

$$F_{\text{Beta}}(x; \gamma, \delta) = I_x(\gamma, \delta+1)$$

**Quantile:**

$$Q_{\text{Beta}}(p; \gamma, \delta) = I_p^{-1}(\gamma, \delta+1)$$

**Moments:**

$$\mathbb{E}(X^r) = \frac{B(\gamma+r, \delta+1)}{B(\gamma, \delta+1)} = \frac{\Gamma(\gamma+r) \, \Gamma(\gamma+\delta+1)}{\Gamma(\gamma) \, \Gamma(\gamma+\delta+1+r)}$$

$$\mathbb{E}(X) = \frac{\gamma}{\gamma+\delta+1}$$

$$\text{Var}(X) = \frac{\gamma(\delta+1)}{(\gamma+\delta+1)^2(\gamma+\delta+2)}$$

------------------------------------------------------------------------

## Hierarchical Structure

                                  GKw(α, β, γ, δ, λ)
                                  /               \
                               λ = 1             γ = 1
                                /                    \
                       BKw(α, β, γ, δ)         KKw(α, β, δ, λ)
                             |                          |
                         α = β = 1                    δ = 0
                             |                          |
                        MC(γ, δ, λ)              EKw(α, β, λ)
                             |                          |
                          λ = 1                       λ = 1
                             |                          |
                        Beta(γ, δ)                   Kw(α, β)

**Note:** The Beta distribution is obtained from MC by setting
$\lambda = 1$, or from GKw by setting $\alpha = \beta = \lambda = 1$.
The Kumaraswamy distribution is obtained from EKw by setting
$\lambda = 1$, or from GKw by setting $\gamma = 1$, $\delta = 0$,
$\lambda = 1$.

------------------------------------------------------------------------

## Usage Examples

### Example 1: Basic Distribution Functions

``` r
library(gkwdist)

# Set parameters
alpha <- 2; beta <- 3; gamma <- 1.5; delta <- 2; lambda <- 1.2
x <- seq(0.01, 0.99, length.out = 100)

# Density
dens <- dgkw(x, alpha, beta, gamma, delta, lambda)

# CDF
cdf <- pgkw(x, alpha, beta, gamma, delta, lambda)

# Quantiles
q <- qgkw(c(0.25, 0.5, 0.75), alpha, beta, gamma, delta, lambda)
print(q)

# Random generation
set.seed(123)
random_sample <- rgkw(1000, alpha, beta, gamma, delta, lambda)

# Visualization
par(mfrow = c(1, 2))
plot(x, dens, type = "l", lwd = 2, col = "blue",
     main = "GKw PDF", xlab = "x", ylab = "Density")
grid()

plot(x, cdf, type = "l", lwd = 2, col = "red",
     main = "GKw CDF", xlab = "x", ylab = "F(x)")
grid()
```

### Example 2: Comparing Distribution Families

``` r
library(gkwdist)

par(mfrow = c(1,1), mar = c(3,3,2,2))

x <- seq(0.001, 0.999, length.out = 500)

# Compute densities for all families
d_gkw  <- dgkw(x, 2, 3, 1.5, 2, 1.2)
d_bkw  <- dbkw(x, 2, 3, 1.5, 2)
d_kkw  <- dkkw(x, 2, 3, 2, 1.2)
d_ekw  <- dekw(x, 2, 3, 1.5)
d_mc   <- dmc(x, 2, 3, 1.2)
d_kw   <- dkw(x, 2, 5)
d_beta <- dbeta_(x, 2, 3)

# Plot comparison
plot(x, d_gkw, type = "l", lwd = 2, col = "black",
     ylim = c(0, max(d_gkw, d_bkw, d_kkw, d_ekw, d_mc, d_kw, d_beta)),
     main = "Distribution Family Comparison",
     xlab = "x", ylab = "Density")
lines(x, d_bkw, lwd = 2, col = "red")
lines(x, d_kkw, lwd = 2, col = "blue")
lines(x, d_ekw, lwd = 2, col = "green")
lines(x, d_mc,  lwd = 2, col = "purple")
lines(x, d_kw,  lwd = 2, col = "orange")
lines(x, d_beta, lwd = 2, col = "brown")
legend("topright",
       legend = c("GKw", "BKw", "KKw", "EKw", "MC", "Kw", "Beta"),
       col = c("black", "red", "blue", "green", "purple", "orange", "brown"),
       lwd = 2, cex = 0.85)
grid()
```

### Example 3: Maximum Likelihood Estimation Using optim()

``` r
library(gkwdist)

# Generate synthetic data from Kumaraswamy distribution
set.seed(2024)

n <- 2000
true_alpha <- 2.5
true_beta  <- 3.5
data <- rkw(n, true_alpha, true_beta)

# Get starting values
par_ini <- gkwgetstartvalues(data, family = "kw", n_starts = 2)

# Optimization using BFGS with analytical gradient
fit <- optim(
  par = par_ini,
  fn = llkw,  # Negative log-likelihood function
  gr = grkw,  # Negative gradient (negative score function)
  data = data,
  method = "BFGS",
  hessian = TRUE
)

# Standard errors from observed information matrix
se <- sqrt(diag(solve(fit$hessian)))

# 95% Confidence intervals
ci <- cbind(
  Lower    = fit$par - 1.96 * se,
  Estimate = fit$par,
  Upper    = fit$par + 1.96 * se
)
rownames(ci) <- c("alpha", "beta")

cat("True parameters:    ", true_alpha, true_beta, "\n")
cat("Estimated parameters:", fit$par, "\n\n")
print(ci)
```

### Example 4: Goodness-of-Fit Diagnostic Plot

``` r
library(gkwdist)

# Using fitted model from Example 3
x_grid <- seq(0.001, 0.999, length.out = 200)

# Fitted density
fitted_dens <- dkw(x_grid, fit$par[1], fit$par[2])

# True density (for comparison)
true_dens <- dkw(x_grid, true_alpha, true_beta)

# Diagnostic plot
hist(data, breaks = 30, probability = TRUE,
     col = "lightgray", border = "white",
     main = "Kumaraswamy Distribution Fit",
     xlab = "Data", ylab = "Density")
lines(x_grid, fitted_dens, col = "red", lwd = 2, lty = 1)
lines(x_grid, true_dens, col = "blue", lwd = 2, lty = 2)
legend("topright",
       legend = c("Observed Data", "Fitted Model", "True Model"),
       col = c("gray", "red", "blue"),
       lwd = c(10, 2, 2), lty = c(1, 1, 2))
grid()
```

### Example 5: Model Selection Using AIC and BIC

``` r
library(gkwdist)

# Generate data from Exponentiated Kumaraswamy
set.seed(456)
n <- 1500
data <- rekw(n, alpha = 2, beta = 3, lambda = 1.5)

# Define competing models
models <- list(
  Beta = list(
    nll = function(par) llbeta(par, data),
    gr = function(par) grbeta(par, data),
    start = gkwgetstartvalues(data, family = "beta"),
    npar = 2
  ),
  Kw = list(
    nll = function(par) llkw(par, data),
    gr = function(par) grkw(par, data),
    start = gkwgetstartvalues(data, family = "kw"),
    npar = 2
  ),
  EKw = list(
    nll = function(par) llekw(par, data),
    gr = function(par) grekw(par, data),
    start = gkwgetstartvalues(data, family = "ekw"),
    npar = 3
  ),
  MC = list(
    nll = function(par) llmc(par, data),
    gr = function(par) grmc(par, data),
    start = gkwgetstartvalues(data, family = "mc"),
    npar = 3
  )
)

# Fit all models using optim with analytical gradients
fits <- lapply(models, function(m) {
  optim(par = m$start, fn = m$nll, gr = m$gr, method = "BFGS")
})

# Extract log-likelihoods and compute information criteria
loglik <- sapply(fits, function(f) -f$value)
k <- sapply(models, `[[`, "npar")

results <- data.frame(
  Model  = names(models),
  Coefs  = sapply(fits, function(f) paste0(round(f$par, 3), collapse = "|")),
  LogLik = round(loglik, 2),
  nPar   = k,
  AIC    = round(-2 * loglik + 2 * k, 2),
  BIC    = round(-2 * loglik + k * log(n), 2)
)

# Sort by AIC
results <- results[order(results$AIC), ]
print(results, row.names = FALSE)
cat("\nBest model by AIC:", results$Model[1], "\n")
```

### Example 6: Using Analytical Functions Directly

``` r
library(gkwdist)

# Generate data from EKw
set.seed(789)
n <- 1000
data <- rekw(n, alpha = 2, beta = 3, lambda = 1.5)
params <- c(2, 3, 1.5)

# Compute analytical functions at true parameters
nll <- llekw(params, data)
neg_score <- grekw(params, data)
neg_hess <- hsekw(params, data)

# Fisher information matrix (negative of negative Hessian = Hessian)
fisher <- -neg_hess

# Asymptotic standard errors
se <- sqrt(diag(solve(fisher)))
names(se) <- c("alpha", "beta", "lambda")

cat("Negative log-likelihood:", nll, "\n")
cat("\nNegative score vector (should be close to zero at true params):\n")
print(neg_score)
cat("\nNegative Hessian matrix:\n")
print(neg_hess)
cat("\nFisher Information matrix:\n")
print(fisher)
cat("\nAsymptotic standard errors:\n")
print(se)
```

### Example 7: Q-Q Plot for Model Validation

``` r
library(gkwdist)

# Generate data and fit model
set.seed(101)
n <- 2000
data <- rkw(n, alpha = 2, beta = 3)

# Fit using optim
fit <- optim(
  par = c(1, 1),
  fn = function(par) llkw(par, data),
  gr = function(par) grkw(par, data),
  method = "BFGS"
)

# Theoretical quantiles
p <- ppoints(n)
theoretical_q <- qkw(p, fit$par[1], fit$par[2])

# Empirical quantiles
empirical_q <- sort(data)

# Q-Q plot
plot(theoretical_q, empirical_q,
     xlab = "Theoretical Quantiles",
     ylab = "Empirical Quantiles",
     main = "Q-Q Plot: Kumaraswamy Distribution",
     pch = 19, col = rgb(0, 0, 1, 0.5))
abline(0, 1, col = "red", lwd = 2, lty = 2)
grid()
```

### Example 8: Fitting GKw (Full Model) Using optim()

``` r
library(gkwdist)

# Generate data from GKw
set.seed(2025)
n <- 10000
true_params <- c(alpha = 2, beta = 2.5, gamma = 1.5, delta = 2, lambda = 1.3)
data <- rgkw(n, true_params[1], true_params[2], true_params[3],
             true_params[4], true_params[5])

# Get starting values
par_ini <- gkwgetstartvalues(data, family = "gkw", n_starts = 5)

# Fit using optim with analytical gradient
fit_gkw <- optim(
  par = par_ini, 
  fn = llgkw,    # Negative log-likelihood
  gr = grgkw,    # Negative gradient
  data = data,
  method = "BFGS",
  hessian = TRUE,
  control = list(maxit = 1000)
)

# Results
se <- sqrt(diag(solve(fit_gkw$hessian)))
estimates <- data.frame(
  Parameter = names(true_params),
  True = true_params,
  Estimate = fit_gkw$par,
  SE = se
)
print(estimates, digits = 4)
```

### Example 9: Comparing Analytical vs Numerical Derivatives

``` r
library(gkwdist)

# Generate data
set.seed(999)
n <- 10000
data <- rkw(n, alpha = 2, beta = 3)
par <- c(2, 3)

# Negative analytical gradient
neg_grad_analytical <- grkw(par, data)

# Numerical gradient (using finite differences)
neg_grad_numerical <- numDeriv::grad(
  func = function(p) llkw(p, data),
  x = par
)

# Compare
comparison <- data.frame(
  Parameter = c("alpha", "beta"),
  Analytical = neg_grad_analytical,
  Numerical = neg_grad_numerical,
  Difference = abs(neg_grad_analytical - neg_grad_numerical)
)
print(comparison, digits = 8)
```

------------------------------------------------------------------------

## Performance Comparison

The C++ implementation provides substantial performance gains:

``` r
library(microbenchmark)

# Generate large dataset
n <- 10000
data <- rkw(n, 2, 3)

# Compare log-likelihood computation
benchmark <- microbenchmark(
  R_sum_log_d = -sum(log(dkw(data, 2, 3))),  # Manual negative log-likelihood
  Cpp_ll = llkw(c(2, 3), data),               # C++ negative log-likelihood
  times = 100
)

print(benchmark)
plot(benchmark)
# Typical results: C++ implementation is 10-50× faster
```

**Why C++ Implementation Matters:**

- **Speed:** Critical for maximum likelihood estimation with large
  datasets
- **Numerical Stability:** Better handling of extreme parameter values
  and edge cases
- **Memory Efficiency:** Optimized memory allocation in tight loops
- **Scalability:** Linear scaling with sample size
- **Precision:** Analytical derivatives are exact (up to floating-point
  precision)

------------------------------------------------------------------------

## When to Use Each Distribution

| Data Characteristics | Recommended Distribution | Rationale |
|:---|:---|:---|
| Unimodal, symmetric | **Beta** | Parsimony; well-studied properties |
| Unimodal, asymmetric | **Kumaraswamy** | Closed-form CDF and quantile |
| Bimodal or U-shaped | **GKw** or **BKw** | Maximum flexibility |
| Extreme skewness | **KKw** or **EKw** | Fine control over tail behavior |
| J-shaped (monotonic) | **Kw** or **Beta** | With appropriate parameter values |
| Power transformations | **McDonald** | Explicit power parameter $\lambda$ |
| Unknown shape | **GKw** | Fit and test nested models |

**Model Selection Workflow:**

1.  Start with **Beta** and **Kumaraswamy** (2 parameters)
2.  Check goodness-of-fit using Q-Q plots and formal tests
3.  If inadequate, try **EKw** or **MC** (3 parameters)
4.  For complex patterns, use **BKw**, **KKw**, or **GKw** (4-5
    parameters)
5.  Use AIC/BIC to balance fit quality and parsimony
6.  Validate final model with residual diagnostics

------------------------------------------------------------------------

## References

- **Carrasco, J. M. F., Ferrari, S. L. P., & Cordeiro, G. M. (2010).** A
  new generalized Kumaraswamy distribution. *arXiv:1004.0911*.
  [arxiv.org/abs/1004.0911](https://arxiv.org/abs/1004.0911)

- **Cordeiro, G. M., & de Castro, M. (2011).** A new family of
  generalized distributions. *Journal of Statistical Computation and
  Simulation*, 81(7), 883-898.
  [doi:10.1080/00949650903530745](https://doi.org/10.1080/00949650903530745)

- **Kumaraswamy, P. (1980).** A generalized probability density function
  for double-bounded random processes. *Journal of Hydrology*, 46(1-2),
  79-88.
  [doi:10.1016/0022-1694(80)90036-0](https://doi.org/10.1016/0022-1694(80)90036-0)

- **Jones, M. C. (2009).** Kumaraswamy’s distribution: A beta-type
  distribution with some tractability advantages. *Statistical
  Methodology*, 6(1), 70-81.
  [doi:10.1016/j.stamet.2008.04.001](https://doi.org/10.1016/j.stamet.2008.04.001)

- **Nadarajah, S., Cordeiro, G. M., & Ortega, E. M. (2012).** The
  exponentiated Kumaraswamy distribution. *Journal of the Franklin
  Institute*, 349(3), 1180-1214.

- **McDonald, J. B. (1984).** Some generalized functions for the size
  distribution of income. *Econometrica*, 52(3), 647-663.
  [doi:10.2307/1913469](https://doi.org/10.2307/1913469)

- **Cordeiro, G. M., & Brito, R. S. (2012).** The beta power
  distribution. *Brazilian Journal of Probability and Statistics*,
  26(1), 88-112.
  [doi:10.1214/10-BJPS124](https://doi.org/10.1214/10-BJPS124)

------------------------------------------------------------------------

## Citation

``` r
citation("gkwdist")
```

``` bibtex
@Manual{gkwdist2025,
  title  = {gkwdist: Generalized Kumaraswamy Distribution Family},
  author = {José Evandeilton Lopes},
  year   = {2025},
  note   = {R package},
  url    = {https://github.com/evandeilton/gkwdist}
}
```

## Author

**José Evandeilton Lopes**  
LEG - Laboratory of Statistics and Geoinformation  
PPGMNE - Graduate Program in Numerical Methods in Engineering  
Federal University of Paraná (UFPR), Brazil  
Email: <evandeilton@gmail.com>

------------------------------------------------------------------------

## License

MIT License. See [LICENSE](LICENSE) file for details.
