/**
 * @file bpmc.cpp
 * @brief Beta-Power / McDonald (BP/Mc) Distribution Functions
 * 
 * @details
 * This file implements the full suite of distribution functions for the
 * three-parameter Beta-Power (BP) or McDonald (Mc) distribution, which is
 * a sub-family of the Generalized Kumaraswamy (GKw) distribution obtained
 * by setting α = 1 and β = 1.
 * 
 * **Relationship to GKw:**
 * \deqn{BP(\gamma, \delta, \lambda) = GKw(1, 1, \gamma, \delta, \lambda)}
 * 
 * The BP distribution has probability density function:
 * \deqn{
 *   f(x; \gamma, \delta, \lambda) = 
 *   \frac{\lambda}{B(\gamma, \delta+1)} x^{\gamma\lambda-1} (1-x^\lambda)^\delta
 * }
 * for \eqn{x \in (0,1)}, where \eqn{B(\cdot,\cdot)} is the Beta function.
 * 
 * **Derivation from GKw:**
 * Setting α=1 and β=1 in the GKw PDF:
 * - \eqn{x^{\alpha-1} = x^0 = 1}
 * - \eqn{(1-x^\alpha)^{\beta-1} = (1-x)^0 = 1}
 * - \eqn{[1-(1-x^\alpha)^\beta] = [1-(1-x)] = x}
 * - \eqn{[1-(1-x^\alpha)^\beta]^{\gamma\lambda-1} = x^{\gamma\lambda-1}}
 * - \eqn{\{1-[1-(1-x^\alpha)^\beta]^\lambda\}^\delta = (1-x^\lambda)^\delta}
 * - The Beta function becomes: \eqn{B(\gamma, \delta+1)}
 * 
 * The cumulative distribution function is:
 * \deqn{
 *   F(x) = I_{x^\lambda}(\gamma, \delta+1)
 * }
 * where \eqn{I_y(a,b)} is the regularized incomplete Beta function.
 * 
 * The quantile function (inverse CDF) is:
 * \deqn{
 *   Q(p) = \left[Q_{Beta}(p; \gamma, \delta+1)\right]^{1/\lambda}
 * }
 * 
 * **Parameter Constraints:**
 * - \eqn{\gamma > 0} (shape parameter)
 * - \eqn{\delta \geq 0} (shape parameter)
 * - \eqn{\lambda > 0} (power parameter)
 * 
 * **Special Cases:**
 * | Distribution | Condition | Relation |
 * |--------------|-----------|----------|
 * | Power function | \eqn{\delta = 0} | BP(γ, 0, λ) |
 * | Beta | \eqn{\lambda = 1} | Beta(γ, δ+1) |
 * 
 * **Random Variate Generation:**
 * Uses transformation method:
 * 1. Generate \eqn{U \sim Beta(\gamma, \delta+1)}
 * 2. Return \eqn{X = U^{1/\lambda}}
 * 
 * **Numerical Stability:**
 * All computations use log-space arithmetic and numerically stable helper
 * functions from utils.h to prevent overflow/underflow.
 * 
 * **Implemented Functions:**
 * - dmc(): Probability density function (PDF)
 * - pmc(): Cumulative distribution function (CDF)
 * - qmc(): Quantile function (inverse CDF)
 * - rmc(): Random variate generation
 * - llmc(): Negative log-likelihood for MLE
 * - grmc(): Gradient of negative log-likelihood
 * - hsmc(): Hessian of negative log-likelihood
 * 
 * **Alternative Names:**
 * This distribution is also known as:
 * - McDonald distribution (McDonald, 1984)
 * - Generalized Beta of the first kind
 * - Libby-Novick Beta
 * 
 * @author Lopes, J. E.
 * @date 2025-01-07
 * 
 * @see gkw.cpp for the parent distribution
 * @see utils.h for numerical stability functions and parameter validators
 * 
 * @note All functions use R's vectorization conventions with parameter recycling.
 * @note Thread-safe: No global state is modified.
 */

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "utils.h"


// ============================================================================
// PROBABILITY DENSITY FUNCTION
// ============================================================================

/**
 * @brief Probability Density Function of the BP/McDonald Distribution
 * 
 * Computes the density (or log-density) for the Beta-Power distribution
 * at specified quantiles.
 * 
 * @param x Vector of quantiles (values in (0,1))
 * @param gamma Shape parameter vector (γ > 0)
 * @param delta Shape parameter vector (δ ≥ 0)
 * @param lambda Power parameter vector (λ > 0)
 * @param log_prob If TRUE, returns log-density; otherwise returns density
 * 
 * @return NumericVector of density values (or log-density if log_prob=TRUE)
 * 
 * @details
 * The log-density is computed as:
 * \deqn{
 *   \log f(x) = \log(\lambda) - \log B(\gamma, \delta+1)
 *   + (\gamma\lambda-1)\log(x) + \delta\log(1-x^\lambda)
 * }
 * 
 * @note Exported as .dmc_cpp for internal package use
 */
// [[Rcpp::export(.dmc_cpp)]]
Rcpp::NumericVector dmc(
    const arma::vec& x,
    const Rcpp::NumericVector& gamma,
    const Rcpp::NumericVector& delta,
    const Rcpp::NumericVector& lambda,
    bool log_prob = false
) {
  // Convert R vectors to Armadillo vectors
  arma::vec g_vec(gamma.begin(), gamma.size());
  arma::vec d_vec(delta.begin(), delta.size());
  arma::vec l_vec(lambda.begin(), lambda.size());
  
  // Determine output length for recycling
  size_t N = std::max({x.n_elem, g_vec.n_elem, d_vec.n_elem, l_vec.n_elem});
  
  // Initialize result with appropriate default
  arma::vec out(N);
  out.fill(log_prob ? R_NegInf : 0.0);
  
  for (size_t i = 0; i < N; i++) {
    // Extract recycled parameters
    double gg = g_vec[i % g_vec.n_elem];
    double dd = d_vec[i % d_vec.n_elem];
    double ll = l_vec[i % l_vec.n_elem];
    double xx = x[i % x.n_elem];
    
    // Validate parameters
    if (!check_bp_pars(gg, dd, ll)) {
      continue;
    }
    
    // Check support: x must be in (0, 1)
    if (xx <= 0.0 || xx >= 1.0 || !R_finite(xx)) {
      continue;
    }
    
    // ---- Log-space computation of density ----
    
    // Normalization constant: log(λ / B(γ, δ+1))
    double logB = R::lbeta(gg, dd + 1.0);
    double logCst = safe_log(ll) - logB;
    
    // Exponent: γλ - 1
    double exponent = gg * ll - 1.0;
    double lx = safe_log(xx);
    
    // Term 1: (γλ - 1) * log(x)
    double term1 = exponent * lx;
    
    // Compute x^λ
    double x_pow_l = safe_pow(xx, ll);
    if (x_pow_l >= 1.0) {
      continue;
    }
    
    // Term 2: δ * log(1 - x^λ)
    double log_1_minus_xpow = safe_log(1.0 - x_pow_l);
    double term2 = dd * log_1_minus_xpow;
    
    // Assemble log-density
    double log_pdf = logCst + term1 + term2;
    
    // Validate result
    if (!R_finite(log_pdf)) {
      continue;
    }
    
    // Return appropriate scale
    out(i) = log_prob ? log_pdf : safe_exp(log_pdf);
  }
  
  return Rcpp::NumericVector(out.memptr(), out.memptr() + out.n_elem);
}


// ============================================================================
// CUMULATIVE DISTRIBUTION FUNCTION
// ============================================================================

/**
 * @brief Cumulative Distribution Function of the BP/McDonald Distribution
 * 
 * Computes the cumulative probability for the Beta-Power distribution
 * at specified quantiles.
 * 
 * @param q Vector of quantiles
 * @param gamma Shape parameter vector (γ > 0)
 * @param delta Shape parameter vector (δ ≥ 0)
 * @param lambda Power parameter vector (λ > 0)
 * @param lower_tail If TRUE, returns P(X ≤ q); otherwise P(X > q)
 * @param log_p If TRUE, returns log-probability
 * 
 * @return NumericVector of cumulative probabilities
 * 
 * @details
 * The CDF is computed as:
 * \deqn{F(x) = I_{x^\lambda}(\gamma, \delta+1)}
 * where \eqn{I_y(a,b)} is the regularized incomplete Beta function.
 * 
 * @note Exported as .pmc_cpp for internal package use
 */
// [[Rcpp::export(.pmc_cpp)]]
Rcpp::NumericVector pmc(
    const arma::vec& q,
    const Rcpp::NumericVector& gamma,
    const Rcpp::NumericVector& delta,
    const Rcpp::NumericVector& lambda,
    bool lower_tail = true,
    bool log_p = false
) {
  // Convert R vectors to Armadillo vectors
  arma::vec g_vec(gamma.begin(), gamma.size());
  arma::vec d_vec(delta.begin(), delta.size());
  arma::vec l_vec(lambda.begin(), lambda.size());
  
  // Determine output length for recycling
  size_t N = std::max({q.n_elem, g_vec.n_elem, d_vec.n_elem, l_vec.n_elem});
  
  arma::vec out(N);
  
  for (size_t i = 0; i < N; i++) {
    // Extract recycled parameters
    double gg = g_vec[i % g_vec.n_elem];
    double dd = d_vec[i % d_vec.n_elem];
    double ll = l_vec[i % l_vec.n_elem];
    double xx = q[i % q.n_elem];
    
    // Validate parameters
    if (!check_bp_pars(gg, dd, ll)) {
      out(i) = NA_REAL;
      continue;
    }
    
    // Handle boundary: q ≤ 0
    if (!R_finite(xx) || xx <= 0.0) {
      double val0 = lower_tail ? 0.0 : 1.0;
      out(i) = log_p ? safe_log(val0) : val0;
      continue;
    }
    
    // Handle boundary: q ≥ 1
    if (xx >= 1.0) {
      double val1 = lower_tail ? 1.0 : 0.0;
      out(i) = log_p ? safe_log(val1) : val1;
      continue;
    }
    
    // ---- Compute CDF ----
    
    // Step 1: x^λ
    double xpow = safe_pow(xx, ll);
    
    // Step 2: F(x) = I_{x^λ}(γ, δ+1) via pbeta
    double val = R::pbeta(xpow, gg, dd + 1.0, true, false);
    
    // Apply tail adjustment
    if (!lower_tail) {
      val = 1.0 - val;
    }
    
    // Apply log transformation
    if (log_p) {
      val = safe_log(val);
    }
    
    out(i) = val;
  }
  
  return Rcpp::NumericVector(out.memptr(), out.memptr() + out.n_elem);
}


// ============================================================================
// QUANTILE FUNCTION
// ============================================================================

/**
 * @brief Quantile Function (Inverse CDF) of the BP/McDonald Distribution
 * 
 * Computes quantiles for the Beta-Power distribution
 * given probability values.
 * 
 * @param p Vector of probabilities (values in [0,1])
 * @param gamma Shape parameter vector (γ > 0)
 * @param delta Shape parameter vector (δ ≥ 0)
 * @param lambda Power parameter vector (λ > 0)
 * @param lower_tail If TRUE, probabilities are P(X ≤ x); otherwise P(X > x)
 * @param log_p If TRUE, probabilities are given as log(p)
 * 
 * @return NumericVector of quantiles
 * 
 * @details
 * The quantile function inverts the CDF:
 * \deqn{Q(p) = \left[Q_{Beta}(p; \gamma, \delta+1)\right]^{1/\lambda}}
 * 
 * @note Exported as .qmc_cpp for internal package use
 */
// [[Rcpp::export(.qmc_cpp)]]
Rcpp::NumericVector qmc(
    const arma::vec& p,
    const Rcpp::NumericVector& gamma,
    const Rcpp::NumericVector& delta,
    const Rcpp::NumericVector& lambda,
    bool lower_tail = true,
    bool log_p = false
) {
  // Convert R vectors to Armadillo vectors
  arma::vec g_vec(gamma.begin(), gamma.size());
  arma::vec d_vec(delta.begin(), delta.size());
  arma::vec l_vec(lambda.begin(), lambda.size());
  
  // Determine output length for recycling
  size_t N = std::max({p.n_elem, g_vec.n_elem, d_vec.n_elem, l_vec.n_elem});
  
  arma::vec out(N);
  
  for (size_t i = 0; i < N; i++) {
    // Extract recycled parameters
    double gg = g_vec[i % g_vec.n_elem];
    double dd = d_vec[i % d_vec.n_elem];
    double ll = l_vec[i % l_vec.n_elem];
    double pp = p[i % p.n_elem];
    
    // Validate parameters
    if (!check_bp_pars(gg, dd, ll)) {
      out(i) = NA_REAL;
      continue;
    }
    
    // ---- Convert probability to linear scale ----
    if (log_p) {
      if (pp > 0.0) {
        out(i) = NA_REAL;
        continue;
      }
      pp = safe_exp(pp);
    }
    
    // Handle upper tail (pp is now always linear scale)
    if (!lower_tail) {
      pp = 1.0 - pp;
    }
    
    // Handle boundary cases
    if (pp <= 0.0) {
      out(i) = 0.0;
      continue;
    }
    if (pp >= 1.0) {
      out(i) = 1.0;
      continue;
    }
    
    // ---- Compute quantile via inverse transformations ----
    
    // Step 1: y = Q_Beta(p, γ, δ+1)
    double y = R::qbeta(pp, gg, dd + 1.0, true, false);
    
    // Step 2: x = y^(1/λ)
    double xval;
    if (ll == 1.0) {
      xval = y;
    } else {
      xval = safe_pow(y, 1.0 / ll);
    }
    
    // Clamp to valid support
    xval = std::max(0.0, std::min(1.0, xval));
    
    out(i) = xval;
  }
  
  return Rcpp::NumericVector(out.memptr(), out.memptr() + out.n_elem);
}


// ============================================================================
// RANDOM NUMBER GENERATION
// ============================================================================

/**
 * @brief Random Variate Generation for the BP/McDonald Distribution
 * 
 * Generates random samples from the Beta-Power distribution
 * using the transformation method.
 * 
 * @param n Number of random variates to generate
 * @param gamma Shape parameter vector (γ > 0)
 * @param delta Shape parameter vector (δ ≥ 0)
 * @param lambda Power parameter vector (λ > 0)
 * 
 * @return NumericVector of n random variates from BP distribution
 * 
 * @details
 * Algorithm:
 * 1. Generate U ~ Beta(γ, δ+1)
 * 2. Return X = U^(1/λ)
 * 
 * @note Exported as .rmc_cpp for internal package use
 */
// [[Rcpp::export(.rmc_cpp)]]
Rcpp::NumericVector rmc(
    int n,
    const Rcpp::NumericVector& gamma,
    const Rcpp::NumericVector& delta,
    const Rcpp::NumericVector& lambda
) {
  if (n <= 0) {
    Rcpp::stop("rmc: n must be positive");
  }
  
  // Convert R vectors to Armadillo vectors
  arma::vec g_vec(gamma.begin(), gamma.size());
  arma::vec d_vec(delta.begin(), delta.size());
  arma::vec l_vec(lambda.begin(), lambda.size());
  
  arma::vec out(n);
  
  for (int i = 0; i < n; i++) {
    // Extract recycled parameters (direct modulo, no intermediate variable)
    double gg = g_vec[i % g_vec.n_elem];
    double dd = d_vec[i % d_vec.n_elem];
    double ll = l_vec[i % l_vec.n_elem];
    
    // Validate parameters
    if (!check_bp_pars(gg, dd, ll)) {
      out(i) = NA_REAL;
      Rcpp::warning("rmc: invalid parameters at index %d", i + 1);
      continue;
    }
    
    // Generate U ~ Beta(γ, δ+1)
    double U = R::rbeta(gg, dd + 1.0);
    
    // Transform: X = U^(1/λ)
    double xval;
    if (ll == 1.0) {
      xval = U;
    } else {
      xval = safe_pow(U, 1.0 / ll);
    }
    
    // Clamp to valid support
    xval = std::max(0.0, std::min(1.0, xval));
    
    out(i) = xval;
  }
  
  return Rcpp::NumericVector(out.memptr(), out.memptr() + out.n_elem);
}


// ============================================================================
// NEGATIVE LOG-LIKELIHOOD FUNCTION
// ============================================================================

/**
 * @brief Negative Log-Likelihood for BP/McDonald Distribution
 * 
 * Computes the negative log-likelihood function for parameter estimation
 * via maximum likelihood.
 * 
 * @param par Parameter vector of length 3: (γ, δ, λ)
 * @param data Vector of observations (must be in (0,1))
 * 
 * @return Negative log-likelihood value (scalar)
 * 
 * @details
 * The log-likelihood for n observations is:
 * \deqn{
 *   \ell(\theta) = n[\ln\lambda - \ln B(\gamma,\delta+1)]
 *   + (\gamma\lambda-1)\sum\ln x_i + \delta\sum\ln(1-x_i^\lambda)
 * }
 * 
 * Returns +Inf for invalid parameters or data outside (0,1).
 * 
 * @note Exported as .llmc_cpp for internal package use
 */
// [[Rcpp::export(.llmc_cpp)]]
double llmc(const Rcpp::NumericVector& par, const Rcpp::NumericVector& data) {
  // Validate parameter vector length
  if (par.size() < 3) return R_PosInf;
  
  // Extract parameters
  double gamma = par[0];
  double delta = par[1];
  double lambda = par[2];
  
  // Validate parameters using consistent checker
  if (!check_bp_pars(gamma, delta, lambda)) return R_PosInf;
  
  // Convert and validate data
  arma::vec x = Rcpp::as<arma::vec>(data);
  if (x.n_elem < 1) return R_PosInf;
  if (arma::any(x <= 0.0) || arma::any(x >= 1.0)) return R_PosInf;
  
  int n = x.n_elem;
  
  // Numerical stability constant
  const double eps = 1e-10;
  
  // Compute log(B(γ, δ+1)) stably
  double log_B;
  if (gamma > 100.0 || delta > 100.0) {
    log_B = lgamma(gamma) + lgamma(delta + 1.0) - lgamma(gamma + delta + 1.0);
  } else {
    log_B = R::lbeta(gamma, delta + 1.0);
  }
  
  // Constant term: n * [log(λ) - log(B(γ, δ+1))]
  double log_lambda = safe_log(lambda);
  double const_term = n * (log_lambda - log_B);
  
  // Calculate γλ - 1
  double gl_minus_1 = gamma * lambda - 1.0;
  
  // Initialize accumulators
  double sum_term1 = 0.0;  // (γλ-1) * Σlog(x)
  double sum_term2 = 0.0;  // δ * Σlog(1-x^λ)
  
  for (int i = 0; i < n; i++) {
    double xi = x(i);
    
    // Handle observations near boundaries
    xi = std::max(eps, std::min(1.0 - eps, xi));
    
    double log_xi = std::log(xi);
    
    // Term 1: (γλ-1) * log(x)
    sum_term1 += gl_minus_1 * log_xi;
    
    // Calculate x^λ stably
    double x_lambda;
    if (lambda * std::abs(log_xi) > 1.0) {
      x_lambda = safe_exp(lambda * log_xi);
    } else {
      x_lambda = std::pow(xi, lambda);
    }
    
    // Term 2: δ * log(1-x^λ)
    double log_1_minus_x_lambda;
    if (x_lambda > 0.9995) {
      log_1_minus_x_lambda = std::log1p(-x_lambda);
    } else {
      log_1_minus_x_lambda = safe_log(1.0 - x_lambda);
    }
    
    // Scale for large δ
    if (delta > 1000.0 && log_1_minus_x_lambda < -0.01) {
      double scaled_term = std::max(log_1_minus_x_lambda, -700.0 / delta);
      sum_term2 += delta * scaled_term;
    } else {
      sum_term2 += delta * log_1_minus_x_lambda;
    }
  }
  
  double loglike = const_term + sum_term1 + sum_term2;
  
  // Check for invalid results
  if (!std::isfinite(loglike)) return R_PosInf;
  
  return -loglike;
}


// ============================================================================
// GRADIENT OF NEGATIVE LOG-LIKELIHOOD
// ============================================================================

/**
 * @brief Gradient of Negative Log-Likelihood for BP/McDonald Distribution
 * 
 * Computes the gradient vector of the negative log-likelihood for
 * optimization-based parameter estimation.
 * 
 * @param par Parameter vector of length 3: (γ, δ, λ)
 * @param data Vector of observations (must be in (0,1))
 * 
 * @return NumericVector of length 3 containing partial derivatives
 *         with respect to (γ, δ, λ)
 * 
 * @details
 * The gradient components are:
 * - ∂ℓ/∂γ = -n[ψ(γ) - ψ(γ+δ+1)] + λ Σlog(x)
 * - ∂ℓ/∂δ = -n[ψ(δ+1) - ψ(γ+δ+1)] + Σlog(1-x^λ)
 * - ∂ℓ/∂λ = n/λ + γ Σlog(x) - δ Σ[x^λ log(x)/(1-x^λ)]
 * 
 * @note Exported as .grmc_cpp for internal package use
 */
// [[Rcpp::export(.grmc_cpp)]]
Rcpp::NumericVector grmc(const Rcpp::NumericVector& par, const Rcpp::NumericVector& data) {
  // Validate parameter vector length
  if (par.size() < 3) {
    return Rcpp::NumericVector(3, R_NaN);
  }
  
  // Extract parameters
  double gamma = par[0];
  double delta = par[1];
  double lambda = par[2];
  
  // Validate parameters using consistent checker
  if (!check_bp_pars(gamma, delta, lambda)) {
    return Rcpp::NumericVector(3, R_NaN);
  }
  
  // Convert and validate data
  arma::vec x = Rcpp::as<arma::vec>(data);
  if (x.n_elem < 1 || arma::any(x <= 0) || arma::any(x >= 1)) {
    return Rcpp::NumericVector(3, R_NaN);
  }
  
  int n = x.n_elem;
  Rcpp::NumericVector grad(3, 0.0);
  
  // Numerical stability constant
  const double eps = 1e-10;
  
  // Calculate digamma terms stably
  double digamma_gamma, digamma_delta_plus_1, digamma_gamma_delta_plus_1;
  
  if (gamma > 100.0) {
    digamma_gamma = std::log(gamma) - 1.0 / (2.0 * gamma);
  } else {
    digamma_gamma = R::digamma(gamma);
  }
  
  if (delta > 100.0) {
    digamma_delta_plus_1 = std::log(delta + 1.0) - 1.0 / (2.0 * (delta + 1.0));
  } else {
    digamma_delta_plus_1 = R::digamma(delta + 1.0);
  }
  
  if (gamma + delta > 100.0) {
    digamma_gamma_delta_plus_1 = std::log(gamma + delta + 1.0) - 1.0 / (2.0 * (gamma + delta + 1.0));
  } else {
    digamma_gamma_delta_plus_1 = R::digamma(gamma + delta + 1.0);
  }
  
  // Initialize accumulators
  double sum_log_x = 0.0;
  double sum_log_v = 0.0;
  double sum_term_lambda = 0.0;
  
  for (int i = 0; i < n; i++) {
    double xi = x(i);
    
    // Handle boundary values
    xi = std::max(eps, std::min(1.0 - eps, xi));
    
    double log_xi = std::log(xi);
    sum_log_x += log_xi;
    
    // Calculate x^λ stably
    double x_lambda;
    if (lambda > 100.0 || lambda * std::abs(log_xi) > 1.0) {
      x_lambda = safe_exp(lambda * log_xi);
    } else {
      x_lambda = std::pow(xi, lambda);
    }
    
    // Calculate v = 1-x^λ with precision
    double v;
    if (x_lambda > 0.9995) {
      v = -std::expm1(lambda * log_xi);
    } else {
      v = 1.0 - x_lambda;
    }
    v = std::max(v, eps);
    double log_v = safe_log(v);
    sum_log_v += log_v;
    
    // Calculate term for λ gradient: (x^λ * log(x)) / (1-x^λ)
    double lambda_term = (x_lambda * log_xi) / v;
    lambda_term = std::min(std::max(lambda_term, -1e6), 1e6);
    sum_term_lambda += lambda_term;
  }
  
  // =========================================================================
  // Compute gradient of LOG-LIKELIHOOD ℓ
  // =========================================================================
  
  // ∂ℓ/∂γ = -n[ψ(γ) - ψ(γ+δ+1)] + λ·Σlog(x)
  double d_gamma = -n * (digamma_gamma - digamma_gamma_delta_plus_1) + lambda * sum_log_x;
  
  // ∂ℓ/∂δ = -n[ψ(δ+1) - ψ(γ+δ+1)] + Σlog(1-x^λ)
  double d_delta = -n * (digamma_delta_plus_1 - digamma_gamma_delta_plus_1) + sum_log_v;
  
  // ∂ℓ/∂λ = n/λ + γ·Σlog(x) - δ·Σ[x^λ*log(x)/(1-x^λ)]
  double d_lambda = n / lambda + gamma * sum_log_x - delta * sum_term_lambda;
  
  // =========================================================================
  // Return NEGATIVE gradient (for minimization of negative log-likelihood)
  // =========================================================================
  grad[0] = -d_gamma;
  grad[1] = -d_delta;
  grad[2] = -d_lambda;
  
  return grad;
}


// ============================================================================
// HESSIAN OF NEGATIVE LOG-LIKELIHOOD
// ============================================================================

/**
 * @brief Hessian Matrix of Negative Log-Likelihood for BP/McDonald Distribution
 * 
 * Computes the Hessian matrix (matrix of second partial derivatives) of
 * the negative log-likelihood for standard error estimation and
 * optimization algorithms.
 * 
 * @param par Parameter vector of length 3: (γ, δ, λ)
 * @param data Vector of observations (must be in (0,1))
 * 
 * @return NumericMatrix of dimension 3×3 containing the Hessian
 * 
 * @details
 * Computes analytical second derivatives. The Hessian is symmetric.
 * Parameter ordering: (γ, δ, λ) → indices (0, 1, 2).
 * 
 * Returns NaN matrix for invalid inputs.
 * 
 * @note Exported as .hsmc_cpp for internal package use
 */
// [[Rcpp::export(.hsmc_cpp)]]
Rcpp::NumericMatrix hsmc(const Rcpp::NumericVector& par, const Rcpp::NumericVector& data) {
  // Initialize NaN matrix for error cases
  Rcpp::NumericMatrix nanHess(3, 3);
  nanHess.fill(R_NaN);
  
  // Validate parameter vector length
  if (par.size() < 3) {
    return nanHess;
  }
  
  // Extract parameters
  double gamma = par[0];
  double delta = par[1];
  double lambda = par[2];
  
  // Validate parameters using consistent checker
  if (!check_bp_pars(gamma, delta, lambda)) {
    return nanHess;
  }
  
  // Convert and validate data
  arma::vec x = Rcpp::as<arma::vec>(data);
  if (x.n_elem < 1 || arma::any(x <= 0) || arma::any(x >= 1)) {
    return nanHess;
  }
  
  int n = x.n_elem;
  Rcpp::NumericMatrix hess(3, 3);
  
  // Numerical stability constants
  const double eps = 1e-10;
  const double max_contrib = 1e6;
  
  // Compute trigamma values stably
  double trigamma_gamma, trigamma_delta_plus_1, trigamma_gamma_plus_delta_plus_1;
  
  if (gamma > 100.0) {
    trigamma_gamma = 1.0 / gamma + 1.0 / (2.0 * gamma * gamma);
  } else {
    trigamma_gamma = R::trigamma(gamma);
  }
  
  if (delta > 100.0) {
    trigamma_delta_plus_1 = 1.0 / (delta + 1.0) + 1.0 / (2.0 * (delta + 1.0) * (delta + 1.0));
  } else {
    trigamma_delta_plus_1 = R::trigamma(delta + 1.0);
  }
  
  if (gamma + delta > 100.0) {
    double z = gamma + delta + 1.0;
    trigamma_gamma_plus_delta_plus_1 = 1.0 / z + 1.0 / (2.0 * z * z);
  } else {
    trigamma_gamma_plus_delta_plus_1 = R::trigamma(gamma + delta + 1.0);
  }
  
  // Initialize accumulators for data-dependent terms
  double sum_log_x = 0.0;
  double sum_x_lambda_log_x_div_v = 0.0;
  double sum_lambda_term = 0.0;
  
  for (int i = 0; i < n; i++) {
    double xi = x(i);
    
    // Handle boundary values
    xi = std::max(eps, std::min(1.0 - eps, xi));
    
    double log_xi = safe_log(xi);
    sum_log_x += log_xi;
    
    // Calculate x^λ stably
    double x_lambda;
    if (lambda > 100.0 || lambda * std::abs(log_xi) > 1.0) {
      double log_x_lambda = lambda * log_xi;
      x_lambda = safe_exp(log_x_lambda);
    } else {
      x_lambda = std::pow(xi, lambda);
    }
    
    // Calculate v = 1-x^λ with precision
    double v;
    if (x_lambda > 0.9995) {
      v = -std::expm1(lambda * log_xi);
    } else {
      v = 1.0 - x_lambda;
    }
    v = std::max(v, eps);
    
    // Term for H[δ,λ]: Σ[x^λ*log(x)/(1-x^λ)]
    double term1 = (x_lambda * log_xi) / v;
    term1 = std::min(std::max(term1, -max_contrib), max_contrib);
    sum_x_lambda_log_x_div_v += term1;
    
    // Term for H[λ,λ]: Σ[x^λ*(log(x))²/(1-x^λ)² * (1-x^λ+x^λ)]
    //                = Σ[x^λ*(log(x))²/(1-x^λ)²]
    double log_xi_squared = log_xi * log_xi;
    double v_squared = v * v;
    
    double lambda_term = x_lambda * log_xi_squared / v_squared;
    lambda_term = std::min(std::max(lambda_term, -max_contrib), max_contrib);
    sum_lambda_term += lambda_term;
  }
  
  // =========================================================================
  // Compute Hessian of LOG-LIKELIHOOD ℓ (not negative log-likelihood)
  // =========================================================================
  
  // H[γ,γ] = ∂²ℓ/∂γ² = -n[ψ'(γ) - ψ'(γ+δ+1)]
  double H_gamma_gamma = -n * (trigamma_gamma - trigamma_gamma_plus_delta_plus_1);
  
  // H[γ,δ] = ∂²ℓ/∂γ∂δ = n·ψ'(γ+δ+1)
  double H_gamma_delta = n * trigamma_gamma_plus_delta_plus_1;
  
  // H[γ,λ] = ∂²ℓ/∂γ∂λ = Σlog(x)
  double H_gamma_lambda = sum_log_x;
  
  // H[δ,δ] = ∂²ℓ/∂δ² = -n[ψ'(δ+1) - ψ'(γ+δ+1)]
  double H_delta_delta = -n * (trigamma_delta_plus_1 - trigamma_gamma_plus_delta_plus_1);
  
  // H[δ,λ] = ∂²ℓ/∂δ∂λ = -Σ[x^λ*log(x)/(1-x^λ)]
  double H_delta_lambda = -sum_x_lambda_log_x_div_v;
  
  // H[λ,λ] = ∂²ℓ/∂λ² = -n/λ² - δ·Σ[x^λ*(log(x))²/(1-x^λ)²]
  double H_lambda_lambda = -n / (lambda * lambda) - delta * sum_lambda_term;
  
  // =========================================================================
  // Fill the Hessian matrix for NEGATIVE log-likelihood: -H
  // =========================================================================
  hess(0, 0) = -H_gamma_gamma;
  hess(0, 1) = hess(1, 0) = -H_gamma_delta;
  hess(0, 2) = hess(2, 0) = -H_gamma_lambda;
  hess(1, 1) = -H_delta_delta;
  hess(1, 2) = hess(2, 1) = -H_delta_lambda;
  hess(2, 2) = -H_lambda_lambda;
  
  return hess;
}

