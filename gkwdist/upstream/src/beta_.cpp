/**
 * @file beta.cpp
 * @brief Beta Distribution Functions (GKw Parameterization)
 * 
 * @details
 * This file implements the full suite of distribution functions for the
 * two-parameter Beta distribution using the GKw family parameterization.
 * This is obtained from the Generalized Kumaraswamy distribution by setting
 * α = 1, β = 1, and λ = 1.
 * 
 * **Relationship to GKw:**
 * \deqn{Beta(\gamma, \delta) = GKw(1, 1, \gamma, \delta, 1)}
 * 
 * **Important Parameterization Note:**
 * Unlike the standard Beta(a, b) parameterization, this implementation uses:
 * \deqn{Beta_{GKw}(\gamma, \delta) \equiv Beta_{standard}(\gamma, \delta+1)}
 * 
 * This parameterization ensures consistency across the entire GKw family of
 * distributions where the second shape parameter can be zero (δ ≥ 0).
 * 
 * The Beta distribution has probability density function:
 * \deqn{
 *   f(x; \gamma, \delta) = 
 *   \frac{1}{B(\gamma, \delta+1)} x^{\gamma-1} (1-x)^\delta
 * }
 * for \eqn{x \in (0,1)}, where \eqn{B(\cdot,\cdot)} is the Beta function.
 * 
 * **Derivation from GKw:**
 * Setting α=1, β=1, λ=1 in the GKw PDF:
 * - \eqn{x^{\alpha-1} = x^0 = 1}
 * - \eqn{(1-x^\alpha)^{\beta-1} = (1-x)^0 = 1}
 * - \eqn{[1-(1-x^\alpha)^\beta] = [1-(1-x)] = x}
 * - \eqn{[1-(1-x^\alpha)^\beta]^{\gamma\lambda-1} = x^{\gamma-1}}
 * - \eqn{\{1-[1-(1-x^\alpha)^\beta]^\lambda\}^\delta = (1-x)^\delta}
 * - The constant becomes: \eqn{1/(B(\gamma,\delta+1))}
 * 
 * The cumulative distribution function is:
 * \deqn{
 *   F(x) = I_x(\gamma, \delta+1)
 * }
 * where \eqn{I_x(a,b)} is the regularized incomplete Beta function.
 * 
 * The quantile function (inverse CDF) is:
 * \deqn{
 *   Q(p) = I_p^{-1}(\gamma, \delta+1)
 * }
 * 
 * **Parameter Constraints:**
 * - \eqn{\gamma > 0} (shape parameter)
 * - \eqn{\delta > 0} (shape parameter - note: must be POSITIVE for Beta)
 * 
 * **Note on δ = 0:**
 * While the GKw family allows δ ≥ 0, the standard Beta distribution requires
 * δ > 0. The case δ = 0 in GKw corresponds to the Power function distribution,
 * not Beta.
 * 
 * **Special Cases:**
 * | Distribution | Condition | Relation |
 * |--------------|-----------|----------|
 * | Uniform(0,1) | \eqn{\gamma = \delta = 1} | Beta(1, 2) in this param. |
 * | Arc-sine | \eqn{\gamma = \delta = 0.5} | Beta(0.5, 1.5) in this param. |
 * 
 * **Random Variate Generation:**
 * Uses R's built-in rbeta with adjusted parameters:
 * \deqn{X \sim rbeta(\gamma, \delta+1)}
 * 
 * **Numerical Stability:**
 * All computations use log-space arithmetic and numerically stable helper
 * functions from utils.h to prevent overflow/underflow.
 * 
 * **Implemented Functions:**
 * - dbeta_(): Probability density function (PDF)
 * - pbeta_(): Cumulative distribution function (CDF)
 * - qbeta_(): Quantile function (inverse CDF)
 * - rbeta_(): Random variate generation
 * - llbeta(): Negative log-likelihood for MLE
 * - grbeta(): Gradient of negative log-likelihood
 * - hsbeta(): Hessian of negative log-likelihood
 * 
 * @author Lopes, J. E.
 * @date 2025-01-07
 * 
 * @references
 * Johnson, N. L., Kotz, S., & Balakrishnan, N. (1995). Continuous Univariate
 * Distributions, Volume 2 (2nd ed.). Wiley.
 * 
 * @see gkw.cpp for the parent distribution
 * @see utils.h for numerical stability functions and parameter validators
 * 
 * @note All functions use R's vectorization conventions with parameter recycling.
 * @note Thread-safe: No global state is modified.
 * @note Function names have trailing underscore to avoid conflicts with R's built-in beta functions.
 */

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "utils.h"


// ============================================================================
// PROBABILITY DENSITY FUNCTION
// ============================================================================

/**
 * @brief Probability Density Function of the Beta Distribution (GKw Parameterization)
 * 
 * Computes the density (or log-density) for the Beta distribution
 * at specified quantiles using GKw-consistent parameterization.
 * 
 * @param x Vector of quantiles (values in (0,1))
 * @param gamma Shape parameter vector (γ > 0)
 * @param delta Shape parameter vector (δ > 0)
 * @param log_prob If TRUE, returns log-density; otherwise returns density
 * 
 * @return NumericVector of density values (or log-density if log_prob=TRUE)
 * 
 * @details
 * The log-density is computed as:
 * \deqn{
 *   \log f(x) = -\log B(\gamma, \delta+1)
 *   + (\gamma-1)\log(x) + \delta\log(1-x)
 * }
 * 
 * This corresponds to Beta(γ, δ+1) in standard parameterization.
 * 
 * @note Exported as .dbeta_cpp for internal package use
 */
// [[Rcpp::export(.dbeta_cpp)]]
Rcpp::NumericVector dbeta_(
    const arma::vec& x,
    const Rcpp::NumericVector& gamma,
    const Rcpp::NumericVector& delta,
    bool log_prob = false
) {
  // Convert R vectors to Armadillo vectors
  arma::vec g_vec(gamma.begin(), gamma.size());
  arma::vec d_vec(delta.begin(), delta.size());
  
  // Determine output length for recycling
  size_t N = std::max({x.n_elem, g_vec.n_elem, d_vec.n_elem});
  
  // Initialize result with appropriate default
  arma::vec out(N);
  out.fill(log_prob ? R_NegInf : 0.0);
  
  for (size_t i = 0; i < N; i++) {
    // Extract recycled parameters
    double g = g_vec[i % g_vec.n_elem];
    double d = d_vec[i % d_vec.n_elem];
    double xx = x[i % x.n_elem];
    
    // Validate parameters (Beta requires δ > 0, not δ ≥ 0)
    if (!check_beta_pars(g, d)) {
      continue;
    }
    
    // Check support: x must be in (0, 1)
    if (xx <= 0.0 || xx >= 1.0 || !R_finite(xx)) {
      continue;
    }
    
    // ---- Log-space computation of density ----
    
    // Normalization constant: -log(B(γ, δ+1))
    double lB = R::lbeta(g, d + 1.0);
    
    // Compute log(x) and log(1-x)
    double lx = safe_log(xx);
    double one_minus_x = 1.0 - xx;
    
    if (one_minus_x <= 0.0) {
      continue;
    }
    double log_1_minus_x = safe_log(one_minus_x);
    
    // Assemble log-density:
    // log(f) = -log(B) + (γ-1)*log(x) + δ*log(1-x)
    double log_pdf = (g - 1.0) * lx + d * log_1_minus_x - lB;
    
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
 * @brief Cumulative Distribution Function of the Beta Distribution (GKw Parameterization)
 * 
 * Computes the cumulative probability for the Beta distribution
 * at specified quantiles using GKw-consistent parameterization.
 * 
 * @param q Vector of quantiles
 * @param gamma Shape parameter vector (γ > 0)
 * @param delta Shape parameter vector (δ > 0)
 * @param lower_tail If TRUE, returns P(X ≤ q); otherwise P(X > q)
 * @param log_p If TRUE, returns log-probability
 * 
 * @return NumericVector of cumulative probabilities
 * 
 * @details
 * The CDF is computed as:
 * \deqn{F(x) = I_x(\gamma, \delta+1)}
 * where \eqn{I_x(a,b)} is the regularized incomplete Beta function.
 * 
 * Uses R's pbeta with adjusted parameters: pbeta(x, γ, δ+1).
 * 
 * @note Exported as .pbeta_cpp for internal package use
 */
// [[Rcpp::export(.pbeta_cpp)]]
Rcpp::NumericVector pbeta_(
    const arma::vec& q,
    const Rcpp::NumericVector& gamma,
    const Rcpp::NumericVector& delta,
    bool lower_tail = true,
    bool log_p = false
) {
  // Convert R vectors to Armadillo vectors
  arma::vec g_vec(gamma.begin(), gamma.size());
  arma::vec d_vec(delta.begin(), delta.size());
  
  // Determine output length for recycling
  size_t N = std::max({q.n_elem, g_vec.n_elem, d_vec.n_elem});
  
  arma::vec out(N);
  
  for (size_t i = 0; i < N; i++) {
    // Extract recycled parameters
    double g = g_vec[i % g_vec.n_elem];
    double d = d_vec[i % d_vec.n_elem];
    double qq = q[i % q.n_elem];
    
    // Validate parameters
    if (!check_beta_pars(g, d)) {
      out(i) = NA_REAL;
      continue;
    }
    
    // Handle boundary: q ≤ 0
    if (!R_finite(qq) || qq <= 0.0) {
      double v0 = lower_tail ? 0.0 : 1.0;
      out(i) = log_p ? safe_log(v0) : v0;
      continue;
    }
    
    // Handle boundary: q ≥ 1
    if (qq >= 1.0) {
      double v1 = lower_tail ? 1.0 : 0.0;
      out(i) = log_p ? safe_log(v1) : v1;
      continue;
    }
    
    // ---- Compute CDF via R's pbeta with adjusted parameters ----
    // Beta_GKw(γ, δ) = Beta_standard(γ, δ+1)
    double val = R::pbeta(qq, g, d + 1.0, lower_tail, false);
    
    // Apply log transformation if requested
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
 * @brief Quantile Function (Inverse CDF) of the Beta Distribution (GKw Parameterization)
 * 
 * Computes quantiles for the Beta distribution given probability values
 * using GKw-consistent parameterization.
 * 
 * @param p Vector of probabilities (values in [0,1])
 * @param gamma Shape parameter vector (γ > 0)
 * @param delta Shape parameter vector (δ > 0)
 * @param lower_tail If TRUE, probabilities are P(X ≤ x); otherwise P(X > x)
 * @param log_p If TRUE, probabilities are given as log(p)
 * 
 * @return NumericVector of quantiles
 * 
 * @details
 * The quantile function inverts the CDF:
 * \deqn{Q(p) = I_p^{-1}(\gamma, \delta+1)}
 * 
 * Uses R's qbeta with adjusted parameters: qbeta(p, γ, δ+1).
 * 
 * @note Exported as .qbeta_cpp for internal package use
 */
// [[Rcpp::export(.qbeta_cpp)]]
Rcpp::NumericVector qbeta_(
    const arma::vec& p,
    const Rcpp::NumericVector& gamma,
    const Rcpp::NumericVector& delta,
    bool lower_tail = true,
    bool log_p = false
) {
  // Convert R vectors to Armadillo vectors
  arma::vec g_vec(gamma.begin(), gamma.size());
  arma::vec d_vec(delta.begin(), delta.size());
  
  // Determine output length for recycling
  size_t N = std::max({p.n_elem, g_vec.n_elem, d_vec.n_elem});
  
  arma::vec out(N);
  
  for (size_t i = 0; i < N; i++) {
    // Extract recycled parameters
    double g = g_vec[i % g_vec.n_elem];
    double d = d_vec[i % d_vec.n_elem];
    double pp = p[i % p.n_elem];
    
    // Validate parameters
    if (!check_beta_pars(g, d)) {
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
    
    // ---- Compute quantile via R's qbeta with adjusted parameters ----
    // Beta_GKw(γ, δ) = Beta_standard(γ, δ+1)
    double val = R::qbeta(pp, g, d + 1.0, true, false);
    
    out(i) = val;
  }
  
  return Rcpp::NumericVector(out.memptr(), out.memptr() + out.n_elem);
}


// ============================================================================
// RANDOM NUMBER GENERATION
// ============================================================================

/**
 * @brief Random Variate Generation for the Beta Distribution (GKw Parameterization)
 * 
 * Generates random samples from the Beta distribution
 * using GKw-consistent parameterization.
 * 
 * @param n Number of random variates to generate
 * @param gamma Shape parameter vector (γ > 0)
 * @param delta Shape parameter vector (δ > 0)
 * 
 * @return NumericVector of n random variates from Beta distribution
 * 
 * @details
 * Uses R's built-in rbeta with adjusted parameters:
 * \deqn{X \sim rbeta(\gamma, \delta+1)}
 * 
 * @note Exported as .rbeta_cpp for internal package use
 */
// [[Rcpp::export(.rbeta_cpp)]]
Rcpp::NumericVector rbeta_(
    int n,
    const Rcpp::NumericVector& gamma,
    const Rcpp::NumericVector& delta
) {
  if (n <= 0) {
    Rcpp::stop("rbeta_: n must be positive");
  }
  
  // Convert R vectors to Armadillo vectors
  arma::vec g_vec(gamma.begin(), gamma.size());
  arma::vec d_vec(delta.begin(), delta.size());
  
  arma::vec out(n);
  
  for (int i = 0; i < n; i++) {
    // Extract recycled parameters (direct modulo, no intermediate variable)
    double g = g_vec[i % g_vec.n_elem];
    double d = d_vec[i % d_vec.n_elem];
    
    // Validate parameters
    if (!check_beta_pars(g, d)) {
      out(i) = NA_REAL;
      Rcpp::warning("rbeta_: invalid parameters at index %d", i + 1);
      continue;
    }
    
    // Generate from Beta(γ, δ+1) using R's rbeta
    double val = R::rbeta(g, d + 1.0);
    out(i) = val;
  }
  
  return Rcpp::NumericVector(out.memptr(), out.memptr() + out.n_elem);
}


// ============================================================================
// NEGATIVE LOG-LIKELIHOOD FUNCTION
// ============================================================================

/**
 * @brief Negative Log-Likelihood for Beta Distribution (GKw Parameterization)
 * 
 * Computes the negative log-likelihood function for parameter estimation
 * via maximum likelihood.
 * 
 * @param par Parameter vector of length 2: (γ, δ)
 * @param data Vector of observations (must be in (0,1))
 * 
 * @return Negative log-likelihood value (scalar)
 * 
 * @details
 * The log-likelihood for n observations is:
 * \deqn{
 *   \ell(\theta) = -n\ln B(\gamma,\delta+1)
 *   + (\gamma-1)\sum\ln x_i + \delta\sum\ln(1-x_i)
 * }
 * 
 * Returns +Inf for invalid parameters or data outside (0,1).
 * 
 * @note Exported as .llbeta_cpp for internal package use
 */
// [[Rcpp::export(.llbeta_cpp)]]
double llbeta(const Rcpp::NumericVector& par, const Rcpp::NumericVector& data) {
  // Validate parameter vector length
  if (par.size() < 2) {
    return R_PosInf;
  }
  
  // Extract parameters
  double gamma = par[0];
  double delta = par[1];
  
  // Validate parameters using consistent checker
  if (!check_beta_pars(gamma, delta)) {
    return R_PosInf;
  }
  
  // Convert and validate data
  arma::vec x = Rcpp::as<arma::vec>(data);
  if (x.n_elem < 1) {
    return R_PosInf;
  }
  if (arma::any(x <= 0.0) || arma::any(x >= 1.0)) {
    return R_PosInf;
  }
  
  int n = x.n_elem;
  
  // Constant term: -n * log(B(γ, δ+1))
  double logB = R::lbeta(gamma, delta + 1.0);
  double cst = -double(n) * logB;
  
  // Term 1: (γ-1) * Σlog(x)
  arma::vec lx = vec_safe_log(x);
  double sum1 = (gamma - 1.0) * arma::sum(lx);
  
  // Term 2: δ * Σlog(1-x)
  arma::vec l1mx = vec_safe_log(1.0 - x);
  double sum2 = delta * arma::sum(l1mx);
  
  double loglike = cst + sum1 + sum2;
  
  return -loglike;
}


// ============================================================================
// GRADIENT OF NEGATIVE LOG-LIKELIHOOD
// ============================================================================

/**
 * @brief Gradient of Negative Log-Likelihood for Beta Distribution (GKw Parameterization)
 * 
 * Computes the gradient vector of the negative log-likelihood for
 * optimization-based parameter estimation.
 * 
 * @param par Parameter vector of length 2: (γ, δ)
 * @param data Vector of observations (must be in (0,1))
 * 
 * @return NumericVector of length 2 containing partial derivatives
 *         with respect to (γ, δ)
 * 
 * @details
 * The gradient components are:
 * - ∂ℓ/∂γ = -n[ψ(γ) - ψ(γ+δ+1)] + Σlog(x)
 * - ∂ℓ/∂δ = -n[ψ(δ+1) - ψ(γ+δ+1)] + Σlog(1-x)
 * 
 * For negative log-likelihood, we return the negation.
 * 
 * @note Exported as .grbeta_cpp for internal package use
 */
// [[Rcpp::export(.grbeta_cpp)]]
Rcpp::NumericVector grbeta(const Rcpp::NumericVector& par, const Rcpp::NumericVector& data) {
  // Validate parameter vector length
  if (par.size() < 2) {
    return Rcpp::NumericVector(2, R_NaN);
  }
  
  // Extract parameters
  double gamma = par[0];
  double delta = par[1];
  
  // Validate parameters using consistent checker
  if (!check_beta_pars(gamma, delta)) {
    return Rcpp::NumericVector(2, R_NaN);
  }
  
  // Convert and validate data
  arma::vec x = Rcpp::as<arma::vec>(data);
  if (x.n_elem < 1) {
    return Rcpp::NumericVector(2, R_NaN);
  }
  if (arma::any(x <= 0.0) || arma::any(x >= 1.0)) {
    return Rcpp::NumericVector(2, R_NaN);
  }
  
  int n = x.n_elem;
  Rcpp::NumericVector grad(2, 0.0);
  
  // Calculate digamma terms with correct parameterization
  double dig_g = R::digamma(gamma);
  double dig_d = R::digamma(delta + 1.0);
  double dig_gd = R::digamma(gamma + delta + 1.0);
  
  // Sum of log terms
  arma::vec lx = vec_safe_log(x);
  arma::vec l1mx = vec_safe_log(1.0 - x);
  double sum_lx = arma::sum(lx);
  double sum_l1mx = arma::sum(l1mx);
  
  // Gradient of log-likelihood
  double d_gamma = -n * (dig_g - dig_gd) + sum_lx;
  double d_delta = -n * (dig_d - dig_gd) + sum_l1mx;
  
  // Return NEGATIVE gradient (for minimization of negative log-likelihood)
  grad[0] = -d_gamma;
  grad[1] = -d_delta;
  
  return grad;
}


// ============================================================================
// HESSIAN OF NEGATIVE LOG-LIKELIHOOD
// ============================================================================

/**
 * @brief Hessian Matrix of Negative Log-Likelihood for Beta Distribution (GKw Parameterization)
 * 
 * Computes the Hessian matrix (matrix of second partial derivatives) of
 * the negative log-likelihood for standard error estimation and
 * optimization algorithms.
 * 
 * @param par Parameter vector of length 2: (γ, δ)
 * @param data Vector of observations (must be in (0,1))
 * 
 * @return NumericMatrix of dimension 2×2 containing the Hessian
 * 
 * @details
 * Computes analytical second derivatives. The Hessian is symmetric.
 * Parameter ordering: (γ, δ) → indices (0, 1).
 * 
 * The Hessian components (of log-likelihood ℓ) are:
 * - H[γ,γ] = ∂²ℓ/∂γ² = -n[ψ₁(γ) - ψ₁(γ+δ+1)]
 * - H[δ,δ] = ∂²ℓ/∂δ² = -n[ψ₁(δ+1) - ψ₁(γ+δ+1)]
 * - H[γ,δ] = ∂²ℓ/∂γ∂δ = n*ψ₁(γ+δ+1)
 * 
 * For negative log-likelihood, we return the negation.
 * 
 * Returns NaN matrix for invalid inputs.
 * 
 * @note Exported as .hsbeta_cpp for internal package use
 */
// [[Rcpp::export(.hsbeta_cpp)]]
Rcpp::NumericMatrix hsbeta(const Rcpp::NumericVector& par, const Rcpp::NumericVector& data) {
  // Initialize NaN matrix for error cases
  Rcpp::NumericMatrix nanHess(2, 2);
  nanHess.fill(R_NaN);
  
  // Validate parameter vector length
  if (par.size() < 2) {
    return nanHess;
  }
  
  // Extract parameters
  double gamma = par[0];
  double delta = par[1];
  
  // Validate parameters using consistent checker
  if (!check_beta_pars(gamma, delta)) {
    return nanHess;
  }
  
  // Convert and validate data
  arma::vec x = Rcpp::as<arma::vec>(data);
  if (x.n_elem < 1) {
    return nanHess;
  }
  if (arma::any(x <= 0.0) || arma::any(x >= 1.0)) {
    return nanHess;
  }
  
  int n = x.n_elem;
  Rcpp::NumericMatrix hess(2, 2);
  
  // Calculate trigamma terms with correct parameterization
  double trig_g = R::trigamma(gamma);
  double trig_d = R::trigamma(delta + 1.0);
  double trig_gd = R::trigamma(gamma + delta + 1.0);
  
  // Hessian components of log-likelihood
  double h_gamma_gamma = -n * (trig_g - trig_gd);
  double h_delta_delta = -n * (trig_d - trig_gd);
  double h_gamma_delta = n * trig_gd;
  
  // Fill the Hessian matrix (symmetric) - NEGATE for negative log-likelihood
  hess(0, 0) = -h_gamma_gamma;
  hess(1, 1) = -h_delta_delta;
  hess(0, 1) = hess(1, 0) = -h_gamma_delta;
  
  return hess;
}
