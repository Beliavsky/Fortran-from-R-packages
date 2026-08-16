/**
 * @file ekw.cpp
 * @brief Exponentiated Kumaraswamy (EKw) Distribution Functions
 * 
 * @details
 * This file implements the full suite of distribution functions for the
 * three-parameter Exponentiated Kumaraswamy (EKw) distribution, which is a
 * sub-family of the Generalized Kumaraswamy (GKw) distribution obtained
 * by setting γ = 1 and δ = 0.
 * 
 * **Relationship to GKw:**
 * \deqn{EKw(\alpha, \beta, \lambda) = GKw(\alpha, \beta, 1, 0, \lambda)}
 * 
 * The EKw distribution has probability density function:
 * \deqn{
 *   f(x; \alpha, \beta, \lambda) = 
 *   \lambda \alpha \beta x^{\alpha-1} (1-x^\alpha)^{\beta-1}
 *   [1-(1-x^\alpha)^\beta]^{\lambda-1}
 * }
 * for \eqn{x \in (0,1)}.
 * 
 * **Derivation from GKw:**
 * Setting γ=1 and δ=0 in the GKw PDF:
 * - The Beta function term becomes: \eqn{B(1, 0+1) = B(1,1) = 1}
 * - The exponent on the outer bracket: \eqn{\gamma\lambda - 1 = 1\cdot\lambda - 1 = \lambda - 1}
 * - The final term: \eqn{\{1-[...]\}^\delta = \{1-[...]\}^0 = 1}
 * 
 * The cumulative distribution function is:
 * \deqn{
 *   F(x) = [1-(1-x^\alpha)^\beta]^\lambda
 * }
 * 
 * The quantile function (inverse CDF) is:
 * \deqn{
 *   Q(p) = \left\{1 - \left[1 - p^{1/\lambda}\right]^{1/\beta}\right\}^{1/\alpha}
 * }
 * 
 * **Parameter Constraints:**
 * - \eqn{\alpha > 0} (shape parameter)
 * - \eqn{\beta > 0} (shape parameter)
 * - \eqn{\lambda > 0} (exponentiation parameter)
 * 
 * **Special Cases:**
 * | Distribution | Condition | Relation |
 * |--------------|-----------|----------|
 * | Kumaraswamy (Kw) | \eqn{\lambda = 1} | Standard Kumaraswamy |
 * | Generalized Rayleigh | \eqn{\alpha = 2, \beta = 1} | EKw(2, 1, λ) |
 * 
 * **Random Variate Generation:**
 * Uses inverse transform method:
 * 1. Generate \eqn{U \sim Uniform(0,1)}
 * 2. Return \eqn{X = Q(U) = \{1 - [1 - U^{1/\lambda}]^{1/\beta}\}^{1/\alpha}}
 * 
 * **Numerical Stability:**
 * Special attention is given to λ ≈ 1, which can cause numerical cancellation.
 * All computations use log-space arithmetic and numerically stable helper
 * functions from utils.h.
 * 
 * **Implemented Functions:**
 * - dekw(): Probability density function (PDF)
 * - pekw(): Cumulative distribution function (CDF)
 * - qekw(): Quantile function (inverse CDF)
 * - rekw(): Random variate generation
 * - llekw(): Negative log-likelihood for MLE
 * - grekw(): Gradient of negative log-likelihood
 * - hsekw(): Hessian of negative log-likelihood
 * 
 * @author Lopes, J. E.
 * @date 2025-01-07
 * 
 * @see gkw.cpp for the parent distribution
 * @see kkw.cpp for the sister distribution with δ free
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
 * @brief Probability Density Function of the EKw Distribution
 * 
 * Computes the density (or log-density) for the Exponentiated Kumaraswamy
 * distribution at specified quantiles.
 * 
 * @param x Vector of quantiles (values in (0,1))
 * @param alpha Shape parameter vector (α > 0)
 * @param beta Shape parameter vector (β > 0)
 * @param lambda Exponentiation parameter vector (λ > 0)
 * @param log_prob If TRUE, returns log-density; otherwise returns density
 * 
 * @return NumericVector of density values (or log-density if log_prob=TRUE)
 * 
 * @details
 * The log-density is computed as:
 * \deqn{
 *   \log f(x) = \log(\lambda) + \log(\alpha) + \log(\beta)
 *   + (\alpha-1)\log(x) + (\beta-1)\log(1-x^\alpha)
 *   + (\lambda-1)\log(1-(1-x^\alpha)^\beta)
 * }
 * 
 * @note Exported as .dekw_cpp for internal package use
 */
// [[Rcpp::export(.dekw_cpp)]]
Rcpp::NumericVector dekw(
    const arma::vec& x,
    const Rcpp::NumericVector& alpha,
    const Rcpp::NumericVector& beta,
    const Rcpp::NumericVector& lambda,
    bool log_prob = false
) {
  // Convert R vectors to Armadillo vectors
  arma::vec a_vec(alpha.begin(), alpha.size());
  arma::vec b_vec(beta.begin(), beta.size());
  arma::vec l_vec(lambda.begin(), lambda.size());
  
  // Determine output length for recycling
  size_t N = std::max({x.n_elem, a_vec.n_elem, b_vec.n_elem, l_vec.n_elem});
  
  // Initialize result with appropriate default
  arma::vec out(N);
  out.fill(log_prob ? R_NegInf : 0.0);
  
  for (size_t i = 0; i < N; i++) {
    // Extract recycled parameters
    double a = a_vec[i % a_vec.n_elem];
    double b = b_vec[i % b_vec.n_elem];
    double l = l_vec[i % l_vec.n_elem];
    double xx = x[i % x.n_elem];
    
    // Validate parameters
    if (!check_ekw_pars(a, b, l)) {
      continue;
    }
    
    // Check support: x must be in (0, 1)
    if (xx <= 0.0 || xx >= 1.0 || !R_finite(xx)) {
      continue;
    }
    
    // ---- Log-space computation of density ----
    
    // Normalization constant: log(λαβ)
    double ll = safe_log(l);
    double la = safe_log(a);
    double lb = safe_log(b);
    double lx = safe_log(xx);
    
    // Compute log(x^α) = α * log(x)
    double log_xalpha = a * lx;
    
    // Compute log(1 - x^α) using stable log1mexp
    double log_v = log1mexp(log_xalpha);
    if (!R_finite(log_v)) {
      continue;
    }
    
    // Term 1: (β-1) * log(1 - x^α)
    double term2 = (b - 1.0) * log_v;
    
    // Compute log((1-x^α)^β) = β * log(1-x^α)
    double log_v_beta = b * log_v;
    
    // Compute log(1 - (1-x^α)^β) = log(w) using log1mexp
    double log_w = log1mexp(log_v_beta);
    if (!R_finite(log_w)) {
      continue;
    }
    
    // Term 2: (λ-1) * log(w)
    double term3 = (l - 1.0) * log_w;
    
    // Assemble log-density:
    // log(f) = log(λαβ) + (α-1)*log(x) + (β-1)*log(v) + (λ-1)*log(w)
    double log_pdf = ll + la + lb + (a - 1.0) * lx + term2 + term3;
    
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
 * @brief Cumulative Distribution Function of the EKw Distribution
 * 
 * Computes the cumulative probability for the Exponentiated Kumaraswamy
 * distribution at specified quantiles.
 * 
 * @param q Vector of quantiles
 * @param alpha Shape parameter vector (α > 0)
 * @param beta Shape parameter vector (β > 0)
 * @param lambda Exponentiation parameter vector (λ > 0)
 * @param lower_tail If TRUE, returns P(X ≤ q); otherwise P(X > q)
 * @param log_p If TRUE, returns log-probability
 * 
 * @return NumericVector of cumulative probabilities
 * 
 * @details
 * The CDF is computed as:
 * \deqn{F(x) = [1-(1-x^\alpha)^\beta]^\lambda}
 * 
 * @note Exported as .pekw_cpp for internal package use
 */
// [[Rcpp::export(.pekw_cpp)]]
Rcpp::NumericVector pekw(
    const arma::vec& q,
    const Rcpp::NumericVector& alpha,
    const Rcpp::NumericVector& beta,
    const Rcpp::NumericVector& lambda,
    bool lower_tail = true,
    bool log_p = false
) {
  // Convert R vectors to Armadillo vectors
  arma::vec a_vec(alpha.begin(), alpha.size());
  arma::vec b_vec(beta.begin(), beta.size());
  arma::vec l_vec(lambda.begin(), lambda.size());
  
  // Determine output length for recycling
  size_t N = std::max({q.n_elem, a_vec.n_elem, b_vec.n_elem, l_vec.n_elem});
  
  arma::vec out(N);
  
  for (size_t i = 0; i < N; i++) {
    // Extract recycled parameters
    double a = a_vec[i % a_vec.n_elem];
    double b = b_vec[i % b_vec.n_elem];
    double l = l_vec[i % l_vec.n_elem];
    double xx = q[i % q.n_elem];
    
    // Validate parameters
    if (!check_ekw_pars(a, b, l)) {
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
    
    // Step 1: x^α
    double lx = safe_log(xx);
    double xalpha = safe_exp(a * lx);
    
    // Step 2: 1 - x^α
    double omx = 1.0 - xalpha;
    if (omx <= 0.0) {
      double val1 = lower_tail ? 1.0 : 0.0;
      out(i) = log_p ? safe_log(val1) : val1;
      continue;
    }
    
    // Step 3: (1 - x^α)^β
    double omx_beta = safe_pow(omx, b);
    
    // Step 4: t = 1 - (1 - x^α)^β
    double t = 1.0 - omx_beta;
    if (t <= 0.0) {
      double val0 = lower_tail ? 0.0 : 1.0;
      out(i) = log_p ? safe_log(val0) : val0;
      continue;
    }
    if (t >= 1.0) {
      double val1 = lower_tail ? 1.0 : 0.0;
      out(i) = log_p ? safe_log(val1) : val1;
      continue;
    }
    
    // Step 5: F(x) = t^λ
    double val = safe_pow(t, l);
    
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
 * @brief Quantile Function (Inverse CDF) of the EKw Distribution
 * 
 * Computes quantiles for the Exponentiated Kumaraswamy distribution
 * given probability values.
 * 
 * @param p Vector of probabilities (values in [0,1])
 * @param alpha Shape parameter vector (α > 0)
 * @param beta Shape parameter vector (β > 0)
 * @param lambda Exponentiation parameter vector (λ > 0)
 * @param lower_tail If TRUE, probabilities are P(X ≤ x); otherwise P(X > x)
 * @param log_p If TRUE, probabilities are given as log(p)
 * 
 * @return NumericVector of quantiles
 * 
 * @details
 * The quantile function inverts the CDF:
 * \deqn{Q(p) = \left\{1 - \left[1 - p^{1/\lambda}\right]^{1/\beta}\right\}^{1/\alpha}}
 * 
 * @note Exported as .qekw_cpp for internal package use
 */
// [[Rcpp::export(.qekw_cpp)]]
Rcpp::NumericVector qekw(
    const arma::vec& p,
    const Rcpp::NumericVector& alpha,
    const Rcpp::NumericVector& beta,
    const Rcpp::NumericVector& lambda,
    bool lower_tail = true,
    bool log_p = false
) {
  // Convert R vectors to Armadillo vectors
  arma::vec a_vec(alpha.begin(), alpha.size());
  arma::vec b_vec(beta.begin(), beta.size());
  arma::vec l_vec(lambda.begin(), lambda.size());
  
  // Determine output length for recycling
  size_t N = std::max({p.n_elem, a_vec.n_elem, b_vec.n_elem, l_vec.n_elem});
  
  arma::vec out(N);
  
  for (size_t i = 0; i < N; i++) {
    // Extract recycled parameters
    double a = a_vec[i % a_vec.n_elem];
    double b = b_vec[i % b_vec.n_elem];
    double l = l_vec[i % l_vec.n_elem];
    double pp = p[i % p.n_elem];
    
    // Validate parameters
    if (!check_ekw_pars(a, b, l)) {
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
    
    // Step 1: p^(1/λ)
    double step1 = safe_pow(pp, 1.0 / l);
    
    // Step 2: 1 - p^(1/λ)
    double step2 = 1.0 - step1;
    step2 = std::max(0.0, step2);
    
    // Step 3: [1 - p^(1/λ)]^(1/β)
    double step3 = safe_pow(step2, 1.0 / b);
    
    // Step 4: 1 - [1 - p^(1/λ)]^(1/β)
    double step4 = 1.0 - step3;
    step4 = std::max(0.0, step4);
    
    // Step 5: {1 - [1 - p^(1/λ)]^(1/β)}^(1/α)
    double x;
    if (a == 1.0) {
      x = step4;
    } else {
      x = safe_pow(step4, 1.0 / a);
    }
    
    // Clamp to valid support
    x = std::max(0.0, std::min(1.0, x));
    
    out(i) = x;
  }
  
  return Rcpp::NumericVector(out.memptr(), out.memptr() + out.n_elem);
}


// ============================================================================
// RANDOM NUMBER GENERATION
// ============================================================================

/**
 * @brief Random Variate Generation for the EKw Distribution
 * 
 * Generates random samples from the Exponentiated Kumaraswamy distribution
 * using the inverse transform method.
 * 
 * @param n Number of random variates to generate
 * @param alpha Shape parameter vector (α > 0)
 * @param beta Shape parameter vector (β > 0)
 * @param lambda Exponentiation parameter vector (λ > 0)
 * 
 * @return NumericVector of n random variates from EKw distribution
 * 
 * @details
 * Algorithm:
 * 1. Generate U ~ Uniform(0,1)
 * 2. Return X = Q(U) = {1 - [1 - U^(1/λ)]^(1/β)}^(1/α)
 * 
 * @note Exported as .rekw_cpp for internal package use
 */
// [[Rcpp::export(.rekw_cpp)]]
Rcpp::NumericVector rekw(
    int n,
    const Rcpp::NumericVector& alpha,
    const Rcpp::NumericVector& beta,
    const Rcpp::NumericVector& lambda
) {
  if (n <= 0) {
    Rcpp::stop("rekw: n must be positive");
  }
  
  // Convert R vectors to Armadillo vectors
  arma::vec a_vec(alpha.begin(), alpha.size());
  arma::vec b_vec(beta.begin(), beta.size());
  arma::vec l_vec(lambda.begin(), lambda.size());
  
  arma::vec out(n);
  
  for (int i = 0; i < n; i++) {
    // Extract recycled parameters (direct modulo, no intermediate variable)
    double a = a_vec[i % a_vec.n_elem];
    double b = b_vec[i % b_vec.n_elem];
    double l = l_vec[i % l_vec.n_elem];
    
    // Validate parameters
    if (!check_ekw_pars(a, b, l)) {
      out(i) = NA_REAL;
      Rcpp::warning("rekw: invalid parameters at index %d", i + 1);
      continue;
    }
    
    // Generate U ~ Uniform(0,1)
    double U = R::runif(0.0, 1.0);
    
    // Step 1: U^(1/λ)
    double step1 = safe_pow(U, 1.0 / l);
    
    // Step 2: 1 - U^(1/λ)
    double step2 = 1.0 - step1;
    step2 = std::max(0.0, step2);
    
    // Step 3: [1 - U^(1/λ)]^(1/β)
    double step3 = safe_pow(step2, 1.0 / b);
    
    // Step 4: 1 - [1 - U^(1/λ)]^(1/β)
    double step4 = 1.0 - step3;
    step4 = std::max(0.0, step4);
    
    // Step 5: {1 - [1 - U^(1/λ)]^(1/β)}^(1/α)
    double x;
    if (a == 1.0) {
      x = step4;
    } else {
      x = safe_pow(step4, 1.0 / a);
      if (!R_finite(x) || x < 0.0) x = 0.0;
      if (x > 1.0) x = 1.0;
    }
    
    out(i) = x;
  }
  
  return Rcpp::NumericVector(out.memptr(), out.memptr() + out.n_elem);
}


// ============================================================================
// NEGATIVE LOG-LIKELIHOOD FUNCTION
// ============================================================================

/**
 * @brief Negative Log-Likelihood for EKw Distribution
 * 
 * Computes the negative log-likelihood function for parameter estimation
 * via maximum likelihood.
 * 
 * @param par Parameter vector of length 3: (α, β, λ)
 * @param data Vector of observations (must be in (0,1))
 * 
 * @return Negative log-likelihood value (scalar)
 * 
 * @details
 * The log-likelihood for n observations is:
 * \deqn{
 *   \ell(\theta) = n[\ln\lambda + \ln\alpha + \ln\beta]
 *   + (\alpha-1)\sum\ln x_i + (\beta-1)\sum\ln v_i
 *   + (\lambda-1)\sum\ln w_i
 * }
 * where:
 * - \eqn{v_i = 1 - x_i^\alpha}
 * - \eqn{w_i = 1 - v_i^\beta}
 * 
 * Returns +Inf for invalid parameters or data outside (0,1).
 * 
 * **Special handling for λ ≈ 1:**
 * When λ is very close to 1, the term (λ-1)*log(w) can suffer from
 * catastrophic cancellation. Special care is taken in this regime.
 * 
 * @note Exported as .llekw_cpp for internal package use
 */
// [[Rcpp::export(.llekw_cpp)]]
double llekw(const Rcpp::NumericVector& par, const Rcpp::NumericVector& data) {
  // Validate parameter vector length
  if (par.size() < 3) return R_PosInf;
  
  // Extract parameters
  double alpha = par[0];
  double beta = par[1];
  double lambda = par[2];
  
  // Validate parameters using consistent checker
  if (!check_ekw_pars(alpha, beta, lambda)) return R_PosInf;
  
  // Convert and validate data
  arma::vec x = Rcpp::as<arma::vec>(data);
  if (x.n_elem < 1) return R_PosInf;
  if (arma::any(x <= 0.0) || arma::any(x >= 1.0)) return R_PosInf;
  
  int n = x.n_elem;
  
  // Constant term: n * [log(λ) + log(α) + log(β)]
  double log_alpha = safe_log(alpha);
  double log_beta = safe_log(beta);
  double log_lambda = safe_log(lambda);
  double const_term = n * (log_lambda + log_alpha + log_beta);
  
  // Initialize accumulators
  double sum_term1 = 0.0;  // (α-1) * Σlog(x)
  double sum_term2 = 0.0;  // (β-1) * Σlog(v)
  double sum_term3 = 0.0;  // (λ-1) * Σlog(w)
  
  for (int i = 0; i < n; i++) {
    double xi = x(i);
    double log_xi = std::log(xi);
    
    // Term 1: (α-1) * log(x)
    sum_term1 += (alpha - 1.0) * log_xi;
    
    // Compute x^α stably
    double x_alpha;
    if (alpha > 100.0 || (alpha * log_xi < -700.0)) {
      x_alpha = safe_exp(alpha * log_xi);
    } else {
      x_alpha = std::pow(xi, alpha);
    }
    
    // Compute v = 1 - x^α with precision for x^α near 1
    double one_minus_x_alpha;
    double log_one_minus_x_alpha;
    
    if (x_alpha > 0.9995) {
      one_minus_x_alpha = -std::expm1(alpha * log_xi);
      log_one_minus_x_alpha = safe_log(one_minus_x_alpha);
    } else {
      one_minus_x_alpha = 1.0 - x_alpha;
      log_one_minus_x_alpha = safe_log(one_minus_x_alpha);
    }
    
    // Term 2: (β-1) * log(v)
    sum_term2 += (beta - 1.0) * log_one_minus_x_alpha;
    
    // Compute v^β stably
    double v_beta;
    if (beta > 100.0 || (beta * log_one_minus_x_alpha < -700.0)) {
      v_beta = safe_exp(beta * log_one_minus_x_alpha);
    } else {
      v_beta = std::pow(one_minus_x_alpha, beta);
    }
    
    // Compute w = 1 - v^β with precision for v^β near 1
    double one_minus_v_beta;
    double log_one_minus_v_beta;
    
    if (v_beta > 0.9995) {
      one_minus_v_beta = -std::expm1(beta * log_one_minus_x_alpha);
    } else {
      one_minus_v_beta = 1.0 - v_beta;
    }
    
    // Prevent extreme underflow
    if (one_minus_v_beta < 1e-300) {
      one_minus_v_beta = 1e-300;
    }
    
    log_one_minus_v_beta = safe_log(one_minus_v_beta);
    
    // Term 3: (λ-1) * log(w); skip when λ is at machine-epsilon distance from 1
    if (std::abs(lambda - 1.0) > 1e-15) {
      sum_term3 += (lambda - 1.0) * log_one_minus_v_beta;
    }
  }
  
  // Combine all terms
  double loglike = const_term + sum_term1 + sum_term2 + sum_term3;
  
  return -loglike;
}


// ============================================================================
// GRADIENT OF NEGATIVE LOG-LIKELIHOOD
// ============================================================================

/**
 * @brief Gradient of Negative Log-Likelihood for EKw Distribution
 * 
 * Computes the gradient vector of the negative log-likelihood for
 * optimization-based parameter estimation.
 * 
 * @param par Parameter vector of length 3: (α, β, λ)
 * @param data Vector of observations (must be in (0,1))
 * 
 * @return NumericVector of length 3 containing partial derivatives
 *         with respect to (α, β, λ)
 * 
 * @details
 * The gradient components are:
 * - ∂ℓ/∂α = n/α + Σlog(x) - Σ[x^α log(x) * ((β-1)/v - (λ-1)βv^(β-1)/w)]
 * - ∂ℓ/∂β = n/β + Σlog(v) - (λ-1)Σ[v^β log(v)/w]
 * - ∂ℓ/∂λ = n/λ + Σlog(w)
 * 
 * @note Exported as .grekw_cpp for internal package use
 */
// [[Rcpp::export(.grekw_cpp)]]
Rcpp::NumericVector grekw(const Rcpp::NumericVector& par, const Rcpp::NumericVector& data) {
  // Validate parameter vector length
  if (par.size() < 3) {
    return Rcpp::NumericVector(3, R_NaN);
  }
  
  // Extract parameters
  double alpha = par[0];
  double beta = par[1];
  double lambda = par[2];
  
  // Validate parameters using consistent checker
  if (!check_ekw_pars(alpha, beta, lambda)) {
    return Rcpp::NumericVector(3, R_NaN);
  }
  
  // Convert and validate data
  arma::vec x = Rcpp::as<arma::vec>(data);
  if (x.n_elem < 1 || arma::any(x <= 0) || arma::any(x >= 1)) {
    return Rcpp::NumericVector(3, R_NaN);
  }
  
  int n = x.n_elem;
  Rcpp::NumericVector grad(3, 0.0);
  
  // Numerical stability constants
  const double min_val = 1e-10;
  
  // Initialize gradient accumulators
  double d_alpha = n / alpha;
  double d_beta = n / beta;
  double d_lambda = n / lambda;
  
  for (int i = 0; i < n; i++) {
    double xi = x(i);
    double log_xi = std::log(xi);
    d_alpha += log_xi;
    
    // Compute x^α stably
    double x_alpha;
    if (alpha > 100.0 || (alpha * log_xi < -700.0)) {
      x_alpha = safe_exp(alpha * log_xi);
    } else {
      x_alpha = std::pow(xi, alpha);
    }
    
    // Compute v = 1 - x^α with precision
    double v;
    if (x_alpha > 0.9995) {
      v = -std::expm1(alpha * log_xi);
    } else {
      v = 1.0 - x_alpha;
    }
    v = std::max(v, min_val);
    double log_v = safe_log(v);
    d_beta += log_v;
    
    // Compute v^β and v^(β-1) stably
    double v_beta, v_beta_m1;
    if (beta > 100.0 || (beta * log_v < -700.0)) {
      double log_v_beta = beta * log_v;
      v_beta = safe_exp(log_v_beta);
      v_beta_m1 = safe_exp((beta - 1.0) * log_v);
    } else {
      v_beta = std::pow(v, beta);
      v_beta_m1 = std::pow(v, beta - 1.0);
    }
    
    // Compute w = 1 - v^β with precision
    double w;
    if (v_beta > 0.9995) {
      w = -std::expm1(beta * log_v);
    } else {
      w = 1.0 - v_beta;
    }
    w = std::max(w, min_val);
    double log_w = safe_log(w);
    d_lambda += log_w;
    
    // ---- Alpha gradient component ----
    double x_alpha_log_x = x_alpha * log_xi;
    
    // Calculate (β-1)/v term
    double alpha_term1 = 0.0;
    if (std::abs(beta - 1.0) > 1e-14) {
      alpha_term1 = (beta - 1.0) / v;
    }
    
    // Calculate (λ-1) * β * v^(β-1) / w term
    double alpha_term2 = 0.0;
    if (std::abs(lambda - 1.0) > 1e-14) {
      alpha_term2 = (lambda - 1.0) * beta * v_beta_m1 / w;
    }

    d_alpha -= x_alpha_log_x * (alpha_term1 - alpha_term2);

    // ---- Beta gradient component ----
    double beta_term = 0.0;
    if (std::abs(lambda - 1.0) > 1e-14) {
      beta_term = v_beta * log_v * (lambda - 1.0) / w;
    }

    d_beta -= beta_term;
  }
  
  // Return NEGATIVE gradient (for minimization)
  grad[0] = -d_alpha;
  grad[1] = -d_beta;
  grad[2] = -d_lambda;
  
  return grad;
}


// ============================================================================
// HESSIAN OF NEGATIVE LOG-LIKELIHOOD
// ============================================================================

/**
 * @brief Hessian Matrix of Negative Log-Likelihood for EKw Distribution
 * 
 * Computes the Hessian matrix (matrix of second partial derivatives) of
 * the negative log-likelihood for standard error estimation and
 * optimization algorithms.
 * 
 * @param par Parameter vector of length 3: (α, β, λ)
 * @param data Vector of observations (must be in (0,1))
 * 
 * @return NumericMatrix of dimension 3×3 containing the Hessian
 * 
 * @details
 * Computes analytical second derivatives. The Hessian is symmetric.
 * Parameter ordering: (α, β, λ) → indices (0, 1, 2).
 * 
 * Returns NaN matrix for invalid inputs.
 * 
 * @note Exported as .hsekw_cpp for internal package use
 */
// [[Rcpp::export(.hsekw_cpp)]]
Rcpp::NumericMatrix hsekw(const Rcpp::NumericVector& par, const Rcpp::NumericVector& data) {
  // Initialize NaN matrix for error cases
  Rcpp::NumericMatrix nanH(3, 3);
  nanH.fill(R_NaN);
  
  // Validate parameter vector length
  if (par.size() < 3) {
    return nanH;
  }
  
  // Extract parameters
  double alpha = par[0];
  double beta = par[1];
  double lambda = par[2];
  
  // Validate parameters using consistent checker
  if (!check_ekw_pars(alpha, beta, lambda)) {
    return nanH;
  }
  
  // Convert and validate data
  arma::vec x = Rcpp::as<arma::vec>(data);
  if (x.n_elem < 1 || arma::any(x <= 0) || arma::any(x >= 1)) {
    return nanH;
  }
  
  int n = x.n_elem;
  
  // Initialize Hessian matrix
  arma::mat H(3, 3, arma::fill::zeros);
  
  // Numerical stability constants
  const double min_v = 1e-10;
  const double min_w = 1e-10;
  
  // Constant diagonal terms
  H(0, 0) = -n / (alpha * alpha);   // -n/α²
  H(1, 1) = -n / (beta * beta);     // -n/β²
  H(2, 2) = -n / (lambda * lambda); // -n/λ²
  
  // Special handling for λ ≈ 1
  bool lambda_near_one = std::abs(lambda - 1.0) < 1e-8;
  
  for (int i = 0; i < n; i++) {
    double xi = x(i);
    double log_xi = safe_log(xi);
    
    // ---- Compute A = x^α and derivatives ----
    double A, dA_dalpha, d2A_dalpha2;
    if (alpha > 100.0 || (alpha * log_xi < -700.0)) {
      double log_A = alpha * log_xi;
      A = safe_exp(log_A);
      dA_dalpha = A * log_xi;
      d2A_dalpha2 = A * log_xi * log_xi;
    } else {
      A = std::pow(xi, alpha);
      dA_dalpha = A * log_xi;
      d2A_dalpha2 = A * log_xi * log_xi;
    }
    
    // ---- Compute v = 1 - A and derivatives ----
    double v;
    if (A > 0.9995) {
      v = -std::expm1(alpha * log_xi);
    } else {
      v = 1.0 - A;
    }
    v = std::max(v, min_v);
    double log_v = safe_log(v);
    
    double dv_dalpha = -dA_dalpha;
    double d2v_dalpha2 = -d2A_dalpha2;
    
    // ---- Derivatives for L5: (β-1)*log(v) ----
    double d2L5_dalpha2 = 0.0;
    double d2L5_dalpha_dbeta = 0.0;
    
    if (beta != 1.0) {
      double v_squared = std::max(v * v, 1e-200);
      d2L5_dalpha2 = (beta - 1.0) * ((d2v_dalpha2 * v - dv_dalpha * dv_dalpha) / v_squared);
      d2L5_dalpha_dbeta = dv_dalpha / v;
    }
    
    // ---- Compute w = 1 - v^β and derivatives ----
    double v_beta, v_beta_m1, v_beta_m2;
    if (beta > 100.0 || (beta * log_v < -700.0)) {
      v_beta = safe_exp(beta * log_v);
      v_beta_m1 = safe_exp((beta - 1.0) * log_v);
      v_beta_m2 = safe_exp((beta - 2.0) * log_v);
    } else {
      v_beta = std::pow(v, beta);
      v_beta_m1 = std::pow(v, beta - 1.0);
      v_beta_m2 = std::pow(v, beta - 2.0);
    }
    
    double w;
    if (v_beta > 0.9995) {
      w = -std::expm1(beta * log_v);
    } else {
      w = 1.0 - v_beta;
    }
    w = std::max(w, min_w);
    double w_squared = std::max(w * w, 1e-200);
    
    // First derivatives of w
    double dw_dv = -beta * v_beta_m1;
    double dw_dalpha = dw_dv * dv_dalpha;
    double dw_dbeta = -v_beta * log_v;
    
    // Second derivatives of w
    double d2w_dalpha2 = -beta * ((beta - 1.0) * v_beta_m2 * (dv_dalpha * dv_dalpha) +
                                  v_beta_m1 * d2v_dalpha2);
    double d2w_dbeta2 = -v_beta * (log_v * log_v);
    double d_dw_dalpha_dbeta = -v_beta_m1 * (1.0 + beta * log_v) * dv_dalpha;
    
    // ---- Derivatives for L6: (λ-1)*log(w) ----
    double d2L6_dalpha2 = 0.0;
    double d2L6_dbeta2 = 0.0;
    double d2L6_dalpha_dbeta = 0.0;
    double d2L6_dalpha_dlambda = 0.0;
    double d2L6_dbeta_dlambda = 0.0;
    
    if (lambda_near_one) {
      // For λ ≈ 1, handle carefully
      if (std::abs(lambda - 1.0) > 1e-15) {
        double factor = lambda - 1.0;
        d2L6_dalpha2 = factor * ((d2w_dalpha2 * w - (dw_dalpha * dw_dalpha)) / w_squared);
        d2L6_dbeta2 = factor * ((d2w_dbeta2 * w - (dw_dbeta * dw_dbeta)) / w_squared);
        d2L6_dalpha_dbeta = factor * ((d_dw_dalpha_dbeta / w) - (dw_dalpha * dw_dbeta) / w_squared);
      }
      d2L6_dalpha_dlambda = dw_dalpha / w;
      d2L6_dbeta_dlambda = dw_dbeta / w;
    } else {
      // Standard case
      d2L6_dalpha2 = (lambda - 1.0) * ((d2w_dalpha2 * w - (dw_dalpha * dw_dalpha)) / w_squared);
      d2L6_dbeta2 = (lambda - 1.0) * ((d2w_dbeta2 * w - (dw_dbeta * dw_dbeta)) / w_squared);
      d2L6_dalpha_dbeta = (lambda - 1.0) * ((d_dw_dalpha_dbeta / w) - (dw_dalpha * dw_dbeta) / w_squared);
      d2L6_dalpha_dlambda = dw_dalpha / w;
      d2L6_dbeta_dlambda = dw_dbeta / w;
    }
    
    // ---- Accumulate Hessian contributions (upper triangle only) ----
    H(0, 0) += d2L5_dalpha2 + d2L6_dalpha2;
    H(0, 1) += d2L5_dalpha_dbeta + d2L6_dalpha_dbeta;
    H(1, 1) += d2L6_dbeta2;
    H(0, 2) += d2L6_dalpha_dlambda;
    H(1, 2) += d2L6_dbeta_dlambda;
  }

  // Symmetrize once after accumulation
  H = arma::symmatu(H);
  
  // Return NEGATIVE Hessian (for minimization)
  return Rcpp::wrap(-H);
}

