/**
 * @file gkw.cpp
 * @brief Generalized Kumaraswamy (GKw) Distribution Functions
 * 
 * @details
 * This file implements the full suite of distribution functions for the
 * five-parameter Generalized Kumaraswamy (GKw) distribution, which serves
 * as the parent distribution for six sub-families in the gkwdist package.
 * 
 * The GKw distribution has probability density function:
 * \deqn{
 *   f(x; \alpha, \beta, \gamma, \delta, \lambda) = 
 *   \frac{\lambda \alpha \beta}{B(\gamma, \delta+1)} x^{\alpha-1} (1-x^\alpha)^{\beta-1}
 *   [1-(1-x^\alpha)^\beta]^{\gamma\lambda-1} \{1-[1-(1-x^\alpha)^\beta]^\lambda\}^\delta
 * }
 * for \eqn{x \in (0,1)}, where \eqn{B(\cdot,\cdot)} is the Beta function.
 * 
 * **Parameter Constraints:**
 * - \eqn{\alpha > 0} (shape parameter)
 * - \eqn{\beta > 0} (shape parameter)
 * - \eqn{\gamma > 0} (shape parameter)
 * - \eqn{\delta \geq 0} (shape parameter)
 * - \eqn{\lambda > 0} (shape parameter)
 * 
 * **Special Cases (Sub-families):**
 * | Distribution | Parameters | Relation to GKw |
 * |--------------|------------|-----------------|
 * | Beta-Kumaraswamy (BKw) | \eqn{\alpha, \beta, \gamma, \delta} | \eqn{\lambda = 1} |
 * | Kumaraswamy-Kumaraswamy (KKw) | \eqn{\alpha, \beta, \delta, \lambda} | \eqn{\gamma = 1} |
 * | Exponentiated Kumaraswamy (EKw) | \eqn{\alpha, \beta, \lambda} | \eqn{\gamma = 1, \delta = 0} |
 * | McDonald/Beta-Power (Mc/BP) | \eqn{\gamma, \delta, \lambda} | \eqn{\alpha = \beta = 1} |
 * | Kumaraswamy (Kw) | \eqn{\alpha, \beta} | \eqn{\gamma = 1, \delta = 0, \lambda = 1} |
 * | Beta | \eqn{\gamma, \delta} | \eqn{\alpha = \beta = \lambda = 1} |
 * 
 * **Numerical Stability:**
 * All computations use log-space arithmetic and numerically stable helper
 * functions from utils.h to prevent overflow/underflow across the full
 * parameter space. Key techniques include:
 * - log1mexp() for computing log(1 - exp(x)) stably
 * - safe_pow(), safe_exp(), safe_log() for protected arithmetic
 * - Vectorized operations via Armadillo with element-wise stability
 * 
 * **Implemented Functions:**
 * - dgkw(): Probability density function (PDF)
 * - pgkw(): Cumulative distribution function (CDF)
 * - qgkw(): Quantile function (inverse CDF)
 * - rgkw(): Random variate generation
 * - llgkw(): Negative log-likelihood for MLE
 * - grgkw(): Gradient of negative log-likelihood
 * - hsgkw(): Hessian of negative log-likelihood
 * 
 * @author Lopes, J. E.
 * @date 2025-01-07
 * 
 * @references
 * Cordeiro, G. M., & de Castro, M. (2011). A new family of generalized
 * distributions. Journal of Statistical Computation and Simulation.
 * 
 * Kumaraswamy, P. (1980). A generalized probability density function for
 * double-bounded random processes. Journal of Hydrology, 46(1-2), 79-88.
 * 
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
 * @brief Probability Density Function (PDF) of the GKw Distribution
 * 
 * Computes the density (or log-density) for the Generalized Kumaraswamy
 * distribution at specified quantiles.
 * 
 * @param x Vector of quantiles (values in (0,1))
 * @param alpha Shape parameter vector (\eqn{\alpha > 0})
 * @param beta Shape parameter vector (\eqn{\beta > 0})
 * @param gamma Shape parameter vector (\eqn{\gamma > 0})
 * @param delta Shape parameter vector (\eqn{\delta \geq 0})
 * @param lambda Shape parameter vector (\eqn{\lambda > 0})
 * @param log_prob If TRUE, returns log-density; otherwise returns density
 * 
 * @return NumericVector of density values (or log-density if log_prob=TRUE)
 * 
 * @details
 * For x outside (0,1), returns 0 (or -Inf for log-density).
 * Uses log-space computation throughout to ensure numerical stability.
 * Parameters are recycled to match the longest input vector.
 * 
 * @note Exported as .dgkw_cpp for internal package use
 */
// [[Rcpp::export(.dgkw_cpp)]]
Rcpp::NumericVector dgkw(
    const arma::vec& x,
    const Rcpp::NumericVector& alpha,
    const Rcpp::NumericVector& beta,
    const Rcpp::NumericVector& gamma,
    const Rcpp::NumericVector& delta,
    const Rcpp::NumericVector& lambda,
    bool log_prob = false
) {
  // Convert R vectors to Armadillo vectors for efficient computation
  arma::vec alpha_vec(alpha.begin(), alpha.size());
  arma::vec beta_vec(beta.begin(), beta.size());
  arma::vec gamma_vec(gamma.begin(), gamma.size());
  arma::vec delta_vec(delta.begin(), delta.size());
  arma::vec lambda_vec(lambda.begin(), lambda.size());
  
  // Determine output length (maximum of all input lengths for recycling)
  size_t n = std::max({x.n_elem, alpha_vec.n_elem, beta_vec.n_elem,
                      gamma_vec.n_elem, delta_vec.n_elem, lambda_vec.n_elem});
  
  // Initialize result vector with appropriate default values
  arma::vec result(n);
  result.fill(log_prob ? R_NegInf : 0.0);
  
  // Process each element with parameter recycling
  for (size_t i = 0; i < n; ++i) {
    // Extract recycled parameter values
    double a = alpha_vec[i % alpha_vec.n_elem];
    double b = beta_vec[i % beta_vec.n_elem];
    double g = gamma_vec[i % gamma_vec.n_elem];
    double d = delta_vec[i % delta_vec.n_elem];
    double l = lambda_vec[i % lambda_vec.n_elem];
    double xi = x[i % x.n_elem];
    
    // Validate parameters
    if (!check_pars(a, b, g, d, l)) {
      Rcpp::warning("dgkw: invalid parameters at index %d (alpha,beta,gamma>0, delta>=0, lambda>0)", i+1);
      continue;
    }
    
    // Check support: x must be in (0, 1)
    if (xi <= 0.0 || xi >= 1.0 || !R_finite(xi)) {
      continue;
    }
    
    // ---- Log-space computation of density ----
    
    // Normalization constant: log(λαβ / B(γ, δ+1))
    double log_beta_val = R::lbeta(g, d + 1.0);
    double log_const = std::log(l) + std::log(a) + std::log(b) - log_beta_val;
    double gamma_lambda = g * l;
    
    // Compute x^α in log-space
    double x_alpha = safe_pow(xi, a);
    if (!R_finite(x_alpha) || x_alpha >= 1.0 - SQRT_EPSILON) {
      continue;
    }
    double log_x_alpha = safe_log(x_alpha);
    
    // Compute log(1 - x^α) using stable log1mexp
    double log_one_minus_x_alpha = log1mexp(log_x_alpha);
    if (!R_finite(log_one_minus_x_alpha)) {
      continue;
    }
    
    // Compute log((1 - x^α)^β) = β * log(1 - x^α)
    double log_one_minus_x_alpha_beta = b * log_one_minus_x_alpha;
    if (!R_finite(log_one_minus_x_alpha_beta)) {
      continue;
    }
    
    // Compute log(1 - (1 - x^α)^β) = log(w)
    double log_term1 = log1mexp(log_one_minus_x_alpha_beta);
    if (!R_finite(log_term1)) {
      continue;
    }
    
    // Compute log([1-(1-x^α)^β]^λ) = λ * log(w)
    double log_term1_lambda = l * log_term1;
    if (!R_finite(log_term1_lambda)) {
      continue;
    }
    
    // Compute log(1 - [1-(1-x^α)^β]^λ) = log(z)
    double log_term2 = log1mexp(log_term1_lambda);
    if (!R_finite(log_term2)) {
      continue;
    }
    
    // Assemble log-density:
    // log(f) = log_const + (α-1)*log(x) + (β-1)*log(v) + (γλ-1)*log(w) + δ*log(z)
    double logdens = log_const +
      (a - 1.0) * std::log(xi) +
      (b - 1.0) * log_one_minus_x_alpha +
      (gamma_lambda - 1.0) * log_term1 +
      d * log_term2;
    
    // Validate final result
    if (!R_finite(logdens)) {
      continue;
    }
    
    // Return appropriate scale
    result(i) = log_prob ? logdens : safe_exp(logdens);
  }
  
  return Rcpp::NumericVector(result.memptr(), result.memptr() + result.n_elem);
}


// ============================================================================
// CUMULATIVE DISTRIBUTION FUNCTION
// ============================================================================

/**
 * @brief Cumulative Distribution Function (CDF) of the GKw Distribution
 * 
 * Computes the cumulative probability for the Generalized Kumaraswamy
 * distribution at specified quantiles.
 * 
 * @param q Vector of quantiles
 * @param alpha Shape parameter vector (\eqn{\alpha > 0})
 * @param beta Shape parameter vector (\eqn{\beta > 0})
 * @param gamma Shape parameter vector (\eqn{\gamma > 0})
 * @param delta Shape parameter vector (\eqn{\delta \geq 0})
 * @param lambda Shape parameter vector (\eqn{\lambda > 0})
 * @param lower_tail If TRUE, returns P(X <= q); otherwise P(X > q)
 * @param log_p If TRUE, returns log-probability
 * 
 * @return NumericVector of cumulative probabilities
 * 
 * @details
 * The CDF is computed as:
 * \deqn{F(x) = I_{[1-(1-x^\alpha)^\beta]^\lambda}(\gamma, \delta+1)}
 * where \eqn{I_y(a,b)} is the regularized incomplete Beta function.
 * 
 * @note Exported as .pgkw_cpp for internal package use
 */
// [[Rcpp::export(.pgkw_cpp)]]
Rcpp::NumericVector pgkw(
    const arma::vec& q,
    const Rcpp::NumericVector& alpha,
    const Rcpp::NumericVector& beta,
    const Rcpp::NumericVector& gamma,
    const Rcpp::NumericVector& delta,
    const Rcpp::NumericVector& lambda,
    bool lower_tail = true,
    bool log_p = false
) {
  // Convert R vectors to Armadillo vectors
  arma::vec alpha_vec(alpha.begin(), alpha.size());
  arma::vec beta_vec(beta.begin(), beta.size());
  arma::vec gamma_vec(gamma.begin(), gamma.size());
  arma::vec delta_vec(delta.begin(), delta.size());
  arma::vec lambda_vec(lambda.begin(), lambda.size());
  
  // Determine output length for recycling
  size_t n = std::max({q.n_elem, alpha_vec.n_elem, beta_vec.n_elem,
                      gamma_vec.n_elem, delta_vec.n_elem, lambda_vec.n_elem});
  
  arma::vec result(n);
  
  for (size_t i = 0; i < n; ++i) {
    // Extract recycled parameter values
    double a = alpha_vec[i % alpha_vec.n_elem];
    double b = beta_vec[i % beta_vec.n_elem];
    double g = gamma_vec[i % gamma_vec.n_elem];
    double d = delta_vec[i % delta_vec.n_elem];
    double l = lambda_vec[i % lambda_vec.n_elem];
    double qi = q[i % q.n_elem];
    
    // Validate parameters
    if (!check_pars(a, b, g, d, l)) {
      result(i) = NA_REAL;
      Rcpp::warning("pgkw: invalid parameters at index %d (alpha,beta,gamma>0, delta>=0, lambda>0)", i+1);
      continue;
    }
    
    // Handle boundary cases: q <= 0
    if (!R_finite(qi) || qi <= 0.0) {
      result(i) = lower_tail ? (log_p ? R_NegInf : 0.0) : (log_p ? 0.0 : 1.0);
      continue;
    }
    
    // Handle boundary cases: q >= 1
    if (qi >= 1.0) {
      result(i) = lower_tail ? (log_p ? 0.0 : 1.0) : (log_p ? R_NegInf : 0.0);
      continue;
    }
    
    // ---- Compute CDF via Beta function ----
    
    // Step 1: q^α
    double qi_alpha = safe_pow(qi, a);
    if (!R_finite(qi_alpha)) {
      result(i) = lower_tail ? (log_p ? R_NegInf : 0.0) : (log_p ? 0.0 : 1.0);
      continue;
    }
    
    // Step 2: log(1 - q^α)
    double log_qi_alpha = safe_log(qi_alpha);
    double log_one_minus_qi_alpha = log1mexp(log_qi_alpha);
    if (!R_finite(log_one_minus_qi_alpha)) {
      result(i) = lower_tail ? (log_p ? 0.0 : 1.0) : (log_p ? R_NegInf : 0.0);
      continue;
    }
    
    // Step 3: (1 - q^α)^β - use log_one_minus_qi_alpha directly (avoid redundant log)
    double log_oma_beta = b * log_one_minus_qi_alpha;
    if (!R_finite(log_oma_beta)) {
      result(i) = lower_tail ? (log_p ? R_NegInf : 0.0) : (log_p ? 0.0 : 1.0);
      continue;
    }
    double oma_beta = safe_exp(log_oma_beta);
    
    // Step 4: 1 - (1 - q^α)^β
    double term = 1.0 - oma_beta;
    if (term <= 0.0) {
      result(i) = lower_tail ? (log_p ? R_NegInf : 0.0) : (log_p ? 0.0 : 1.0);
      continue;
    }
    if (term >= 1.0) {
      result(i) = lower_tail ? (log_p ? 0.0 : 1.0) : (log_p ? R_NegInf : 0.0);
      continue;
    }
    
    // Step 5: [1 - (1 - q^α)^β]^λ
    double log_term = safe_log(term);
    double log_y = l * log_term;
    if (!R_finite(log_y)) {
      result(i) = lower_tail ? (log_p ? R_NegInf : 0.0) : (log_p ? 0.0 : 1.0);
      continue;
    }
    double y = safe_exp(log_y);
    
    // Boundary validation for y
    if (y <= 0.0) {
      result(i) = lower_tail ? (log_p ? R_NegInf : 0.0) : (log_p ? 0.0 : 1.0);
      continue;
    }
    if (y >= 1.0) {
      result(i) = lower_tail ? (log_p ? 0.0 : 1.0) : (log_p ? R_NegInf : 0.0);
      continue;
    }
    
    // Step 6: Compute I_y(γ, δ+1) via R's pbeta
    double prob = R::pbeta(y, g, d + 1.0, true, false);
    
    // Apply tail adjustment
    if (!lower_tail) {
      prob = 1.0 - prob;
    }
    
    // Apply log transformation if requested
    if (log_p) {
      if (prob <= 0.0) {
        prob = R_NegInf;
      } else if (prob >= 1.0) {
        prob = 0.0;
      } else {
        prob = std::log(prob);
      }
    }
    
    result(i) = prob;
  }
  
  return Rcpp::NumericVector(result.memptr(), result.memptr() + result.n_elem);
}


// ============================================================================
// QUANTILE FUNCTION
// ============================================================================

/**
 * @brief Quantile Function (Inverse CDF) of the GKw Distribution
 * 
 * Computes quantiles for the Generalized Kumaraswamy distribution
 * given probability values.
 * 
 * @param p Vector of probabilities (values in [0,1])
 * @param alpha Shape parameter vector (\eqn{\alpha > 0})
 * @param beta Shape parameter vector (\eqn{\beta > 0})
 * @param gamma Shape parameter vector (\eqn{\gamma > 0})
 * @param delta Shape parameter vector (\eqn{\delta \geq 0})
 * @param lambda Shape parameter vector (\eqn{\lambda > 0})
 * @param lower_tail If TRUE, probabilities are P(X <= x); otherwise P(X > x)
 * @param log_p If TRUE, probabilities are given as log(p)
 * 
 * @return NumericVector of quantiles
 * 
 * @details
 * The quantile function is computed by inverting the CDF:
 * \deqn{Q(p) = \{1 - [1 - y^{1/\lambda}]^{1/\beta}\}^{1/\alpha}}
 * where \eqn{y = Q_{Beta}(p; \gamma, \delta+1)}.
 * 
 * @note Exported as .qgkw_cpp for internal package use
 */
// [[Rcpp::export(.qgkw_cpp)]]
Rcpp::NumericVector qgkw(
    const arma::vec& p,
    const Rcpp::NumericVector& alpha,
    const Rcpp::NumericVector& beta,
    const Rcpp::NumericVector& gamma,
    const Rcpp::NumericVector& delta,
    const Rcpp::NumericVector& lambda,
    bool lower_tail = true,
    bool log_p = false
) {
  // Convert R vectors to Armadillo vectors
  arma::vec alpha_vec(alpha.begin(), alpha.size());
  arma::vec beta_vec(beta.begin(), beta.size());
  arma::vec gamma_vec(gamma.begin(), gamma.size());
  arma::vec delta_vec(delta.begin(), delta.size());
  arma::vec lambda_vec(lambda.begin(), lambda.size());
  
  // Determine output length for recycling
  size_t n = std::max({p.n_elem, alpha_vec.n_elem, beta_vec.n_elem,
                      gamma_vec.n_elem, delta_vec.n_elem, lambda_vec.n_elem});
  
  arma::vec result(n);
  
  for (size_t i = 0; i < n; ++i) {
    // Extract recycled parameter values
    double a = alpha_vec[i % alpha_vec.n_elem];
    double b = beta_vec[i % beta_vec.n_elem];
    double g = gamma_vec[i % gamma_vec.n_elem];
    double d = delta_vec[i % delta_vec.n_elem];
    double l = lambda_vec[i % lambda_vec.n_elem];
    double pp = p[i % p.n_elem];
    
    // Validate parameters
    if (!check_pars(a, b, g, d, l)) {
      result(i) = NA_REAL;
      Rcpp::warning("qgkw: invalid parameters at index %d (alpha,beta,gamma>0, delta>=0, lambda>0)", i+1);
      continue;
    }
    
    // ---- Convert probability to linear scale ----
    // Handle log_p: convert from log-probability to probability
    if (log_p) {
      if (pp > 0.0) {
        // log(p) > 0 implies p > 1, which is invalid
        result(i) = NA_REAL;
        continue;
      }
      pp = safe_exp(pp);
      if (!R_finite(pp)) {
        // Handle extreme values: exp(-Inf) = 0, edge cases
        result(i) = (pp == 0.0) ? 0.0 : 1.0;
        continue;
      }
    }
    
    // Handle lower_tail: convert upper-tail to lower-tail probability
    // NOTE: At this point pp is ALWAYS in linear scale [0, 1]
    if (!lower_tail) {
      pp = 1.0 - pp;
    }
    
    // Validate probability bounds
    if (!R_finite(pp)) {
      result(i) = NA_REAL;
      continue;
    }
    
    if (pp <= 0.0) {
      result(i) = 0.0;
      continue;
    }
    
    if (pp >= 1.0) {
      result(i) = 1.0;
      continue;
    }
    
    // ---- Compute quantile via inverse transformations ----
    
    // Step 1: y = Q_Beta(p, γ, δ+1)
    double y = R::qbeta(pp, g, d + 1.0, true, false);
    
    if (!R_finite(y)) {
      result(i) = (y == R_PosInf) ? 1.0 : 0.0;
      continue;
    }
    
    if (y <= 0.0) {
      result(i) = 0.0;
      continue;
    }
    
    if (y >= 1.0) {
      result(i) = 1.0;
      continue;
    }
    
    // Step 2: v = y^(1/λ)
    double v = (l == 1.0) ? y : safe_pow(y, 1.0/l);
    if (!R_finite(v)) {
      result(i) = (v == R_PosInf) ? 1.0 : 0.0;
      continue;
    }
    
    // Step 3: tmp = 1 - v
    double tmp = 1.0 - v;
    if (!R_finite(tmp)) {
      result(i) = (tmp == R_PosInf) ? 0.0 : 1.0;
      continue;
    }
    
    if (tmp <= 0.0) {
      result(i) = 1.0;
      continue;
    }
    
    if (tmp >= 1.0) {
      result(i) = 0.0;
      continue;
    }
    
    // Step 4: tmp2 = tmp^(1/β)
    double tmp2 = (b == 1.0) ? tmp : safe_pow(tmp, 1.0/b);
    if (!R_finite(tmp2)) {
      result(i) = (tmp2 == R_PosInf) ? 0.0 : 1.0;
      continue;
    }
    
    if (tmp2 <= 0.0) {
      result(i) = 1.0;
      continue;
    }
    
    if (tmp2 >= 1.0) {
      result(i) = 0.0;
      continue;
    }
    
    // Step 5: q = (1 - tmp2)^(1/α)
    double one_minus_tmp2 = 1.0 - tmp2;
    if (!R_finite(one_minus_tmp2)) {
      result(i) = (one_minus_tmp2 == R_PosInf) ? 0.0 : 1.0;
      continue;
    }
    
    double qq = (a == 1.0) ? one_minus_tmp2 : safe_pow(one_minus_tmp2, 1.0/a);
    if (!R_finite(qq)) {
      result(i) = (qq == R_PosInf) ? 1.0 : 0.0;
      continue;
    }
    
    // Clamp result to valid support [0, 1]
    if (qq < 0.0) {
      qq = 0.0;
    } else if (qq > 1.0) {
      qq = 1.0;
    }
    
    result(i) = qq;
  }
  
  return Rcpp::NumericVector(result.memptr(), result.memptr() + result.n_elem);
}


// ============================================================================
// RANDOM NUMBER GENERATION
// ============================================================================

/**
 * @brief Random Variate Generation for the GKw Distribution
 * 
 * Generates random samples from the Generalized Kumaraswamy distribution
 * using the inverse transform method.
 * 
 * @param n Number of random variates to generate
 * @param alpha Shape parameter vector (\eqn{\alpha > 0})
 * @param beta Shape parameter vector (\eqn{\beta > 0})
 * @param gamma Shape parameter vector (\eqn{\gamma > 0})
 * @param delta Shape parameter vector (\eqn{\delta \geq 0})
 * @param lambda Shape parameter vector (\eqn{\lambda > 0})
 * 
 * @return NumericVector of n random variates from GKw distribution
 * 
 * @details
 * Uses the representation: if V ~ Beta(γ, δ+1), then
 * \deqn{X = \{1 - [1 - V^{1/\lambda}]^{1/\beta}\}^{1/\alpha} \sim GKw}
 * 
 * @note Exported as .rgkw_cpp for internal package use
 */
// [[Rcpp::export(.rgkw_cpp)]]
Rcpp::NumericVector rgkw(
    int n,
    const Rcpp::NumericVector& alpha,
    const Rcpp::NumericVector& beta,
    const Rcpp::NumericVector& gamma,
    const Rcpp::NumericVector& delta,
    const Rcpp::NumericVector& lambda
) {
  // Convert R vectors to Armadillo vectors
  arma::vec alpha_vec(alpha.begin(), alpha.size());
  arma::vec beta_vec(beta.begin(), beta.size());
  arma::vec gamma_vec(gamma.begin(), gamma.size());
  arma::vec delta_vec(delta.begin(), delta.size());
  arma::vec lambda_vec(lambda.begin(), lambda.size());
  
  arma::vec result(n);
  
  for (int i = 0; i < n; ++i) {
    // Extract recycled parameter values
    double a = alpha_vec[i % alpha_vec.n_elem];
    double b = beta_vec[i % beta_vec.n_elem];
    double g = gamma_vec[i % gamma_vec.n_elem];
    double d = delta_vec[i % delta_vec.n_elem];
    double l = lambda_vec[i % lambda_vec.n_elem];
    
    // Validate parameters
    if (!check_pars(a, b, g, d, l)) {
      result(i) = NA_REAL;
      Rcpp::warning("rgkw: invalid parameters at index %d (alpha,beta,gamma>0, delta>=0, lambda>0)", i+1);
      continue;
    }
    
    // Generate V ~ Beta(γ, δ+1)
    double vi = R::rbeta(g, d + 1.0);
    
    // Handle boundary cases
    if (vi <= 0.0) {
      result(i) = 0.0;
      continue;
    }
    
    if (vi >= 1.0) {
      result(i) = 1.0;
      continue;
    }
    
    // Transform: v = V^(1/λ)
    double vl = (l == 1.0) ? vi : safe_pow(vi, 1.0/l);
    if (!R_finite(vl)) {
      result(i) = (vl == R_PosInf) ? 1.0 : 0.0;
      continue;
    }
    
    // Transform: tmp = 1 - v
    double tmp = 1.0 - vl;
    if (!R_finite(tmp)) {
      result(i) = (tmp == R_PosInf) ? 0.0 : 1.0;
      continue;
    }
    
    if (tmp <= 0.0) {
      result(i) = 1.0;
      continue;
    }
    
    if (tmp >= 1.0) {
      result(i) = 0.0;
      continue;
    }
    
    // Transform: tmp2 = tmp^(1/β)
    double tmp2 = (b == 1.0) ? tmp : safe_pow(tmp, 1.0/b);
    if (!R_finite(tmp2)) {
      result(i) = (tmp2 == R_PosInf) ? 0.0 : 1.0;
      continue;
    }
    
    if (tmp2 <= 0.0) {
      result(i) = 1.0;
      continue;
    }
    
    if (tmp2 >= 1.0) {
      result(i) = 0.0;
      continue;
    }
    
    // Transform: x = (1 - tmp2)^(1/α)
    double one_minus_tmp2 = 1.0 - tmp2;
    if (!R_finite(one_minus_tmp2)) {
      result(i) = (one_minus_tmp2 == R_PosInf) ? 0.0 : 1.0;
      continue;
    }
    
    double xx = (a == 1.0) ? one_minus_tmp2 : safe_pow(one_minus_tmp2, 1.0/a);
    if (!R_finite(xx)) {
      result(i) = (xx == R_PosInf) ? 1.0 : 0.0;
      continue;
    }
    
    // Clamp to valid support
    if (xx < 0.0) {
      xx = 0.0;
    } else if (xx > 1.0) {
      xx = 1.0;
    }
    
    result(i) = xx;
  }
  
  return Rcpp::NumericVector(result.memptr(), result.memptr() + result.n_elem);
}


// ============================================================================
// NEGATIVE LOG-LIKELIHOOD FUNCTION
// ============================================================================

/**
 * @brief Negative Log-Likelihood for GKw Distribution
 * 
 * Computes the negative log-likelihood function for parameter estimation
 * via maximum likelihood.
 * 
 * @param par Parameter vector of length 5: (α, β, γ, δ, λ)
 * @param data Vector of observations (must be in (0,1))
 * 
 * @return Negative log-likelihood value (scalar)
 * 
 * @details
 * The log-likelihood for n observations is:
 * \deqn{\ell(\theta) = n[\ln\lambda + \ln\alpha + \ln\beta - \ln B(\gamma,\delta+1)]
 *       + (\alpha-1)\sum\ln x_i + (\beta-1)\sum\ln(1-x_i^\alpha)
 *       + (\gamma\lambda-1)\sum\ln w_i + \delta\sum\ln z_i}
 * where \eqn{w_i = 1-(1-x_i^\alpha)^\beta} and \eqn{z_i = 1-w_i^\lambda}.
 * 
 * Returns +Inf for invalid parameters or data outside (0,1).
 * 
 * @note Exported as .llgkw_cpp for internal package use
 */
// [[Rcpp::export(.llgkw_cpp)]]
double llgkw(const Rcpp::NumericVector& par, const Rcpp::NumericVector& data) {
  // Extract parameters
  double alpha = par[0];
  double beta = par[1];
  double gamma = par[2];
  double delta = par[3];
  double lambda = par[4];
  
  // Validate parameters
  if (!check_pars(alpha, beta, gamma, delta, lambda)) {
    return R_PosInf;
  }
  
  // Convert data to Armadillo vector
  arma::vec x = Rcpp::as<arma::vec>(data);
  
  // Validate data support
  if (arma::any(x <= 0) || arma::any(x >= 1)) {
    return R_NegInf;
  }
  
  int n = x.n_elem;
  
  // ---- Compute log-likelihood terms ----
  
  // Constant term: n * log(λαβ / B(γ, δ+1))
  double log_beta_term = R::lbeta(gamma, delta + 1);
  double constant_term = n * (std::log(lambda) + std::log(alpha) + std::log(beta) - log_beta_term);
  
  // Term 1: (α-1) * Σ log(x_i)
  arma::vec log_x = vec_safe_log(x);
  double term1 = arma::sum((alpha - 1.0) * log_x);
  
  // Compute v = 1 - x^α in log-space
  arma::vec x_alpha = vec_safe_pow(x, alpha);
  arma::vec log_x_alpha = vec_safe_log(x_alpha);
  arma::vec log_v = vec_log1mexp(log_x_alpha);
  
  // Term 2: (β-1) * Σ log(v_i)
  double term2 = arma::sum((beta - 1.0) * log_v);
  
  // Compute w = 1 - v^β in log-space
  arma::vec log_v_beta = beta * log_v;
  arma::vec log_w = vec_log1mexp(log_v_beta);
  
  // Term 3: (γλ-1) * Σ log(w_i)
  double term3 = arma::sum((gamma * lambda - 1.0) * log_w);
  
  // Compute z = 1 - w^λ in log-space
  arma::vec log_w_lambda = lambda * log_w;
  arma::vec log_z = vec_log1mexp(log_w_lambda);
  
  // Term 4: δ * Σ log(z_i)
  double term4 = arma::sum(delta * log_z);
  
  // Return negative log-likelihood
  return -(constant_term + term1 + term2 + term3 + term4);
}


// ============================================================================
// GRADIENT OF NEGATIVE LOG-LIKELIHOOD
// ============================================================================

/**
 * @brief Gradient of Negative Log-Likelihood for GKw Distribution
 * 
 * Computes the gradient vector of the negative log-likelihood for
 * optimization-based parameter estimation.
 * 
 * @param par Parameter vector of length 5: (α, β, γ, δ, λ)
 * @param data Vector of observations (must be in (0,1))
 * 
 * @return NumericVector of length 5 containing partial derivatives
 *         with respect to (α, β, γ, δ, λ)
 * 
 * @details
 * Computes analytical gradients using chain rule and log-space arithmetic
 * for numerical stability. Returns NaN vector for invalid inputs.
 * 
 * @note Exported as .grgkw_cpp for internal package use
 */
// [[Rcpp::export(.grgkw_cpp)]]
Rcpp::NumericVector grgkw(const Rcpp::NumericVector& par, const Rcpp::NumericVector& data) {
  // Extract parameters
  double alpha = par[0];
  double beta = par[1];
  double gamma = par[2];
  double delta = par[3];
  double lambda = par[4];
  
  // Validate parameters
  if (!check_pars(alpha, beta, gamma, delta, lambda)) {
    return Rcpp::NumericVector(5, R_NaN);
  }
  
  // Convert and validate data
  arma::vec x = Rcpp::as<arma::vec>(data);
  
  if (arma::any(x <= 0) || arma::any(x >= 1)) {
    return Rcpp::NumericVector(5, R_NaN);
  }
  
  int n = x.n_elem;
  Rcpp::NumericVector grad(5, 0.0);
  
  // ---- Compute intermediate quantities in log-space ----
  
  arma::vec log_x = vec_safe_log(x);
  arma::vec x_alpha = vec_safe_pow(x, alpha);
  arma::vec log_x_alpha = vec_safe_log(x_alpha);
  
  // v = 1 - x^α (computed in log-space, but we only need log_v)
  arma::vec log_v = vec_log1mexp(log_x_alpha);
  
  // v^β and v^(β-1)
  arma::vec log_v_beta = beta * log_v;
  arma::vec v_beta = vec_safe_exp(log_v_beta);
  arma::vec log_v_beta_m1 = (beta - 1.0) * log_v;
  arma::vec v_beta_m1 = vec_safe_exp(log_v_beta_m1);
  
  // w = 1 - v^β (computed in log-space)
  arma::vec log_w = vec_log1mexp(log_v_beta);
  
  // w^λ and w^(λ-1)
  arma::vec log_w_lambda = lambda * log_w;
  arma::vec w_lambda = vec_safe_exp(log_w_lambda);
  arma::vec log_w_lambda_m1 = (lambda - 1.0) * log_w;
  arma::vec w_lambda_m1 = vec_safe_exp(log_w_lambda_m1);
  
  // z = 1 - w^λ (computed in log-space)
  arma::vec log_z = vec_log1mexp(log_w_lambda);
  
  // Validate intermediate calculations
  if (!log_v.is_finite() || !log_w.is_finite() || !log_z.is_finite()) {
    return Rcpp::NumericVector(5, R_NaN);
  }
  
  // ---- Compute gradient components ----
  
  // ∂ℓ/∂α = n/α + Σ log(x_i) - complex_term
  double d_alpha = n / alpha + arma::sum(log_x);
  
  // Complex term for α gradient
  arma::vec x_alpha_log_x = x_alpha % log_x;  // x^α * log(x)
  
  // Term 1: (β-1) / v
  arma::vec alpha_term1 = (beta - 1.0) * vec_safe_exp(-log_v);
  
  // Term 2: (γλ-1) * β * v^(β-1) / w
  double coeff2 = (gamma * lambda - 1.0) * beta;
  arma::vec alpha_term2 = coeff2 * v_beta_m1 % vec_safe_exp(-log_w);
  
  // Term 3: δ * λ * β * v^(β-1) * w^(λ-1) / z
  double coeff3 = delta * lambda * beta;
  arma::vec alpha_term3 = coeff3 * v_beta_m1 % w_lambda_m1 % vec_safe_exp(-log_z);
  
  d_alpha -= arma::sum(x_alpha_log_x % (alpha_term1 - alpha_term2 + alpha_term3));
  
  // ∂ℓ/∂β = n/β + Σ log(v) - complex_term
  double d_beta = n / beta + arma::sum(log_v);
  
  arma::vec v_beta_log_v = v_beta % log_v;  // v^β * log(v)
  
  // Term 1: (γλ-1) / w
  double coeff_b1 = gamma * lambda - 1.0;
  arma::vec beta_term1 = coeff_b1 * vec_safe_exp(-log_w);
  
  // Term 2: δ * λ * w^(λ-1) / z
  double coeff_b2 = delta * lambda;
  arma::vec beta_term2 = coeff_b2 * w_lambda_m1 % vec_safe_exp(-log_z);
  
  d_beta -= arma::sum(v_beta_log_v % (beta_term1 - beta_term2));
  
  // ∂ℓ/∂γ = -n[ψ(γ) - ψ(γ+δ+1)] + λ Σ log(w)
  double d_gamma = -n * (R::digamma(gamma) - R::digamma(gamma + delta + 1.0)) + 
    lambda * arma::sum(log_w);
  
  // ∂ℓ/∂δ = -n[ψ(δ+1) - ψ(γ+δ+1)] + Σ log(z)
  double d_delta = -n * (R::digamma(delta + 1.0) - R::digamma(gamma + delta + 1.0)) + 
    arma::sum(log_z);
  
  // ∂ℓ/∂λ = n/λ + γ Σ log(w) - δ Σ [(w^λ * log(w)) / z]
  double d_lambda = n / lambda + gamma * arma::sum(log_w);
  
  if (delta > 0.0) {
    arma::vec w_lambda_log_w = w_lambda % log_w;
    d_lambda -= delta * arma::sum(w_lambda_log_w % vec_safe_exp(-log_z));
  }
  
  // Validate gradient components
  if (!R_finite(d_alpha) || !R_finite(d_beta) || !R_finite(d_gamma) || 
      !R_finite(d_delta) || !R_finite(d_lambda)) {
      return Rcpp::NumericVector(5, R_NaN);
  }
  
  // Return NEGATIVE gradient (for minimization of negative log-likelihood)
  grad[0] = -d_alpha;
  grad[1] = -d_beta;
  grad[2] = -d_gamma;
  grad[3] = -d_delta;
  grad[4] = -d_lambda;
  
  return grad;
}


// ============================================================================
// HESSIAN OF NEGATIVE LOG-LIKELIHOOD
// ============================================================================

/**
 * @brief Hessian Matrix of Negative Log-Likelihood for GKw Distribution
 * 
 * Computes the Hessian matrix (matrix of second partial derivatives) of
 * the negative log-likelihood for standard error estimation and
 * optimization algorithms.
 * 
 * @param par Parameter vector of length 5: (α, β, γ, δ, λ)
 * @param data Vector of observations (must be in (0,1))
 * 
 * @return NumericMatrix of dimension 5×5 containing the Hessian
 * 
 * @details
 * Computes analytical second derivatives. The Hessian is symmetric,
 * so only unique elements are computed. Returns NaN matrix for invalid inputs.
 * 
 * Parameter ordering in matrix: (α, β, γ, δ, λ) corresponding to
 * indices (0, 1, 2, 3, 4).
 * 
 * @note Exported as .hsgkw_cpp for internal package use
 */
// [[Rcpp::export(.hsgkw_cpp)]]
Rcpp::NumericMatrix hsgkw(const Rcpp::NumericVector& par, const Rcpp::NumericVector& data) {
  // Extract parameters
  double alpha  = par[0];
  double beta   = par[1];
  double gamma  = par[2];
  double delta  = par[3];
  double lambda = par[4];
  
  // Validate parameters
  if (!check_pars(alpha, beta, gamma, delta, lambda)) {
    Rcpp::NumericMatrix nanH(5, 5);
    nanH.fill(R_NaN);
    return nanH;
  }
  
  // Convert and validate data
  arma::vec x = Rcpp::as<arma::vec>(data);
  if (arma::any(x <= 0) || arma::any(x >= 1)) {
    Rcpp::NumericMatrix nanH(5, 5);
    nanH.fill(R_NaN);
    return nanH;
  }
  
  int n = x.n_elem;
  arma::mat H(5, 5, arma::fill::zeros);
  
  // ---- Constant terms (independent of observations) ----
  
  // H(λ,λ) from n*ln(λ): ∂²/∂λ² = -n/λ²
  H(4, 4) += -n / (lambda * lambda);
  
  // H(α,α) from n*ln(α): ∂²/∂α² = -n/α²
  H(0, 0) += -n / (alpha * alpha);
  
  // H(β,β) from n*ln(β): ∂²/∂β² = -n/β²
  H(1, 1) += -n / (beta * beta);
  
  // H(γ,γ) from -n*ln(B(γ,δ+1)): involves trigamma
  H(2, 2) += -n * (R::trigamma(gamma) - R::trigamma(gamma + delta + 1));
  
  // H(δ,δ) from -n*ln(B(γ,δ+1)): involves trigamma
  H(3, 3) += -n * (R::trigamma(delta + 1) - R::trigamma(gamma + delta + 1));
  
  // H(γ,δ) = H(δ,γ): mixed derivative
  H(2, 3) += n * R::trigamma(gamma + delta + 1);
  H(3, 2) = H(2, 3);
  
  // Accumulators for mixed derivatives involving λ
  double acc_gamma_lambda = 0.0;
  double acc_delta_lambda = 0.0;
  double acc_alpha_lambda = 0.0;
  double acc_beta_lambda = 0.0;
  
  // ---- Observation-dependent terms ----
  for (int i = 0; i < n; i++) {
    double xi = x(i);
    
    // Compute A = x^α and derivatives
    double ln_xi = safe_log(xi);
    double A = safe_pow(xi, alpha);
    double dA_dalpha = A * ln_xi;
    double d2A_dalpha2 = A * ln_xi * ln_xi;
    
    // v = 1 - A and derivatives (using log-space for v)
    double log_A = alpha * ln_xi;
    double log_v = log1mexp(log_A);
    if (!R_finite(log_v)) continue;
    double v = safe_exp(log_v);
    double ln_v = log_v;
    double dv_dalpha = -dA_dalpha;
    double d2v_dalpha2 = -d2A_dalpha2;
    
    // --- L6: (β-1) ln(v) contributions ---
    double d2L6_dalpha2 = (beta - 1.0) * ((d2v_dalpha2 * v - dv_dalpha * dv_dalpha) / (v * v));
    double d2L6_dalpha_dbeta = dv_dalpha / v;
    
    // --- L7: (γλ - 1) ln(w), where w = 1 - v^β ---
    double log_v_beta = beta * log_v;
    double log_w = log1mexp(log_v_beta);
    if (!R_finite(log_w)) continue;
    double w = safe_exp(log_w);
    double ln_w = log_w;
    
    // Derivatives of w
    double v_beta_m1 = safe_pow(v, beta - 1.0);
    double dw_dv = -beta * v_beta_m1;
    double dw_dalpha = dw_dv * dv_dalpha;
    
    double d2w_dalpha2 = -beta * ((beta - 1.0) * safe_pow(v, beta - 2.0) * (dv_dalpha * dv_dalpha)
                                    + v_beta_m1 * d2v_dalpha2);
    double d2L7_dalpha2 = (gamma * lambda - 1.0) * ((d2w_dalpha2 * w - (dw_dalpha * dw_dalpha)) / (w * w));
    
    double dw_dbeta = -safe_pow(v, beta) * ln_v;
    double d2w_dbeta2 = -safe_pow(v, beta) * (ln_v * ln_v);
    double d2L7_dbeta2 = (gamma * lambda - 1.0) * ((d2w_dbeta2 * w - (dw_dbeta * dw_dbeta)) / (w * w));
    
    double d_dw_dalpha_dbeta = -safe_pow(v, beta - 1.0) * (1.0 + beta * ln_v) * dv_dalpha;
    double d2L7_dalpha_dbeta = (gamma * lambda - 1.0) * ((d_dw_dalpha_dbeta / w) - (dw_dalpha * dw_dbeta) / (w * w));
    
    // --- L8: δ ln(z), where z = 1 - w^λ ---
    double log_w_lambda = lambda * log_w;
    double log_z = log1mexp(log_w_lambda);
    if (!R_finite(log_z)) continue;
    double z = safe_exp(log_z);
    
    double w_lambda_m1 = safe_pow(w, lambda - 1.0);
    double dz_dalpha = -lambda * w_lambda_m1 * dw_dalpha;
    
    double d2z_dalpha2 = -lambda * ((lambda - 1.0) * safe_pow(w, lambda - 2.0) * (dw_dalpha * dw_dalpha)
                                      + w_lambda_m1 * d2w_dalpha2);
    double d2L8_dalpha2 = delta * ((d2z_dalpha2 * z - dz_dalpha * dz_dalpha) / (z * z));
    
    double dz_dbeta = -lambda * w_lambda_m1 * dw_dbeta;
    double d2z_dbeta2 = -lambda * ((lambda - 1.0) * safe_pow(w, lambda - 2.0) * (dw_dbeta * dw_dbeta)
                                     + w_lambda_m1 * d2w_dbeta2);
    double d2L8_dbeta2 = delta * ((d2z_dbeta2 * z - dz_dbeta * dz_dbeta) / (z * z));
    
    double d_dw_dalpha_dbeta_2 = -lambda * ((lambda - 1.0) * safe_pow(w, lambda - 2.0) * dw_dbeta * dw_dalpha
                                              + w_lambda_m1 * d_dw_dalpha_dbeta);
    double d2L8_dalpha_dbeta = delta * ((d_dw_dalpha_dbeta_2 / z) - (dz_dalpha * dz_dbeta) / (z * z));
    
    double dz_dlambda = -safe_pow(w, lambda) * ln_w;
    double d2z_dlambda2 = -safe_pow(w, lambda) * (ln_w * ln_w);
    double d2L8_dlambda2 = delta * ((d2z_dlambda2 * z - dz_dlambda * dz_dlambda) / (z * z));
    
    double d_dalpha_dz_dlambda = -w_lambda_m1 * dw_dalpha - lambda * ln_w * w_lambda_m1 * dw_dalpha;
    double d2L8_dalpha_dlambda = delta * ((d_dalpha_dz_dlambda / z) - (dz_dlambda * dz_dalpha) / (z * z));
    
    double d_dbeta_dz_dlambda = -w_lambda_m1 * dw_dbeta - lambda * ln_w * w_lambda_m1 * dw_dbeta;
    double d2L8_dbeta_dlambda = delta * ((d_dbeta_dz_dlambda / z) - (dz_dlambda * dz_dbeta) / (z * z));
    
    // Validate intermediate results
    if (!R_finite(d2L6_dalpha2) || !R_finite(d2L7_dalpha2) || !R_finite(d2L8_dalpha2) ||
        !R_finite(d2L6_dalpha_dbeta) || !R_finite(d2L7_dalpha_dbeta) || !R_finite(d2L8_dalpha_dbeta) ||
        !R_finite(d2L7_dbeta2) || !R_finite(d2L8_dbeta2) ||
        !R_finite(d2L8_dlambda2) ||
        !R_finite(dw_dalpha) || !R_finite(dw_dbeta) ||
        !R_finite(dz_dalpha) || !R_finite(dz_dbeta) ||
        !R_finite(dz_dlambda)) {
        Rcpp::NumericMatrix nanH(5, 5);
      nanH.fill(R_NaN);
      return nanH;
    }
    
    // ---- Accumulate upper-triangle Hessian contributions ----
    H(0, 0) += d2L6_dalpha2 + d2L7_dalpha2 + d2L8_dalpha2;
    H(0, 1) += d2L6_dalpha_dbeta + d2L7_dalpha_dbeta + d2L8_dalpha_dbeta;
    H(1, 1) += d2L7_dbeta2 + d2L8_dbeta2;
    H(4, 4) += d2L8_dlambda2;
    H(0, 2) += lambda * (dw_dalpha / w);
    H(1, 2) += lambda * (dw_dbeta / w);
    H(0, 3) += dz_dalpha / z;
    H(1, 3) += dz_dbeta / z;

    // λ mixed derivatives
    acc_alpha_lambda += gamma * (dw_dalpha / w) + d2L8_dalpha_dlambda;
    acc_beta_lambda  += gamma * (dw_dbeta  / w) + d2L8_dbeta_dlambda;
    acc_gamma_lambda += ln_w;
    acc_delta_lambda += dz_dlambda / z;
  }

  // Apply accumulated λ mixed derivatives (upper triangle)
  H(0, 4) = acc_alpha_lambda;
  H(1, 4) = acc_beta_lambda;
  H(2, 4) = acc_gamma_lambda;
  H(3, 4) = acc_delta_lambda;

  // Symmetrize once after all accumulations
  H = arma::symmatu(H);

  // Return NEGATIVE Hessian (for minimization of negative log-likelihood)
  return Rcpp::wrap(-H);
}

