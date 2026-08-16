/**
 * @file kkw.cpp
 * @brief Kumaraswamy-Kumaraswamy (KKw) Distribution Functions
 * 
 * @details
 * This file implements the full suite of distribution functions for the
 * four-parameter Kumaraswamy-Kumaraswamy (KKw) distribution, which is a
 * sub-family of the Generalized Kumaraswamy (GKw) distribution obtained
 * by setting γ = 1.
 * 
 * **Relationship to GKw:**
 * \deqn{KKw(\alpha, \beta, \delta, \lambda) = GKw(\alpha, \beta, 1, \delta, \lambda)}
 * 
 * The KKw distribution has probability density function:
 * \deqn{
 *   f(x; \alpha, \beta, \delta, \lambda) = 
 *   \lambda \alpha \beta (\delta+1) x^{\alpha-1} (1-x^\alpha)^{\beta-1}
 *   [1-(1-x^\alpha)^\beta]^{\lambda-1} \{1-[1-(1-x^\alpha)^\beta]^\lambda\}^\delta
 * }
 * for \eqn{x \in (0,1)}.
 * 
 * The cumulative distribution function is:
 * \deqn{
 *   F(x) = 1 - \{1 - [1-(1-x^\alpha)^\beta]^\lambda\}^{\delta+1}
 * }
 * 
 * The quantile function (inverse CDF) is:
 * \deqn{
 *   Q(p) = \left\{1 - \left[1 - \left(1 - (1-p)^{1/(\delta+1)}\right)^{1/\lambda}\right]^{1/\beta}\right\}^{1/\alpha}
 * }
 * 
 * **Parameter Constraints:**
 * - \eqn{\alpha > 0} (shape parameter)
 * - \eqn{\beta > 0} (shape parameter)
 * - \eqn{\delta \geq 0} (shape parameter)
 * - \eqn{\lambda > 0} (shape parameter)
 * 
 * **Special Cases:**
 * | Distribution | Condition | Relation |
 * |--------------|-----------|----------|
 * | Exponentiated Kumaraswamy (EKw) | \eqn{\delta = 0} | KKw with δ=0 |
 * | Kumaraswamy (Kw) | \eqn{\delta = 0, \lambda = 1} | Standard Kumaraswamy |
 * 
 * **Random Variate Generation:**
 * Uses inverse transform method:
 * 1. Generate \eqn{V \sim Uniform(0,1)}
 * 2. Compute \eqn{U = 1 - (1-V)^{1/(\delta+1)}}
 * 3. Return \eqn{X = \{1 - [1 - U^{1/\lambda}]^{1/\beta}\}^{1/\alpha}}
 * 
 * **Numerical Stability:**
 * All computations use log-space arithmetic and numerically stable helper
 * functions from utils.h to prevent overflow/underflow.
 * 
 * **Implemented Functions:**
 * - dkkw(): Probability density function (PDF)
 * - pkkw(): Cumulative distribution function (CDF)
 * - qkkw(): Quantile function (inverse CDF)
 * - rkkw(): Random variate generation
 * - llkkw(): Negative log-likelihood for MLE
 * - grkkw(): Gradient of negative log-likelihood
 * - hskkw(): Hessian of negative log-likelihood
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
 * @brief Probability Density Function of the KKw Distribution
 * 
 * Computes the density (or log-density) for the Kumaraswamy-Kumaraswamy
 * distribution at specified quantiles.
 * 
 * @param x Vector of quantiles (values in (0,1))
 * @param alpha Shape parameter vector (α > 0)
 * @param beta Shape parameter vector (β > 0)
 * @param delta Shape parameter vector (δ ≥ 0)
 * @param lambda Shape parameter vector (λ > 0)
 * @param log_prob If TRUE, returns log-density; otherwise returns density
 * 
 * @return NumericVector of density values (or log-density if log_prob=TRUE)
 * 
 * @details
 * The log-density is computed as:
 * \deqn{
 *   \log f(x) = \log(\lambda) + \log(\alpha) + \log(\beta) + \log(\delta+1)
 *   + (\alpha-1)\log(x) + (\beta-1)\log(1-x^\alpha)
 *   + (\lambda-1)\log(1-(1-x^\alpha)^\beta)
 *   + \delta\log(1-[1-(1-x^\alpha)^\beta]^\lambda)
 * }
 * 
 * @note Exported as .dkkw_cpp for internal package use
 */
// [[Rcpp::export(.dkkw_cpp)]]
Rcpp::NumericVector dkkw(
    const arma::vec& x,
    const Rcpp::NumericVector& alpha,
    const Rcpp::NumericVector& beta,
    const Rcpp::NumericVector& delta,
    const Rcpp::NumericVector& lambda,
    bool log_prob = false
) {
  // Convert R vectors to Armadillo vectors
  arma::vec a_vec(alpha.begin(), alpha.size());
  arma::vec b_vec(beta.begin(), beta.size());
  arma::vec d_vec(delta.begin(), delta.size());
  arma::vec l_vec(lambda.begin(), lambda.size());
  
  // Determine output length for recycling
  size_t N = std::max({x.n_elem, a_vec.n_elem, b_vec.n_elem, 
                      d_vec.n_elem, l_vec.n_elem});
  
  // Initialize result with appropriate default
  arma::vec out(N);
  out.fill(log_prob ? R_NegInf : 0.0);
  
  for (size_t i = 0; i < N; ++i) {
    // Extract recycled parameters
    double a = a_vec[i % a_vec.n_elem];
    double b = b_vec[i % b_vec.n_elem];
    double dd = d_vec[i % d_vec.n_elem];
    double ll = l_vec[i % l_vec.n_elem];
    double xx = x[i % x.n_elem];
    
    // Validate parameters
    if (!check_kkw_pars(a, b, dd, ll)) {
      continue;
    }
    
    // Check support: x must be in (0, 1)
    if (xx <= 0.0 || xx >= 1.0 || !R_finite(xx)) {
      continue;
    }
    
    // ---- Log-space computation of density ----
    
    // Normalization constant: log(λαβ(δ+1))
    double logCst = safe_log(ll) + safe_log(a) + safe_log(b) + safe_log(dd + 1.0);
    
    // Compute log(x) and log(x^α)
    double lx = safe_log(xx);
    double log_xalpha = a * lx;
    
    // Compute log(1 - x^α) using stable log1mexp
    double log_1_minus_xalpha = log1mexp(log_xalpha);
    if (!R_finite(log_1_minus_xalpha)) {
      continue;
    }
    
    // Term: (β-1) * log(1 - x^α)
    double term1 = (b - 1.0) * log_1_minus_xalpha;
    
    // Compute A = (1 - x^α)^β → log(A) = β * log(1 - x^α)
    double logA = b * log_1_minus_xalpha;
    
    // Compute log(1 - A) = log(1 - (1-x^α)^β) using log1mexp
    double log_1_minusA = log1mexp(logA);
    if (!R_finite(log_1_minusA)) {
      continue;
    }
    
    // Term: (λ-1) * log(1 - A)
    double term2 = (ll - 1.0) * log_1_minusA;
    
    // Compute B = [1 - (1-x^α)^β]^λ → log(B) = λ * log(1 - A)
    double logB = ll * log_1_minusA;
    
    // Compute log(1 - B) using log1mexp
    double log_1_minus_B = log1mexp(logB);
    if (!R_finite(log_1_minus_B)) {
      continue;
    }
    
    // Term: δ * log(1 - B)
    double term3 = dd * log_1_minus_B;
    
    // Assemble log-density
    double log_pdf = logCst + (a - 1.0) * lx + term1 + term2 + term3;
    
    // Return appropriate scale
    out(i) = log_prob ? log_pdf : safe_exp(log_pdf);
  }
  
  return Rcpp::NumericVector(out.memptr(), out.memptr() + out.n_elem);
}


// ============================================================================
// CUMULATIVE DISTRIBUTION FUNCTION
// ============================================================================

/**
 * @brief Cumulative Distribution Function of the KKw Distribution
 * 
 * Computes the cumulative probability for the Kumaraswamy-Kumaraswamy
 * distribution at specified quantiles.
 * 
 * @param q Vector of quantiles
 * @param alpha Shape parameter vector (α > 0)
 * @param beta Shape parameter vector (β > 0)
 * @param delta Shape parameter vector (δ ≥ 0)
 * @param lambda Shape parameter vector (λ > 0)
 * @param lower_tail If TRUE, returns P(X ≤ q); otherwise P(X > q)
 * @param log_p If TRUE, returns log-probability
 * 
 * @return NumericVector of cumulative probabilities
 * 
 * @details
 * The CDF is computed as:
 * \deqn{F(x) = 1 - \{1 - [1-(1-x^\alpha)^\beta]^\lambda\}^{\delta+1}}
 * 
 * @note Exported as .pkkw_cpp for internal package use
 */
// [[Rcpp::export(.pkkw_cpp)]]
Rcpp::NumericVector pkkw(
    const arma::vec& q,
    const Rcpp::NumericVector& alpha,
    const Rcpp::NumericVector& beta,
    const Rcpp::NumericVector& delta,
    const Rcpp::NumericVector& lambda,
    bool lower_tail = true,
    bool log_p = false
) {
  // Convert R vectors to Armadillo vectors
  arma::vec a_vec(alpha.begin(), alpha.size());
  arma::vec b_vec(beta.begin(), beta.size());
  arma::vec d_vec(delta.begin(), delta.size());
  arma::vec l_vec(lambda.begin(), lambda.size());
  
  // Determine output length for recycling
  size_t N = std::max({q.n_elem, a_vec.n_elem, b_vec.n_elem, 
                      d_vec.n_elem, l_vec.n_elem});
  
  arma::vec out(N);
  
  for (size_t i = 0; i < N; ++i) {
    // Extract recycled parameters
    double a = a_vec[i % a_vec.n_elem];
    double b = b_vec[i % b_vec.n_elem];
    double dd = d_vec[i % d_vec.n_elem];
    double ll = l_vec[i % l_vec.n_elem];
    double xx = q[i % q.n_elem];
    
    // Validate parameters
    if (!check_kkw_pars(a, b, dd, ll)) {
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
    double log_xalpha = a * safe_log(xx);
    double xalpha = safe_exp(log_xalpha);
    
    // Step 2: 1 - x^α
    double one_minus_xalpha = 1.0 - xalpha;
    if (one_minus_xalpha <= 0.0) {
      double val1 = lower_tail ? 1.0 : 0.0;
      out(i) = log_p ? safe_log(val1) : val1;
      continue;
    }
    
    // Step 3: (1 - x^α)^β
    double vbeta = safe_pow(one_minus_xalpha, b);
    
    // Step 4: y = 1 - (1 - x^α)^β
    double y = 1.0 - vbeta;
    if (y <= 0.0) {
      double val0 = lower_tail ? 0.0 : 1.0;
      out(i) = log_p ? safe_log(val0) : val0;
      continue;
    }
    if (y >= 1.0) {
      double val1 = lower_tail ? 1.0 : 0.0;
      out(i) = log_p ? safe_log(val1) : val1;
      continue;
    }
    
    // Step 5: y^λ = [1-(1-x^α)^β]^λ
    double ylambda = safe_pow(y, ll);
    if (ylambda <= 0.0) {
      double val0 = lower_tail ? 0.0 : 1.0;
      out(i) = log_p ? safe_log(val0) : val0;
      continue;
    }
    if (ylambda >= 1.0) {
      double val1 = lower_tail ? 1.0 : 0.0;
      out(i) = log_p ? safe_log(val1) : val1;
      continue;
    }
    
    // Step 6: F(x) = 1 - (1 - y^λ)^(δ+1)
    double outer = 1.0 - ylambda;
    double cdfval = 1.0 - safe_pow(outer, dd + 1.0);
    
    // Apply tail adjustment
    if (!lower_tail) {
      cdfval = 1.0 - cdfval;
    }
    
    // Apply log transformation
    if (log_p) {
      cdfval = safe_log(cdfval);
    }
    
    out(i) = cdfval;
  }
  
  return Rcpp::NumericVector(out.memptr(), out.memptr() + out.n_elem);
}


// ============================================================================
// QUANTILE FUNCTION
// ============================================================================

/**
 * @brief Quantile Function (Inverse CDF) of the KKw Distribution
 * 
 * Computes quantiles for the Kumaraswamy-Kumaraswamy distribution
 * given probability values.
 * 
 * @param p Vector of probabilities (values in [0,1])
 * @param alpha Shape parameter vector (α > 0)
 * @param beta Shape parameter vector (β > 0)
 * @param delta Shape parameter vector (δ ≥ 0)
 * @param lambda Shape parameter vector (λ > 0)
 * @param lower_tail If TRUE, probabilities are P(X ≤ x); otherwise P(X > x)
 * @param log_p If TRUE, probabilities are given as log(p)
 * 
 * @return NumericVector of quantiles
 * 
 * @details
 * The quantile function inverts the CDF:
 * \deqn{Q(p) = \left\{1 - \left[1 - \left(1 - (1-p)^{1/(\delta+1)}\right)^{1/\lambda}\right]^{1/\beta}\right\}^{1/\alpha}}
 * 
 * @note Exported as .qkkw_cpp for internal package use
 */
// [[Rcpp::export(.qkkw_cpp)]]
Rcpp::NumericVector qkkw(
    const arma::vec& p,
    const Rcpp::NumericVector& alpha,
    const Rcpp::NumericVector& beta,
    const Rcpp::NumericVector& delta,
    const Rcpp::NumericVector& lambda,
    bool lower_tail = true,
    bool log_p = false
) {
  // Convert R vectors to Armadillo vectors
  arma::vec a_vec(alpha.begin(), alpha.size());
  arma::vec b_vec(beta.begin(), beta.size());
  arma::vec d_vec(delta.begin(), delta.size());
  arma::vec l_vec(lambda.begin(), lambda.size());
  
  // Determine output length for recycling
  size_t N = std::max({p.n_elem, a_vec.n_elem, b_vec.n_elem, 
                      d_vec.n_elem, l_vec.n_elem});
  
  arma::vec out(N);
  
  for (size_t i = 0; i < N; ++i) {
    // Extract recycled parameters
    double a = a_vec[i % a_vec.n_elem];
    double b = b_vec[i % b_vec.n_elem];
    double dd = d_vec[i % d_vec.n_elem];
    double ll = l_vec[i % l_vec.n_elem];
    double pp = p[i % p.n_elem];
    
    // Validate parameters
    if (!check_kkw_pars(a, b, dd, ll)) {
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
    
    // Step 1: tmp1 = 1 - (1-p)^(1/(δ+1))
    double tmp1 = 1.0 - safe_pow(1.0 - pp, 1.0 / (dd + 1.0));
    tmp1 = std::max(0.0, std::min(1.0, tmp1));
    
    // Step 2: T = tmp1^(1/λ)
    double T = (ll == 1.0) ? tmp1 : safe_pow(tmp1, 1.0 / ll);
    
    // Step 3: M = 1 - T  →  (1 - x^α)^β = M
    double M = 1.0 - T;
    M = std::max(0.0, std::min(1.0, M));
    
    // Step 4: Mpow = M^(1/β)  →  1 - x^α = Mpow
    double Mpow = safe_pow(M, 1.0 / b);
    
    // Step 5: xalpha = 1 - Mpow  →  x^α = xalpha
    double xalpha = 1.0 - Mpow;
    xalpha = std::max(0.0, std::min(1.0, xalpha));
    
    // Step 6: x = xalpha^(1/α)
    double xx = (a == 1.0) ? xalpha : safe_pow(xalpha, 1.0 / a);
    
    // Clamp to valid support
    xx = std::max(0.0, std::min(1.0, xx));
    
    out(i) = xx;
  }
  
  return Rcpp::NumericVector(out.memptr(), out.memptr() + out.n_elem);
}


// ============================================================================
// RANDOM NUMBER GENERATION
// ============================================================================

/**
 * @brief Random Variate Generation for the KKw Distribution
 * 
 * Generates random samples from the Kumaraswamy-Kumaraswamy distribution
 * using the inverse transform method.
 * 
 * @param n Number of random variates to generate
 * @param alpha Shape parameter vector (α > 0)
 * @param beta Shape parameter vector (β > 0)
 * @param delta Shape parameter vector (δ ≥ 0)
 * @param lambda Shape parameter vector (λ > 0)
 * 
 * @return NumericVector of n random variates from KKw distribution
 * 
 * @details
 * Algorithm:
 * 1. Generate V ~ Uniform(0,1)
 * 2. U = 1 - (1-V)^(1/(δ+1))
 * 3. X = {1 - [1 - U^(1/λ)]^(1/β)}^(1/α)
 * 
 * @note Exported as .rkkw_cpp for internal package use
 */
// [[Rcpp::export(.rkkw_cpp)]]
Rcpp::NumericVector rkkw(
    int n,
    const Rcpp::NumericVector& alpha,
    const Rcpp::NumericVector& beta,
    const Rcpp::NumericVector& delta,
    const Rcpp::NumericVector& lambda
) {
  if (n <= 0) {
    Rcpp::stop("rkkw: n must be positive");
  }
  
  // Convert R vectors to Armadillo vectors
  arma::vec a_vec(alpha.begin(), alpha.size());
  arma::vec b_vec(beta.begin(), beta.size());
  arma::vec d_vec(delta.begin(), delta.size());
  arma::vec l_vec(lambda.begin(), lambda.size());
  
  arma::vec out(n);
  
  for (int i = 0; i < n; i++) {
    // Extract recycled parameters (direct modulo, no intermediate variable)
    double a = a_vec[i % a_vec.n_elem];
    double b = b_vec[i % b_vec.n_elem];
    double dd = d_vec[i % d_vec.n_elem];
    double ll = l_vec[i % l_vec.n_elem];
    
    // Validate parameters
    if (!check_kkw_pars(a, b, dd, ll)) {
      out(i) = NA_REAL;
      Rcpp::warning("rkkw: invalid parameters at index %d", i + 1);
      continue;
    }
    
    // Generate V ~ Uniform(0,1)
    double V = R::runif(0.0, 1.0);
    
    // Step 1: U = 1 - (1-V)^(1/(δ+1))
    double U = 1.0 - safe_pow(1.0 - V, 1.0 / (dd + 1.0));
    U = std::max(0.0, std::min(1.0, U));
    
    // Step 2: u_pow = U^(1/λ)
    double u_pow = (ll == 1.0) ? U : safe_pow(U, 1.0 / ll);
    
    // Step 3: bracket = 1 - u_pow
    double bracket = 1.0 - u_pow;
    bracket = std::max(0.0, std::min(1.0, bracket));
    
    // Step 4: bracket2 = bracket^(1/β)
    double bracket2 = safe_pow(bracket, 1.0 / b);
    
    // Step 5: xalpha = 1 - bracket2
    double xalpha = 1.0 - bracket2;
    xalpha = std::max(0.0, std::min(1.0, xalpha));
    
    // Step 6: x = xalpha^(1/α)
    double xx;
    if (a == 1.0) {
      xx = xalpha;
    } else {
      xx = safe_pow(xalpha, 1.0 / a);
      if (!R_finite(xx) || xx < 0.0) xx = 0.0;
      if (xx > 1.0) xx = 1.0;
    }
    
    out(i) = xx;
  }
  
  return Rcpp::NumericVector(out.memptr(), out.memptr() + out.n_elem);
}


// ============================================================================
// NEGATIVE LOG-LIKELIHOOD FUNCTION
// ============================================================================

/**
 * @brief Negative Log-Likelihood for KKw Distribution
 * 
 * Computes the negative log-likelihood function for parameter estimation
 * via maximum likelihood.
 * 
 * @param par Parameter vector of length 4: (α, β, δ, λ)
 * @param data Vector of observations (must be in (0,1))
 * 
 * @return Negative log-likelihood value (scalar)
 * 
 * @details
 * The log-likelihood for n observations is:
 * \deqn{
 *   \ell(\theta) = n[\ln\lambda + \ln\alpha + \ln\beta + \ln(\delta+1)]
 *   + (\alpha-1)\sum\ln x_i + (\beta-1)\sum\ln v_i
 *   + (\lambda-1)\sum\ln w_i + \delta\sum\ln z_i
 * }
 * where:
 * - \eqn{v_i = 1 - x_i^\alpha}
 * - \eqn{w_i = 1 - v_i^\beta}
 * - \eqn{z_i = 1 - w_i^\lambda}
 * 
 * Returns +Inf for invalid parameters or data outside (0,1).
 * 
 * @note Exported as .llkkw_cpp for internal package use
 */
// [[Rcpp::export(.llkkw_cpp)]]
double llkkw(const Rcpp::NumericVector& par, const Rcpp::NumericVector& data) {
  // Validate parameter vector length
  if (par.size() < 4) return R_PosInf;
  
  // Extract parameters
  double alpha = par[0];
  double beta = par[1];
  double delta = par[2];
  double lambda = par[3];
  
  // Validate parameters using consistent checker
  if (!check_kkw_pars(alpha, beta, delta, lambda)) return R_PosInf;
  
  // Convert and validate data
  arma::vec x = Rcpp::as<arma::vec>(data);
  if (x.n_elem < 1 || arma::any(x <= 0.0) || arma::any(x >= 1.0)) return R_PosInf;
  
  int n = x.n_elem;
  
  // Numerical stability constants
  const double eps = 1e-10;
  
  // Special case: δ ≈ 0 reduces to EKw
  bool is_ekw = (delta < eps);
  
  // Constant term: n * [log(λ) + log(α) + log(β) + log(δ+1)]
  double log_alpha = safe_log(alpha);
  double log_beta = safe_log(beta);
  double log_lambda = safe_log(lambda);
  double log_delta_plus_1 = is_ekw ? 0.0 : std::log1p(delta);
  
  double const_term = n * (log_lambda + log_alpha + log_beta + log_delta_plus_1);
  
  // Initialize accumulators
  double sum_term1 = 0.0;  // (α-1) * Σlog(x)
  double sum_term2 = 0.0;  // (β-1) * Σlog(v)
  double sum_term3 = 0.0;  // (λ-1) * Σlog(w)
  double sum_term4 = 0.0;  // δ * Σlog(z)
  
  for (int i = 0; i < n; i++) {
    double xi = x(i);
    
    // Clamp to safe bounds
    xi = std::max(eps, std::min(1.0 - eps, xi));
    
    double log_xi = std::log(xi);
    sum_term1 += (alpha - 1.0) * log_xi;
    
    // Compute x^α stably
    double log_x_alpha = alpha * log_xi;
    double x_alpha = safe_exp(log_x_alpha);
    
    // Compute v = 1 - x^α
    double v, log_v;
    if (x_alpha > 0.9999) {
      v = -std::expm1(log_x_alpha);
      log_v = std::log(std::max(v, eps));
    } else {
      v = 1.0 - x_alpha;
      log_v = std::log1p(-x_alpha);
    }
    v = std::max(v, eps);
    
    sum_term2 += (beta - 1.0) * log_v;
    
    // Compute v^β stably
    double log_v_beta = beta * log_v;
    double v_beta = safe_exp(log_v_beta);
    
    // Compute w = 1 - v^β
    double w, log_w;
    if (v_beta > 0.9999) {
      w = -std::expm1(log_v_beta);
      log_w = std::log(std::max(w, eps));
    } else {
      w = 1.0 - v_beta;
      log_w = std::log1p(-v_beta);
    }
    w = std::max(w, eps);
    
    sum_term3 += (lambda - 1.0) * log_w;
    
    // Compute z term only if δ > 0
    if (!is_ekw) {
      // Compute w^λ stably
      double log_w_lambda = lambda * log_w;
      double w_lambda = safe_exp(log_w_lambda);
      
      // Compute z = 1 - w^λ
      double z, log_z;
      if (w_lambda > 0.9999) {
        z = -std::expm1(log_w_lambda);
        log_z = std::log(std::max(z, eps));
      } else {
        z = 1.0 - w_lambda;
        log_z = std::log1p(-w_lambda);
      }
      z = std::max(z, eps);
      log_z = std::log(z);  // Recompute after clamping
      
      // Scale for very large delta
      double effective_delta = std::min(delta, 1000.0);
      sum_term4 += effective_delta * log_z;
    }
  }
  
  double loglike = const_term + sum_term1 + sum_term2 + sum_term3 + sum_term4;
  
  // Guard against NaN/Inf
  if (!std::isfinite(loglike)) {
    return R_PosInf;
  }
  
  return -loglike;
}


// ============================================================================
// GRADIENT OF NEGATIVE LOG-LIKELIHOOD
// ============================================================================

/**
 * @brief Gradient of Negative Log-Likelihood for KKw Distribution
 * 
 * Computes the gradient vector of the negative log-likelihood for
 * optimization-based parameter estimation.
 * 
 * @param par Parameter vector of length 4: (α, β, δ, λ)
 * @param data Vector of observations (must be in (0,1))
 * 
 * @return NumericVector of length 4 containing partial derivatives
 *         with respect to (α, β, δ, λ)
 * 
 * @details
 * Computes analytical gradients using chain rule and log-space arithmetic
 * for numerical stability. Returns NaN vector for invalid inputs.
 * 
 * @note Exported as .grkkw_cpp for internal package use
 */
// [[Rcpp::export(.grkkw_cpp)]]
Rcpp::NumericVector grkkw(const Rcpp::NumericVector& par, const Rcpp::NumericVector& data) {
  // Validate parameter vector length
  if (par.size() < 4) {
    return Rcpp::NumericVector(4, R_NaN);
  }
  
  // Extract parameters
  double alpha = par[0];
  double beta = par[1];
  double delta = par[2];
  double lambda = par[3];
  
  // Validate parameters using consistent checker
  if (!check_kkw_pars(alpha, beta, delta, lambda)) {
    return Rcpp::NumericVector(4, R_NaN);
  }
  
  // Convert and validate data
  arma::vec x = Rcpp::as<arma::vec>(data);
  if (x.n_elem < 1 || arma::any(x <= 0) || arma::any(x >= 1)) {
    return Rcpp::NumericVector(4, R_NaN);
  }
  
  int n = x.n_elem;
  Rcpp::NumericVector grad(4, 0.0);
  
  // Numerical stability constants
  const double eps = 1e-10;
  
  // Initialize gradient accumulators
  double d_alpha = n / alpha;
  double d_beta = n / beta;
  double d_delta = n / (delta + 1.0);
  double d_lambda = n / lambda;
  
  // Special case: δ ≈ 0 reduces to EKw
  bool delta_near_zero = (delta < eps);
  
  for (int i = 0; i < n; i++) {
    double xi = x(i);
    
    // Clamp to safe bounds
    xi = std::max(eps, std::min(1.0 - eps, xi));
    
    double log_xi = std::log(xi);
    d_alpha += log_xi;
    
    // Compute x^α stably
    double log_x_alpha = alpha * log_xi;
    double x_alpha = safe_exp(log_x_alpha);
    double x_alpha_log_x = x_alpha * log_xi;
    
    // Compute v = 1 - x^α
    double v, log_v;
    if (x_alpha > 0.9995) {
      v = -std::expm1(log_x_alpha);
      log_v = std::log(std::max(v, eps));
    } else {
      v = 1.0 - x_alpha;
      log_v = std::log1p(-x_alpha);
    }
    v = std::max(v, eps);
    d_beta += log_v;
    
    // Compute v^β and v^(β-1) stably
    double log_v_beta = beta * log_v;
    double v_beta = safe_exp(log_v_beta);
    double v_beta_m1 = safe_exp((beta - 1.0) * log_v);
    double v_beta_log_v = v_beta * log_v;
    
    // Compute w = 1 - v^β
    double w, log_w;
    if (v_beta > 0.9995) {
      w = -std::expm1(log_v_beta);
      log_w = std::log(std::max(w, eps));
    } else {
      w = 1.0 - v_beta;
      log_w = std::log1p(-v_beta);
    }
    w = std::max(w, eps);
    d_lambda += log_w;
    
    // Compute gradient terms
    double lambda_factor = (std::abs(lambda - 1.0) > eps) ? (lambda - 1.0) : 0.0;
    
    // Term for α gradient: (β-1)/v
    double alpha_term1 = (std::abs(beta - 1.0) > eps) ? ((beta - 1.0) / v) : 0.0;
    
    // Term for α gradient: (λ-1)*β*v^(β-1)/w
    double alpha_term2 = (lambda_factor != 0.0) ? (lambda_factor * beta * v_beta_m1 / w) : 0.0;
    
    // For δ > 0, compute additional terms
    if (!delta_near_zero) {
      // Compute w^λ and w^(λ-1) stably
      double log_w_lambda = lambda * log_w;
      double w_lambda = safe_exp(log_w_lambda);
      double w_lambda_m1 = safe_exp((lambda - 1.0) * log_w);
      double w_lambda_log_w = w_lambda * log_w;
      
      // Compute z = 1 - w^λ
      double z, log_z;
      if (w_lambda > 0.9995) {
        z = -std::expm1(log_w_lambda);
      } else {
        z = 1.0 - w_lambda;
      }
      z = std::max(z, eps);
      log_z = std::log(z);  // Recompute after clamping (BUG FIX)
      
      d_delta += log_z;
      
      // Term for α gradient: δ*λ*β*v^(β-1)*w^(λ-1)/z
      double alpha_term3 = std::min(delta * lambda * beta * v_beta_m1 * w_lambda_m1 / z, 1000.0);
      
      // Combine terms for α gradient
      d_alpha -= x_alpha_log_x * (alpha_term1 - alpha_term2 + alpha_term3);
      
      // Term for β gradient: (λ-1)/w
      double beta_term1 = (lambda_factor != 0.0) ? (lambda_factor / w) : 0.0;
      
      // Term for β gradient: δ*λ*w^(λ-1)/z
      double beta_term2 = std::min(delta * lambda * w_lambda_m1 / z, 1000.0);
      
      // Combine terms for β gradient
      d_beta -= v_beta_log_v * (beta_term1 - beta_term2);
      
      // λ gradient term: δ*(w^λ*log(w))/z
      d_lambda -= delta * w_lambda_log_w / z;
    } else {
      // Simplified for δ ≈ 0 (EKw case)
      d_alpha -= x_alpha_log_x * (alpha_term1 - alpha_term2);
      
      if (lambda_factor != 0.0) {
        d_beta -= v_beta_log_v * (lambda_factor / w);
      }
    }
  }
  
  // Return NEGATIVE gradient (for minimization)
  grad[0] = -d_alpha;
  grad[1] = -d_beta;
  grad[2] = -d_delta;
  grad[3] = -d_lambda;
  
  return grad;
}


// ============================================================================
// HESSIAN OF NEGATIVE LOG-LIKELIHOOD
// ============================================================================

/**
 * @brief Hessian Matrix of Negative Log-Likelihood for KKw Distribution
 * 
 * Computes the Hessian matrix (matrix of second partial derivatives) of
 * the negative log-likelihood for standard error estimation and
 * optimization algorithms.
 * 
 * @param par Parameter vector of length 4: (α, β, δ, λ)
 * @param data Vector of observations (must be in (0,1))
 * 
 * @return NumericMatrix of dimension 4×4 containing the Hessian
 * 
 * @details
 * Computes analytical second derivatives. The Hessian is symmetric.
 * Parameter ordering: (α, β, δ, λ) → indices (0, 1, 2, 3).
 * 
 * Returns NaN matrix for invalid inputs.
 * 
 * @note Exported as .hskkw_cpp for internal package use
 */
// [[Rcpp::export(.hskkw_cpp)]]
Rcpp::NumericMatrix hskkw(const Rcpp::NumericVector& par, const Rcpp::NumericVector& data) {
  // Validate parameter vector length
  if (par.size() < 4) {
    Rcpp::NumericMatrix nanH(4, 4);
    nanH.fill(R_NaN);
    return nanH;
  }
  
  // Extract parameters
  double alpha = par[0];
  double beta = par[1];
  double delta = par[2];
  double lambda = par[3];
  
  // Validate parameters using consistent checker
  if (!check_kkw_pars(alpha, beta, delta, lambda)) {
    Rcpp::NumericMatrix nanH(4, 4);
    nanH.fill(R_NaN);
    return nanH;
  }
  
  // Convert and validate data
  arma::vec x = Rcpp::as<arma::vec>(data);
  if (x.n_elem < 1 || arma::any(x <= 0) || arma::any(x >= 1)) {
    Rcpp::NumericMatrix nanH(4, 4);
    nanH.fill(R_NaN);
    return nanH;
  }
  
  int n = x.n_elem;
  
  // Numerical stability constants
  const double eps = 1e-10;
  const double max_contrib = 1e6;
  
  // Initialize Hessian matrix
  arma::mat H(4, 4, arma::fill::zeros);
  
  // Special case: δ ≈ 0 reduces to EKw
  bool delta_near_zero = (delta < eps);
  
  // Constant diagonal terms
  H(0, 0) = -n / (alpha * alpha);
  H(1, 1) = -n / (beta * beta);
  H(3, 3) = -n / (lambda * lambda);
  H(2, 2) = delta_near_zero ? -n : -n / std::pow(delta + 1.0, 2.0);
  
  // Accumulators for mixed derivatives with λ
  double acc_alpha_lambda = 0.0;
  double acc_beta_lambda = 0.0;
  double acc_delta_lambda = 0.0;
  
  for (int i = 0; i < n; i++) {
    double xi = x(i);
    
    // Clamp to safe bounds
    xi = std::max(eps, std::min(1.0 - eps, xi));
    
    double log_xi = std::log(xi);
    
    // --- Compute A = x^α and derivatives ---
    double log_A = alpha * log_xi;
    double A = safe_exp(log_A);
    double dA_dalpha = A * log_xi;
    double d2A_dalpha2 = A * log_xi * log_xi;
    
    // --- Compute v = 1 - A and derivatives ---
    double v, log_v;
    if (A > 0.9995) {
      v = -std::expm1(log_A);
      log_v = std::log(std::max(v, eps));
    } else {
      v = 1.0 - A;
      log_v = std::log1p(-A);
    }
    v = std::max(v, eps);
    
    double dv_dalpha = -dA_dalpha;
    double d2v_dalpha2 = -d2A_dalpha2;
    
    // --- Terms for (β-1)ln(v) ---
    double d2L6_dalpha2 = 0.0;
    double d2L6_dalpha_dbeta = 0.0;
    
    if (std::abs(beta - 1.0) > eps) {
      double v_squared = v * v;
      d2L6_dalpha2 = (beta - 1.0) * ((d2v_dalpha2 * v - dv_dalpha * dv_dalpha) / v_squared);
      d2L6_dalpha_dbeta = dv_dalpha / v;
    }
    
    // --- Compute w = 1 - v^β and derivatives ---
    double log_v_beta = beta * log_v;
    double v_beta = safe_exp(log_v_beta);
    double v_beta_m1 = safe_exp((beta - 1.0) * log_v);
    double v_beta_m2 = safe_exp((beta - 2.0) * log_v);
    
    double w, log_w;
    if (v_beta > 0.9995) {
      w = -std::expm1(log_v_beta);
      log_w = std::log(std::max(w, eps));
    } else {
      w = 1.0 - v_beta;
      log_w = std::log1p(-v_beta);
    }
    w = std::max(w, eps);
    double w_squared = w * w;
    
    // Derivatives of w
    double dw_dv = -beta * v_beta_m1;
    double dw_dalpha = dw_dv * dv_dalpha;
    double dw_dbeta = -v_beta * log_v;
    
    // Second derivatives of w
    double d2w_dalpha2 = -beta * ((beta - 1.0) * v_beta_m2 * (dv_dalpha * dv_dalpha) +
                                  v_beta_m1 * d2v_dalpha2);
    double d2w_dbeta2 = -v_beta * (log_v * log_v);
    double d_dw_dalpha_dbeta = -v_beta_m1 * (1.0 + beta * log_v) * dv_dalpha;
    
    // --- Terms for (λ-1)ln(w) ---
    double d2L7_dalpha2 = 0.0;
    double d2L7_dbeta2 = 0.0;
    double d2L7_dalpha_dbeta = 0.0;
    
    double lambda_minus_1 = lambda - 1.0;
    if (std::abs(lambda_minus_1) > eps) {
      d2L7_dalpha2 = lambda_minus_1 * ((d2w_dalpha2 * w - (dw_dalpha * dw_dalpha)) / w_squared);
      d2L7_dbeta2 = lambda_minus_1 * ((d2w_dbeta2 * w - (dw_dbeta * dw_dbeta)) / w_squared);
      d2L7_dalpha_dbeta = lambda_minus_1 * ((d_dw_dalpha_dbeta / w) - (dw_dalpha * dw_dbeta) / w_squared);
      
      // Clamp extreme values
      d2L7_dalpha2 = std::min(std::max(d2L7_dalpha2, -max_contrib), max_contrib);
      d2L7_dbeta2 = std::min(std::max(d2L7_dbeta2, -max_contrib), max_contrib);
      d2L7_dalpha_dbeta = std::min(std::max(d2L7_dalpha_dbeta, -max_contrib), max_contrib);
    }
    
    // For δ ≈ 0, skip z calculations (EKw case)
    if (delta_near_zero) {
      // Update Hessian elements
      H(0, 0) += d2L6_dalpha2 + d2L7_dalpha2;
      H(0, 1) += d2L6_dalpha_dbeta + d2L7_dalpha_dbeta;
      H(1, 0) = H(0, 1);
      H(1, 1) += d2L7_dbeta2;
      
      // Mixed derivatives with λ (from L7 only)
      acc_alpha_lambda += dw_dalpha / w;
      acc_beta_lambda += dw_dbeta / w;
    } else {
      // --- Compute z = 1 - w^λ and derivatives ---
      double log_w_lambda = lambda * log_w;
      double w_lambda = safe_exp(log_w_lambda);
      double w_lambda_m1 = safe_exp((lambda - 1.0) * log_w);
      double w_lambda_m2 = safe_exp((lambda - 2.0) * log_w);
      
      double z;
      if (w_lambda > 0.9995) {
        z = -std::expm1(log_w_lambda);
      } else {
        z = 1.0 - w_lambda;
      }
      z = std::max(z, eps);
      double z_squared = z * z;
      
      // First derivatives of z
      double dz_dalpha = -lambda * w_lambda_m1 * dw_dalpha;
      double dz_dbeta = -lambda * w_lambda_m1 * dw_dbeta;
      double dz_dlambda = -w_lambda * log_w;
      
      // Second derivatives of z
      double d2z_dalpha2 = -lambda * ((lambda - 1.0) * w_lambda_m2 * (dw_dalpha * dw_dalpha) +
                                      w_lambda_m1 * d2w_dalpha2);
      double d2z_dbeta2 = -lambda * ((lambda - 1.0) * w_lambda_m2 * (dw_dbeta * dw_dbeta) +
                                     w_lambda_m1 * d2w_dbeta2);
      double d2z_dlambda2 = -w_lambda * (log_w * log_w);
      
      // Mixed derivatives of z
      double d_dw_dalpha_dbeta_2 = -lambda * ((lambda - 1.0) * w_lambda_m2 * dw_dbeta * dw_dalpha +
                                              w_lambda_m1 * d_dw_dalpha_dbeta);
      double d_dalpha_dz_dlambda = -lambda * w_lambda_m1 * dw_dalpha * log_w - w_lambda * (dw_dalpha / w);
      double d_dbeta_dz_dlambda = -lambda * w_lambda_m1 * dw_dbeta * log_w - w_lambda * (dw_dbeta / w);
      
      // Terms for δln(z)
      double d2L8_dalpha2 = delta * ((d2z_dalpha2 * z - dz_dalpha * dz_dalpha) / z_squared);
      double d2L8_dbeta2 = delta * ((d2z_dbeta2 * z - dz_dbeta * dz_dbeta) / z_squared);
      double d2L8_dlambda2 = delta * ((d2z_dlambda2 * z - dz_dlambda * dz_dlambda) / z_squared);
      double d2L8_dalpha_dbeta = delta * ((d_dw_dalpha_dbeta_2 / z) - (dz_dalpha * dz_dbeta) / z_squared);
      double d2L8_dalpha_dlambda = delta * ((d_dalpha_dz_dlambda / z) - (dz_dlambda * dz_dalpha) / z_squared);
      double d2L8_dbeta_dlambda = delta * ((d_dbeta_dz_dlambda / z) - (dz_dlambda * dz_dbeta) / z_squared);
      double d2L8_ddelta_dlambda = dz_dlambda / z;
      
      // Clamp extreme values
      d2L8_dalpha2 = std::min(std::max(d2L8_dalpha2, -max_contrib), max_contrib);
      d2L8_dbeta2 = std::min(std::max(d2L8_dbeta2, -max_contrib), max_contrib);
      d2L8_dlambda2 = std::min(std::max(d2L8_dlambda2, -max_contrib), max_contrib);
      d2L8_dalpha_dbeta = std::min(std::max(d2L8_dalpha_dbeta, -max_contrib), max_contrib);
      d2L8_dalpha_dlambda = std::min(std::max(d2L8_dalpha_dlambda, -max_contrib), max_contrib);
      d2L8_dbeta_dlambda = std::min(std::max(d2L8_dbeta_dlambda, -max_contrib), max_contrib);
      
      // Update Hessian elements
      H(0, 0) += d2L6_dalpha2 + d2L7_dalpha2 + d2L8_dalpha2;
      H(0, 1) += d2L6_dalpha_dbeta + d2L7_dalpha_dbeta + d2L8_dalpha_dbeta;
      H(1, 0) = H(0, 1);
      H(1, 1) += d2L7_dbeta2 + d2L8_dbeta2;
      H(3, 3) += d2L8_dlambda2;
      
      // Mixed derivatives with δ
      H(0, 2) += dz_dalpha / z;
      H(2, 0) = H(0, 2);
      H(1, 2) += dz_dbeta / z;
      H(2, 1) = H(1, 2);
      
      // Accumulators for mixed derivatives with λ
      acc_alpha_lambda += (dw_dalpha / w) + d2L8_dalpha_dlambda;
      acc_beta_lambda += (dw_dbeta / w) + d2L8_dbeta_dlambda;
      acc_delta_lambda += d2L8_ddelta_dlambda;
    }
  }
  
  // Apply mixed derivatives with λ
  H(0, 3) = acc_alpha_lambda;
  H(3, 0) = H(0, 3);
  H(1, 3) = acc_beta_lambda;
  H(3, 1) = H(1, 3);
  H(2, 3) = acc_delta_lambda;
  H(3, 2) = H(2, 3);
  
  // Enforce symmetry
  for (int i = 0; i < 4; i++) {
    for (int j = i + 1; j < 4; j++) {
      double avg = (H(i, j) + H(j, i)) / 2.0;
      H(i, j) = H(j, i) = avg;
    }
  }
  
  // Return NEGATIVE Hessian (for minimization)
  return Rcpp::wrap(-H);
}

