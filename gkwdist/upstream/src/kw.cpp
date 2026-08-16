/**
 * @file kw.cpp
 * @brief Kumaraswamy (Kw) Distribution Functions
 * 
 * @details
 * This file implements the full suite of distribution functions for the
 * two-parameter Kumaraswamy (Kw) distribution, which is the base case
 * of the entire GKw family (setting γ=1, δ=0, λ=1).
 * 
 * **Relationship to GKw:**
 * \deqn{Kw(\alpha, \beta) = GKw(\alpha, \beta, 1, 0, 1)}
 * 
 * The Kumaraswamy distribution has probability density function:
 * \deqn{
 *   f(x; \alpha, \beta) = \alpha \beta x^{\alpha-1} (1-x^\alpha)^{\beta-1}
 * }
 * for \eqn{x \in (0,1)}.
 * 
 * **Historical Context:**
 * Introduced by Kumaraswamy (1980) as a tractable alternative to the Beta
 * distribution for modeling double-bounded random processes. Unlike the Beta
 * distribution, both the PDF and CDF have closed-form expressions.
 * 
 * **Derivation from GKw:**
 * Setting γ=1, δ=0, λ=1 in the GKw PDF:
 * - The Beta function becomes: \eqn{B(1,0+1) = B(1,1) = 1}
 * - The exponent simplifies: \eqn{\gamma\lambda - 1 = 1 \cdot 1 - 1 = 0}
 * - The final term vanishes: \eqn{\{...\}^\delta = \{...\}^0 = 1}
 * 
 * The cumulative distribution function is:
 * \deqn{
 *   F(x) = 1 - (1-x^\alpha)^\beta
 * }
 * 
 * The quantile function (inverse CDF) is:
 * \deqn{
 *   Q(p) = \left\{1 - (1-p)^{1/\beta}\right\}^{1/\alpha}
 * }
 * 
 * **Parameter Constraints:**
 * - \eqn{\alpha > 0} (shape parameter)
 * - \eqn{\beta > 0} (shape parameter)
 * 
 * **Special Cases:**
 * | Distribution | Condition | Relation |
 * |--------------|-----------|----------|
 * | Uniform(0,1) | \eqn{\alpha = \beta = 1} | Kw(1, 1) |
 * | Power function | \eqn{\beta = 1} | Kw(α, 1) |
 * 
 * **Random Variate Generation:**
 * Uses inverse transform method:
 * 1. Generate \eqn{U \sim Uniform(0,1)}
 * 2. Return \eqn{X = Q(U) = \{1 - (1-U)^{1/\beta}\}^{1/\alpha}}
 * 
 * **Advantages over Beta:**
 * - Closed-form CDF: \eqn{F(x) = 1-(1-x^\alpha)^\beta}
 * - Simpler moments formulas
 * - Efficient random number generation
 * - No special functions required (no Beta function in CDF/quantile)
 * 
 * **Numerical Stability:**
 * All computations use log-space arithmetic and numerically stable helper
 * functions from utils.h to prevent overflow/underflow.
 * 
 * **Implemented Functions:**
 * - dkw(): Probability density function (PDF)
 * - pkw(): Cumulative distribution function (CDF)
 * - qkw(): Quantile function (inverse CDF)
 * - rkw(): Random variate generation
 * - llkw(): Negative log-likelihood for MLE
 * - grkw(): Gradient of negative log-likelihood
 * - hskw(): Hessian of negative log-likelihood
 * 
 * @author Lopes, J. E.
 * @date 2025-01-07
 * 
 * @references
 * Kumaraswamy, P. (1980). A generalized probability density function for
 * double-bounded random processes. Journal of Hydrology, 46(1-2), 79-88.
 * 
 * Jones, M. C. (2009). Kumaraswamy's distribution: A beta-type distribution
 * with some tractability advantages. Statistical Methodology, 6(1), 70-81.
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
 * @brief Probability Density Function of the Kumaraswamy Distribution
 * 
 * Computes the density (or log-density) for the Kumaraswamy distribution
 * at specified quantiles.
 * 
 * @param x Vector of quantiles (values in (0,1))
 * @param alpha Shape parameter vector (α > 0)
 * @param beta Shape parameter vector (β > 0)
 * @param log_prob If TRUE, returns log-density; otherwise returns density
 * 
 * @return NumericVector of density values (or log-density if log_prob=TRUE)
 * 
 * @details
 * The log-density is computed as:
 * \deqn{
 *   \log f(x) = \log(\alpha) + \log(\beta)
 *   + (\alpha-1)\log(x) + (\beta-1)\log(1-x^\alpha)
 * }
 * 
 * @note Exported as .dkw_cpp for internal package use
 */
// [[Rcpp::export(.dkw_cpp)]]
Rcpp::NumericVector dkw(
    const arma::vec& x,
    const Rcpp::NumericVector& alpha,
    const Rcpp::NumericVector& beta,
    bool log_prob = false
) {
  // Convert R vectors to Armadillo vectors
  arma::vec a_vec(alpha.begin(), alpha.size());
  arma::vec b_vec(beta.begin(), beta.size());
  
  // Determine output length for recycling
  size_t N = std::max({x.n_elem, a_vec.n_elem, b_vec.n_elem});
  
  // Initialize result with appropriate default
  arma::vec out(N);
  out.fill(log_prob ? R_NegInf : 0.0);
  
  for (size_t i = 0; i < N; i++) {
    // Extract recycled parameters
    double a = a_vec[i % a_vec.n_elem];
    double b = b_vec[i % b_vec.n_elem];
    double xx = x[i % x.n_elem];
    
    // Validate parameters
    if (!check_kw_pars(a, b)) {
      continue;
    }
    
    // Check support: x must be in (0, 1)
    if (xx <= 0.0 || xx >= 1.0 || !R_finite(xx)) {
      continue;
    }
    
    // ---- Log-space computation of density ----
    
    // Normalization constant: log(αβ)
    double la = safe_log(a);
    double lb = safe_log(b);
    
    // Compute log(x) and log(x^α)
    double lx = safe_log(xx);
    double log_xalpha = a * lx;  // log(x^α)
    
    // Compute log(1 - x^α) using stable log1mexp
    double log_v = log1mexp(log_xalpha);
    if (!R_finite(log_v)) {
      continue;
    }
    
    // Assemble log-density:
    // log(f) = log(α) + log(β) + (α-1)*log(x) + (β-1)*log(1-x^α)
    double log_pdf = la + lb + (a - 1.0) * lx + (b - 1.0) * log_v;
    
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
 * @brief Cumulative Distribution Function of the Kumaraswamy Distribution
 * 
 * Computes the cumulative probability for the Kumaraswamy distribution
 * at specified quantiles.
 * 
 * @param q Vector of quantiles
 * @param alpha Shape parameter vector (α > 0)
 * @param beta Shape parameter vector (β > 0)
 * @param lower_tail If TRUE, returns P(X ≤ q); otherwise P(X > q)
 * @param log_p If TRUE, returns log-probability
 * 
 * @return NumericVector of cumulative probabilities
 * 
 * @details
 * The CDF has a closed form (unlike the Beta distribution):
 * \deqn{F(x) = 1 - (1-x^\alpha)^\beta}
 * 
 * @note Exported as .pkw_cpp for internal package use
 */
// [[Rcpp::export(.pkw_cpp)]]
Rcpp::NumericVector pkw(
    const arma::vec& q,
    const Rcpp::NumericVector& alpha,
    const Rcpp::NumericVector& beta,
    bool lower_tail = true,
    bool log_p = false
) {
  // Convert R vectors to Armadillo vectors
  arma::vec a_vec(alpha.begin(), alpha.size());
  arma::vec b_vec(beta.begin(), beta.size());
  
  // Determine output length for recycling
  size_t N = std::max({q.n_elem, a_vec.n_elem, b_vec.n_elem});
  
  arma::vec out(N);
  
  for (size_t i = 0; i < N; i++) {
    // Extract recycled parameters
    double a = a_vec[i % a_vec.n_elem];
    double b = b_vec[i % b_vec.n_elem];
    double xx = q[i % q.n_elem];
    
    // Validate parameters
    if (!check_kw_pars(a, b)) {
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
    double xalpha = safe_pow(xx, a);
    
    // Step 2: (1 - x^α)^β
    double one_minus_xalpha_beta = safe_pow(1.0 - xalpha, b);
    
    // Step 3: F(x) = 1 - (1 - x^α)^β
    double tmp = 1.0 - one_minus_xalpha_beta;
    
    // Boundary checks after computation
    if (tmp <= 0.0) {
      double val0 = lower_tail ? 0.0 : 1.0;
      out(i) = log_p ? safe_log(val0) : val0;
      continue;
    }
    if (tmp >= 1.0) {
      double val1 = lower_tail ? 1.0 : 0.0;
      out(i) = log_p ? safe_log(val1) : val1;
      continue;
    }
    
    double val = tmp;
    
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
 * @brief Quantile Function (Inverse CDF) of the Kumaraswamy Distribution
 * 
 * Computes quantiles for the Kumaraswamy distribution
 * given probability values.
 * 
 * @param p Vector of probabilities (values in [0,1])
 * @param alpha Shape parameter vector (α > 0)
 * @param beta Shape parameter vector (β > 0)
 * @param lower_tail If TRUE, probabilities are P(X ≤ x); otherwise P(X > x)
 * @param log_p If TRUE, probabilities are given as log(p)
 * 
 * @return NumericVector of quantiles
 * 
 * @details
 * The quantile function has a closed form (unlike the Beta distribution):
 * \deqn{Q(p) = \left\{1 - (1-p)^{1/\beta}\right\}^{1/\alpha}}
 * 
 * @note Exported as .qkw_cpp for internal package use
 */
// [[Rcpp::export(.qkw_cpp)]]
Rcpp::NumericVector qkw(
    const arma::vec& p,
    const Rcpp::NumericVector& alpha,
    const Rcpp::NumericVector& beta,
    bool lower_tail = true,
    bool log_p = false
) {
  // Convert R vectors to Armadillo vectors
  arma::vec a_vec(alpha.begin(), alpha.size());
  arma::vec b_vec(beta.begin(), beta.size());
  
  // Determine output length for recycling
  size_t N = std::max({p.n_elem, a_vec.n_elem, b_vec.n_elem});
  
  arma::vec out(N);
  
  for (size_t i = 0; i < N; i++) {
    // Extract recycled parameters
    double a = a_vec[i % a_vec.n_elem];
    double b = b_vec[i % b_vec.n_elem];
    double pp = p[i % p.n_elem];
    
    // Validate parameters
    if (!check_kw_pars(a, b)) {
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
    
    // ---- Compute quantile via closed-form formula ----
    
    // Step 1: 1 - p
    double step1 = 1.0 - pp;
    step1 = std::max(0.0, step1);
    
    // Step 2: (1 - p)^(1/β)
    double step2 = safe_pow(step1, 1.0 / b);
    
    // Step 3: 1 - (1 - p)^(1/β)
    double step3 = 1.0 - step2;
    step3 = std::max(0.0, step3);
    
    // Step 4: {1 - (1 - p)^(1/β)}^(1/α)
    double xval;
    if (a == 1.0) {
      xval = step3;
    } else {
      xval = safe_pow(step3, 1.0 / a);
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
 * @brief Random Variate Generation for the Kumaraswamy Distribution
 * 
 * Generates random samples from the Kumaraswamy distribution
 * using the inverse transform method.
 * 
 * @param n Number of random variates to generate
 * @param alpha Shape parameter vector (α > 0)
 * @param beta Shape parameter vector (β > 0)
 * 
 * @return NumericVector of n random variates from Kw distribution
 * 
 * @details
 * Algorithm (extremely efficient due to closed-form quantile):
 * 1. Generate U ~ Uniform(0,1)
 * 2. Return X = Q(U) = {1 - (1-U)^(1/β)}^(1/α)
 * 
 * This is faster than Beta distribution random generation which requires
 * rejection sampling or more complex algorithms.
 * 
 * @note Exported as .rkw_cpp for internal package use
 */
// [[Rcpp::export(.rkw_cpp)]]
Rcpp::NumericVector rkw(
    int n,
    const Rcpp::NumericVector& alpha,
    const Rcpp::NumericVector& beta
) {
  if (n <= 0) {
    Rcpp::stop("rkw: n must be positive");
  }
  
  // Convert R vectors to Armadillo vectors
  arma::vec a_vec(alpha.begin(), alpha.size());
  arma::vec b_vec(beta.begin(), beta.size());
  
  arma::vec out(n);
  
  for (int i = 0; i < n; i++) {
    // Extract recycled parameters (direct modulo, no intermediate variable)
    double a = a_vec[i % a_vec.n_elem];
    double b = b_vec[i % b_vec.n_elem];
    
    // Validate parameters
    if (!check_kw_pars(a, b)) {
      out(i) = NA_REAL;
      Rcpp::warning("rkw: invalid parameters at index %d", i + 1);
      continue;
    }
    
    // Generate U ~ Uniform(0,1)
    double U = R::runif(0.0, 1.0);
    
    // Step 1: 1 - U
    double step1 = 1.0 - U;
    step1 = std::max(0.0, step1);
    
    // Step 2: (1 - U)^(1/β)
    double step2 = safe_pow(step1, 1.0 / b);
    
    // Step 3: 1 - (1 - U)^(1/β)
    double step3 = 1.0 - step2;
    step3 = std::max(0.0, step3);
    
    // Step 4: {1 - (1 - U)^(1/β)}^(1/α)
    double x;
    if (a == 1.0) {
      x = step3;
    } else {
      x = safe_pow(step3, 1.0 / a);
    }
    
    // Clamp to valid support
    x = std::max(0.0, std::min(1.0, x));
    
    out(i) = x;
  }
  
  return Rcpp::NumericVector(out.memptr(), out.memptr() + out.n_elem);
}


// ============================================================================
// NEGATIVE LOG-LIKELIHOOD FUNCTION
// ============================================================================

/**
 * @brief Negative Log-Likelihood for Kumaraswamy Distribution
 * 
 * Computes the negative log-likelihood function for parameter estimation
 * via maximum likelihood.
 * 
 * @param par Parameter vector of length 2: (α, β)
 * @param data Vector of observations (must be in (0,1))
 * 
 * @return Negative log-likelihood value (scalar)
 * 
 * @details
 * The log-likelihood for n observations is:
 * \deqn{
 *   \ell(\theta) = n[\ln\alpha + \ln\beta]
 *   + (\alpha-1)\sum\ln x_i + (\beta-1)\sum\ln(1-x_i^\alpha)
 * }
 * 
 * Returns +Inf for invalid parameters or data outside (0,1).
 * 
 * @note Exported as .llkw_cpp for internal package use
 */
// [[Rcpp::export(.llkw_cpp)]]
double llkw(const Rcpp::NumericVector& par, const Rcpp::NumericVector& data) {
  // Validate parameter vector length
  if (par.size() < 2) {
    return R_PosInf;
  }
  
  // Extract parameters
  double a = par[0];
  double b = par[1];
  
  // Validate parameters using consistent checker
  if (!check_kw_pars(a, b)) {
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
  
  // Constant term: n * [log(α) + log(β)]
  double cst = n * (safe_log(a) + safe_log(b));
  
  // Term 1: (α-1) * Σ log(x)
  arma::vec lx = vec_safe_log(x);
  double sum1 = (a - 1.0) * arma::sum(lx);
  
  // Term 2: (β-1) * Σ log(1-x^α)
  // Use log-space for numerical stability
  arma::vec log_xalpha = a * lx;
  arma::vec log_v = vec_log1mexp(log_xalpha);
  double sum2 = (b - 1.0) * arma::sum(log_v);
  
  double loglike = cst + sum1 + sum2;
  
  return -loglike;
}


// ============================================================================
// GRADIENT OF NEGATIVE LOG-LIKELIHOOD
// ============================================================================

/**
 * @brief Gradient of Negative Log-Likelihood for Kumaraswamy Distribution
 * 
 * Computes the gradient vector of the negative log-likelihood for
 * optimization-based parameter estimation.
 * 
 * @param par Parameter vector of length 2: (α, β)
 * @param data Vector of observations (must be in (0,1))
 * 
 * @return NumericVector of length 2 containing partial derivatives
 *         with respect to (α, β)
 * 
 * @details
 * The gradient components are:
 * - ∂ℓ/∂α = n/α + Σlog(x) - (β-1)Σ[x^α log(x)/(1-x^α)]
 * - ∂ℓ/∂β = n/β + Σlog(1-x^α)
 * 
 * @note Exported as .grkw_cpp for internal package use
 */
// [[Rcpp::export(.grkw_cpp)]]
Rcpp::NumericVector grkw(const Rcpp::NumericVector& par, const Rcpp::NumericVector& data) {
  // Validate parameter vector length
  if (par.size() < 2) {
    return Rcpp::NumericVector(2, R_NaN);
  }
  
  // Extract parameters
  double alpha = par[0];
  double beta = par[1];
  
  // Validate parameters using consistent checker
  if (!check_kw_pars(alpha, beta)) {
    return Rcpp::NumericVector(2, R_NaN);
  }
  
  // Convert and validate data
  arma::vec x = Rcpp::as<arma::vec>(data);
  if (arma::any(x <= 0) || arma::any(x >= 1)) {
    return Rcpp::NumericVector(2, R_NaN);
  }
  
  int n = x.n_elem;
  Rcpp::NumericVector grad(2, 0.0);
  
  // Numerical stability constant
  const double eps = std::numeric_limits<double>::epsilon() * 100;
  
  // ---- Compute intermediate quantities ----
  
  arma::vec log_x = vec_safe_log(x);
  arma::vec x_alpha = vec_safe_pow(x, alpha);
  arma::vec x_alpha_log_x = x_alpha % log_x;
  
  // v = 1 - x^α (with clamping for numerical stability)
  arma::vec v = 1.0 - x_alpha;
  v = arma::clamp(v, eps, 1.0 - eps);
  
  arma::vec log_v = vec_safe_log(v);
  
  // ---- Calculate gradient components ----
  
  // ∂ℓ/∂α = n/α + Σlog(x) - (β-1)Σ[x^α log(x)/(1-x^α)]
  double d_alpha = n / alpha + arma::sum(log_x);
  arma::vec alpha_term = (beta - 1.0) * x_alpha_log_x / v;
  d_alpha -= arma::sum(alpha_term);
  
  // ∂ℓ/∂β = n/β + Σlog(1-x^α)
  double d_beta = n / beta + arma::sum(log_v);
  
  // Return NEGATIVE gradient (for minimization of negative log-likelihood)
  grad[0] = -d_alpha;
  grad[1] = -d_beta;
  
  return grad;
}


// ============================================================================
// HESSIAN OF NEGATIVE LOG-LIKELIHOOD
// ============================================================================

/**
 * @brief Hessian Matrix of Negative Log-Likelihood for Kumaraswamy Distribution
 * 
 * Computes the Hessian matrix (matrix of second partial derivatives) of
 * the negative log-likelihood for standard error estimation and
 * optimization algorithms.
 * 
 * @param par Parameter vector of length 2: (α, β)
 * @param data Vector of observations (must be in (0,1))
 * 
 * @return NumericMatrix of dimension 2×2 containing the Hessian
 * 
 * @details
 * Computes analytical second derivatives. The Hessian is symmetric.
 * Parameter ordering: (α, β) → indices (0, 1).
 * 
 * The Hessian components are:
 * - H[α,α] = -n/α² - (β-1)Σ[x^α(log x)²(1+x^α/(1-x^α))/(1-x^α)]
 * - H[α,β] = -Σ[x^α log(x)/(1-x^α)]
 * - H[β,β] = -n/β²
 * 
 * Returns NaN matrix for invalid inputs.
 * 
 * @note Exported as .hskw_cpp for internal package use
 */
// [[Rcpp::export(.hskw_cpp)]]
Rcpp::NumericMatrix hskw(const Rcpp::NumericVector& par, const Rcpp::NumericVector& data) {
  // Initialize NaN matrix for error cases
  Rcpp::NumericMatrix nanHess(2, 2);
  nanHess.fill(R_NaN);
  
  // Validate parameter vector length
  if (par.size() < 2) {
    return nanHess;
  }
  
  // Extract parameters
  double alpha = par[0];
  double beta = par[1];
  
  // Validate parameters using consistent checker
  if (!check_kw_pars(alpha, beta)) {
    return nanHess;
  }
  
  // Convert and validate data
  arma::vec x = Rcpp::as<arma::vec>(data);
  if (arma::any(x <= 0) || arma::any(x >= 1)) {
    return nanHess;
  }
  
  int n = x.n_elem;
  Rcpp::NumericMatrix hess(2, 2);
  
  // Numerical stability constant
  const double eps = std::numeric_limits<double>::epsilon() * 100;
  
  // ---- Compute intermediate quantities ----
  
  arma::vec log_x = vec_safe_log(x);
  arma::vec log_x_squared = arma::square(log_x);
  arma::vec x_alpha = vec_safe_pow(x, alpha);
  arma::vec x_alpha_log_x = x_alpha % log_x;
  
  // v = 1 - x^α (with clamping for numerical stability)
  arma::vec v = 1.0 - x_alpha;
  v = arma::clamp(v, eps, 1.0 - eps);
  
  // Additional terms for second derivatives
  arma::vec term_ratio = x_alpha / v;              // x^α / (1-x^α)
  arma::vec term_combined = 1.0 + term_ratio;      // 1 + x^α/(1-x^α) = 1/(1-x^α)
  
  // ---- Calculate Hessian components (of log-likelihood ℓ) ----
  
  // H[α,α] = ∂²ℓ/∂α² = -n/α² - (β-1)Σ[x^α(log x)²(1+x^α/(1-x^α))/(1-x^α)]
  double h_alpha_alpha = -n / (alpha * alpha);
  arma::vec d2a_terms = (beta - 1.0) * x_alpha % log_x_squared % term_combined / v;
  h_alpha_alpha -= arma::sum(d2a_terms);
  
  // H[α,β] = H[β,α] = ∂²ℓ/∂α∂β = -Σ[x^α log(x)/(1-x^α)]
  double h_alpha_beta = -arma::sum(x_alpha_log_x / v);
  
  // H[β,β] = ∂²ℓ/∂β² = -n/β²
  double h_beta_beta = -n / (beta * beta);
  
  // Fill the Hessian matrix (symmetric) - NEGATE for negative log-likelihood
  hess(0, 0) = -h_alpha_alpha;
  hess(0, 1) = hess(1, 0) = -h_alpha_beta;
  hess(1, 1) = -h_beta_beta;
  
  return hess;
}
