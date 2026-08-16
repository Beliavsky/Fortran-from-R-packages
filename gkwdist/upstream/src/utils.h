// utils.h
// Utility functions for Generalized Kumaraswamy distributions package
// Author: Lopes, J. E.
// Date: 2025-10-07
//
// This header provides numerically stable functions and parameter validators
// for various Kumaraswamy-based distributions implemented in the package.

#ifndef GKWDIST_UTILS_H
#define GKWDIST_UTILS_H

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <limits>
#include <algorithm>
#include <string>
#include <functional>
#include <vector>
#include <random>

/*
 * ===========================================================================
 * COMPILE-TIME MATHEMATICAL CONSTANTS
 * ===========================================================================
 * All constants are computed at compile-time (constexpr) to eliminate
 * runtime initialization overhead. Values are provided with maximum
 * precision available in IEEE 754 double precision (53-bit mantissa).
 */

namespace {  // Anonymous namespace prevents ODR violations across translation units

// Machine precision constants
constexpr double EPSILON      = std::numeric_limits<double>::epsilon();        // ~2.220446e-16
  constexpr double DBL_MIN_SAFE = std::numeric_limits<double>::min() * 10.0;   // ~2.225074e-307
  constexpr double DBL_MAX_SAFE = std::numeric_limits<double>::max() / 10.0;   // ~1.797693e+307
  
  // Logarithmic bounds (pre-calculated for performance)
  // NOTE: LOG_DBL_MIN = log(DBL_MIN) = log(numeric_limits<double>::min()), NOT log(DBL_MIN_SAFE)
  constexpr double LOG_DBL_MIN  = -708.3964185322641;  // log(std::numeric_limits<double>::min())
  constexpr double LOG_DBL_MAX  = 308.2547155599167;   // log(DBL_MAX_SAFE)
  
  // Mathematical constants (maximum precision)
  constexpr double LN2          = 0.6931471805599453094172321214581766;  // log(2)
  constexpr double SQRT_EPSILON = 1.4901161193847656e-08;  // sqrt(EPSILON) for double
  
  // Optimized thresholds for numerical stability (based on Mächler 2012)
  constexpr double LOG1MEXP_CROSSOVER = -0.6931471805599453;  // -log(2)
  constexpr double LOG1MEXP_TINY      = -1.0e-14;  // Threshold for Taylor expansion
  constexpr double LOG1PEXP_LOWER     = -37.0;     // Below: exp(x) alone suffices
  constexpr double LOG1PEXP_MEDIUM    = 18.0;      // Transition to log1p formulation
  constexpr double LOG1PEXP_UPPER     = 33.3;      // Transition to x + exp(-x)
  constexpr double LOG1PEXP_LARGE     = 700.0;     // Above: x alone suffices
  
  // Parameter validation bounds (strict mode)
  constexpr double STRICT_MIN_PARAM   = 1.0e-8;    // Minimum parameter value (strict)
  constexpr double STRICT_MAX_PARAM   = 1.0e8;     // Maximum parameter value (strict)
  
  // Exponent limits for safe_pow overflow detection
  constexpr double EXTREME_EXPONENT   = 1.0e10;    // Threshold for special exponent handling
  
  // Tolerance for integer detection in negative base powers
  constexpr double INTEGER_TOLERANCE  = 1.0e-12;   // Tolerance for y ≈ round(y)
  
} // namespace

/*
 * ===========================================================================
 * CORE NUMERICAL STABILITY FUNCTIONS
 * ===========================================================================
 * These functions implement numerically stable computations for logarithmic
 * operations that are prone to catastrophic cancellation or overflow.
 * Implementations follow best practices from:
 * - Mächler, M. (2012). "Accurately Computing log(1-exp(-|a|))"
 * - R Core Team. R Mathlib implementation
 */

/**
 * log1mexp: Compute log(1 - exp(u)) with enhanced numerical stability
 * 
 * For u <= 0, computes log(1 - exp(u)) using different approximations
 * depending on the magnitude of u to avoid catastrophic cancellation.
 * 
 * Method selection (Mächler 2012):
 *   u > -1e-14        : Taylor expansion log(-u) + corrections
 *   -log(2) < u <= 0  : log(-expm1(u))
 *   u <= -log(2)      : log1p(-exp(u))
 * 
 * @param u Non-positive value (u <= 0)
 * @return log(1 - exp(u)), or NaN if u > 0
 * 
 * @note Time complexity: O(1)
 * @note Relative error: < 2*EPSILON for all u <= 0
 */
inline double log1mexp(double u) {
  // Input validation: u must be non-positive
  if (u > 0.0) {
    return R_NaN;
  }
  
  // Region 1: Very small |u| - use Taylor series with correction
  // For u ≈ 0⁻, 1 - exp(u) ≈ -u - u²/2 + O(u³)
  // Therefore log(1 - exp(u)) ≈ log(-u) + log(1 + u/2) ≈ log(-u) - u/2
  if (u > LOG1MEXP_TINY) {
    double neg_u = -u;
    // Second-order correction for improved accuracy
    return std::log(neg_u) - 0.5 * u;
  }
  
  // Region 2: -log(2) < u <= -1e-14 - use expm1 formulation
  // expm1(u) = exp(u) - 1, so -expm1(u) = 1 - exp(u)
  if (u > LOG1MEXP_CROSSOVER) {
    return std::log(-std::expm1(u));
  }
  
  // Region 3: u <= -log(2) - use log1p formulation
  // Most numerically stable for large |u|
  return std::log1p(-std::exp(u));
}

/**
 * log1pexp: Compute log(1 + exp(x)) with protection against overflow
 * 
 * Handles various regimes of x with appropriate approximations to maintain
 * numerical stability across the entire real line.
 * 
 * Method selection:
 *   x > 700     : x (overflow protection)
 *   x > 33.3    : x + exp(-x) (asymptotic expansion)
 *   x > 18      : x + log1p(exp(-x)) (better for moderate x)
 *   x > -37     : log1p(exp(x)) (standard range)
 *   x > -700    : exp(x) (exp(x) << 1)
 *   x <= -700   : 0 (complete underflow)
 * 
 * @param x Input value (unrestricted)
 * @return log(1 + exp(x)) calculated with numerical stability
 * 
 * @note Time complexity: O(1)
 * @note Relative error: < 2*EPSILON for all x
 */
inline double log1pexp(double x) {
  // Region 1: Very large x - asymptotic to x
  if (x > LOG1PEXP_LARGE) {
    return x;
  }
  
  // Region 2: Large x - first-order correction
  if (x > LOG1PEXP_UPPER) {
    return x + std::exp(-x);
  }
  
  // Region 3: Moderately large x - use log1p with negative exponent
  if (x > LOG1PEXP_MEDIUM) {
    return x + std::log1p(std::exp(-x));
  }
  
  // Region 4: Standard range - direct log1p
  if (x > LOG1PEXP_LOWER) {
    return std::log1p(std::exp(x));
  }
  
  // Region 5: Large negative x - exp(x) dominates
  if (x > -LOG1PEXP_LARGE) {
    return std::exp(x);
  }
  
  // Region 6: Extreme negative x - complete underflow
  return 0.0;
}

/**
 * safe_log: Compute log(x) with comprehensive error handling
 * 
 * @param x Input value
 * @return log(x), -Inf for x=0, NaN for x<0, or scaled result for tiny x
 * 
 * @note Handles underflow gracefully for very small positive x
 * @note Time complexity: O(1)
 */
inline double safe_log(double x) {
  // Handle invalid inputs
  if (x < 0.0) {
    return R_NaN;
  }
  
  if (x == 0.0) {
    return R_NegInf;
  }
  
  // Handle potential underflow with scaled computation
  // For x = ε * DBL_MIN_SAFE where ε is small,
  // log(x) = log(ε) + log(DBL_MIN_SAFE) = log(ε) + LOG_DBL_MIN
  if (x < DBL_MIN_SAFE) {
    return LOG_DBL_MIN + std::log(x / DBL_MIN_SAFE);
  }
  
  return std::log(x);
}

/**
 * safe_exp: Compute exp(x) with protection against overflow/underflow
 * 
 * @param x Input value
 * @return exp(x), +Inf for overflow, 0 or scaled result for underflow
 * 
 * @note Uses scaled arithmetic near underflow threshold for gradual transition
 * @note Time complexity: O(1)
 */
inline double safe_exp(double x) {
  // Handle overflow
  if (x > LOG_DBL_MAX) {
    return R_PosInf;
  }
  
  // Handle severe underflow
  if (x < LOG_DBL_MIN - 10.0) {
    return 0.0;
  }
  
  // Handle moderate underflow with scaling:
  // DBL_MIN * exp(x - log(DBL_MIN)) = exp(x)  — exact, avoids abrupt flush to zero
  if (x < LOG_DBL_MIN) {
    return std::numeric_limits<double>::min() * std::exp(x - LOG_DBL_MIN);
  }
  
  return std::exp(x);
}

/**
 * safe_pow: Compute x^y with robust error handling and numerical stability
 * 
 * Handles special cases comprehensively:
 * - x = 0: Returns 0 (y>0), 1 (y=0), +Inf (y<0)
 * - x = 1 or y = 0: Returns 1
 * - y = 1: Returns x
 * - x < 0: Requires y to be effectively integer; handles sign correctly
 * - Extreme exponents: Prevents overflow/underflow with early detection
 * 
 * For positive x, uses logarithmic transformation: x^y = exp(y * log(x))
 * This provides better numerical stability than direct pow() for extreme values.
 * 
 * @param x Base value
 * @param y Exponent value
 * @return x^y calculated with numerical stability and comprehensive edge case handling
 * 
 * @note Time complexity: O(1)
 * @note For x < 0, y must satisfy |y - round(y)| < INTEGER_TOLERANCE
 */
inline double safe_pow(double x, double y) {
  // Handle NaN propagation
  if (std::isnan(x) || std::isnan(y)) {
    return R_NaN;
  }
  
  // ===== Handle x = 0 cases =====
  if (x == 0.0) {
    if (y > 0.0)  return 0.0;        // 0^positive = 0
    if (y == 0.0) return 1.0;        // 0^0 = 1 (standard convention in probability)
    return R_PosInf;                  // 0^negative = +Inf
  }
  
  // ===== Trivial cases =====
  if (x == 1.0 || y == 0.0) return 1.0;  // 1^y = 1, x^0 = 1
  if (y == 1.0) return x;                // x^1 = x
  
  // ===== Handle negative base =====
  if (x < 0.0) {
    // Check if y is effectively an integer
    double y_rounded = std::round(y);
    if (std::abs(y - y_rounded) > INTEGER_TOLERANCE) {
      return R_NaN;  // Non-integer power of negative number is undefined in reals
    }
    
    // y is integer - compute |x|^|y| then apply sign
    int y_int = static_cast<int>(y_rounded);
    bool y_is_odd = (y_int % 2 != 0);
    double abs_x = -x;  // x is negative, so -x is positive
    
    // Compute absolute result using logarithmic method for stability
    double log_abs_x = std::log(abs_x);
    double log_result = std::abs(y) * log_abs_x;
    
    // Check for overflow/underflow
    if (log_result > LOG_DBL_MAX) {
      return y_is_odd ? R_NegInf : R_PosInf;
    }
    if (log_result < LOG_DBL_MIN) {
      return 0.0;
    }
    
    double abs_result = std::exp(log_result);
    
    // Apply sign based on whether y is odd and whether we're inverting
    if (y < 0) {
      // Negative exponent: invert result
      if (abs_result == 0.0) return y_is_odd ? R_NegInf : R_PosInf;
      abs_result = 1.0 / abs_result;
    }
    
    return y_is_odd ? -abs_result : abs_result;
  }
  
  // ===== Positive base: use logarithmic transformation =====
  
  // For extreme exponents, check bounds before computation
  if (std::abs(y) > EXTREME_EXPONENT) {
    double log_x = std::log(x);
    double log_result = y * log_x;
    
    // Early overflow/underflow detection
    if (log_result > LOG_DBL_MAX) {
      return R_PosInf;
    }
    if (log_result < LOG_DBL_MIN) {
      return 0.0;
    }
    
    return std::exp(log_result);
  }
  
  // Standard case: compute via logarithm for better stability
  double log_x = std::log(x);
  double log_result = y * log_x;
  
  // Use safe_exp for final result
  return safe_exp(log_result);
}

/*
 * ===========================================================================
 * VECTORIZED NUMERICAL STABILITY FUNCTIONS
 * ===========================================================================
 * Element-wise operations on Armadillo vectors with optimized implementation.
 * These functions maintain numerical stability while leveraging SIMD when possible.
 */

/**
 * vec_log1mexp: Vectorized log(1 - exp(u)) computation
 * 
 * @param u Vector of non-positive values
 * @return Vector of log(1 - exp(u)) values
 * 
 * @note Time complexity: O(n)
 * @note Memory complexity: O(n)
 */
inline arma::vec vec_log1mexp(const arma::vec& u) {
  const size_t n = u.n_elem;
  arma::vec result(n);
  
  // Element-wise processing for maximum numerical reliability
  // Each element may fall in different numerical regime
  for (size_t i = 0; i < n; ++i) {
    result(i) = log1mexp(u(i));
  }
  
  return result;
}

/**
 * vec_log1pexp: Vectorized log(1 + exp(x)) computation
 * 
 * @param x Vector of input values
 * @return Vector of log(1 + exp(x)) values
 * 
 * @note Time complexity: O(n)
 * @note Memory complexity: O(n)
 */
inline arma::vec vec_log1pexp(const arma::vec& x) {
  const size_t n = x.n_elem;
  arma::vec result(n);
  
  for (size_t i = 0; i < n; ++i) {
    result(i) = log1pexp(x(i));
  }
  
  return result;
}

/**
 * vec_safe_log: Vectorized safe logarithm computation
 * 
 * @param x Vector of input values
 * @return Vector of safe_log(x) values
 * 
 * @note Time complexity: O(n)
 * @note Memory complexity: O(n)
 */
inline arma::vec vec_safe_log(const arma::vec& x) {
  const size_t n = x.n_elem;
  arma::vec result(n);
  
  for (size_t i = 0; i < n; ++i) {
    result(i) = safe_log(x(i));
  }
  
  return result;
}

/**
 * vec_safe_exp: Vectorized safe exponential computation
 * 
 * @param x Vector of input values
 * @return Vector of safe_exp(x) values
 * 
 * @note Time complexity: O(n)
 * @note Memory complexity: O(n)
 */
inline arma::vec vec_safe_exp(const arma::vec& x) {
  const size_t n = x.n_elem;
  arma::vec result(n);
  
  for (size_t i = 0; i < n; ++i) {
    result(i) = safe_exp(x(i));
  }
  
  return result;
}

/**
 * vec_safe_pow: Vectorized safe power computation (scalar exponent)
 * 
 * @param x Vector of base values
 * @param y Scalar exponent value
 * @return Vector of x[i]^y values
 * 
 * @note Time complexity: O(n)
 * @note Memory complexity: O(n)
 * @note Optimized for case where y is constant across all elements
 */
inline arma::vec vec_safe_pow(const arma::vec& x, double y) {
  const size_t n = x.n_elem;
  arma::vec result(n);
  
  // Optimize for common trivial cases (vectorizable)
  if (y == 0.0) {
    result.ones();
    return result;
  }
  
  if (y == 1.0) {
    return x;
  }

  // Fast SIMD path: y > 0, all x positive — typical distribution case (data in (0,1)).
  // arma::exp(y*arma::log(x)) is fully auto-vectorizable; the branch-heavy scalar loop below is not.
  if (y > 0.0 && y < EXTREME_EXPONENT && x.min() > 0.0) {
    return arma::exp(y * arma::log(x));
  }

  // Check if y is effectively an integer (for negative base handling)
  double y_rounded = std::round(y);
  bool y_is_integer = (std::abs(y - y_rounded) <= INTEGER_TOLERANCE);
  // Guard against UB: static_cast<int> is UB when y_rounded > INT_MAX
  bool y_fits_in_int = (std::abs(y_rounded) <= static_cast<double>(std::numeric_limits<int>::max()));
  int y_int = y_fits_in_int ? static_cast<int>(y_rounded) : 0;
  bool y_is_odd = y_is_integer && y_fits_in_int && (y_int % 2 != 0);
  
  // Element-wise computation with shared exponent logic
  for (size_t i = 0; i < n; ++i) {
    double xi = x(i);
    
    // Handle NaN
    if (std::isnan(xi)) {
      result(i) = R_NaN;
      continue;
    }
    
    // Handle xi = 0
    if (xi == 0.0) {
      if (y > 0.0) {
        result(i) = 0.0;
      } else if (y == 0.0) {
        result(i) = 1.0;
      } else {
        result(i) = R_PosInf;
      }
      continue;
    }
    
    // Handle xi = 1
    if (xi == 1.0) {
      result(i) = 1.0;
      continue;
    }
    
    // Handle negative base
    if (xi < 0.0) {
      if (!y_is_integer) {
        result(i) = R_NaN;
      } else {
        double abs_xi = -xi;
        double log_abs_xi = std::log(abs_xi);
        double log_result = std::abs(y) * log_abs_xi;
        
        if (log_result > LOG_DBL_MAX) {
          result(i) = y_is_odd ? R_NegInf : R_PosInf;
        } else if (log_result < LOG_DBL_MIN) {
          result(i) = 0.0;
        } else {
          double abs_result = std::exp(log_result);
          if (y < 0) abs_result = 1.0 / abs_result;
          result(i) = y_is_odd ? -abs_result : abs_result;
        }
      }
      continue;
    }
    
    // Positive base: logarithmic computation
    double log_xi = std::log(xi);
    double log_result = y * log_xi;
    result(i) = safe_exp(log_result);
  }
  
  return result;
}

/**
 * vec_safe_pow: Vectorized safe power computation (vector exponents)
 * 
 * @param x Vector of base values
 * @param y Vector of exponent values (must match size of x)
 * @return Vector of x[i]^y[i] values
 * 
 * @note Time complexity: O(n)
 * @note Memory complexity: O(n)
 */
inline arma::vec vec_safe_pow(const arma::vec& x, const arma::vec& y) {
  const size_t n = x.n_elem;
  
  // Input validation
  if (y.n_elem != n) {
    Rcpp::stop("vec_safe_pow: vectors must have same length (x: %d, y: %d)", n, y.n_elem);
  }
  
  arma::vec result(n);
  
  // Element-wise computation
  for (size_t i = 0; i < n; ++i) {
    result(i) = safe_pow(x(i), y(i));
  }
  
  return result;
}

/*
 * ===========================================================================
 * PARAMETER VALIDATION FUNCTIONS
 * ===========================================================================
 * These functions verify that distribution parameters satisfy required
 * constraints. Each distribution family has specific requirements.
 * 
 * The 'strict' parameter enables additional bounds checking to prevent
 * numerical instabilities that can arise from extreme parameter values.
 */

/**
 * check_pars: Validate parameters for Generalized Kumaraswamy (GKw) distribution
 * 
 * Parameter constraints:
 *   alpha > 0   (shape parameter)
 *   beta > 0    (shape parameter)
 *   gamma > 0   (shape parameter)
 *   delta >= 0  (shape parameter, allows zero)
 *   lambda > 0  (shape parameter)
 * 
 * Strict mode additionally enforces:
 *   All parameters in [1e-8, 1e8] to prevent numerical issues
 * 
 * @param alpha Shape parameter (must be > 0)
 * @param beta Shape parameter (must be > 0)
 * @param gamma Shape parameter (must be > 0)
 * @param delta Shape parameter (must be >= 0)
 * @param lambda Shape parameter (must be > 0)
 * @param strict Enable strict bounds checking
 * @return true if parameters are valid, false otherwise
 */
inline bool check_pars(double alpha,
                       double beta,
                       double gamma,
                       double delta,
                       double lambda,
                       bool strict = false) {
  // Check for NaN values
  if (std::isnan(alpha) || std::isnan(beta) || std::isnan(gamma) ||
      std::isnan(delta) || std::isnan(lambda)) {
    return false;
  }
  
  // Check for Inf values
  if (std::isinf(alpha) || std::isinf(beta) || std::isinf(gamma) ||
      std::isinf(delta) || std::isinf(lambda)) {
    return false;
  }
  
  // Basic parameter constraints
  if (alpha <= 0.0 || beta <= 0.0 || gamma <= 0.0 || delta < 0.0 || lambda <= 0.0) {
    return false;
  }
  
  // Strict bounds for numerical stability
  if (strict) {
    if (alpha < STRICT_MIN_PARAM || beta < STRICT_MIN_PARAM || 
        gamma < STRICT_MIN_PARAM || lambda < STRICT_MIN_PARAM) {
      return false;
    }
    
    if (alpha > STRICT_MAX_PARAM || beta > STRICT_MAX_PARAM || 
        gamma > STRICT_MAX_PARAM || delta > STRICT_MAX_PARAM || 
        lambda > STRICT_MAX_PARAM) {
      return false;
    }
    
    // Delta can be zero, but if non-zero must satisfy bounds
    if (delta > 0.0 && delta < STRICT_MIN_PARAM) {
      return false;
    }
  }
  
  return true;
}

/**
 * check_pars_vec: Vectorized parameter validation for GKw distribution
 * 
 * Validates all combinations of parameter values using R-style recycling.
 * If vectors have different lengths, shorter ones are recycled.
 * 
 * @param alpha Vector of alpha values
 * @param beta Vector of beta values
 * @param gamma Vector of gamma values
 * @param delta Vector of delta values
 * @param lambda Vector of lambda values
 * @param strict Enable strict bounds checking
 * @return Vector of boolean values (0/1) indicating parameter validity
 * 
 * @note Return type is arma::uvec for compatibility with Armadillo indexing
 */
inline arma::uvec check_pars_vec(const arma::vec& alpha,
                                 const arma::vec& beta,
                                 const arma::vec& gamma,
                                 const arma::vec& delta,
                                 const arma::vec& lambda,
                                 bool strict = false) {
  // Find maximum length for recycling
  size_t n = std::max({alpha.n_elem, beta.n_elem, gamma.n_elem,
                      delta.n_elem, lambda.n_elem});
  
  arma::uvec valid(n);
  
  // Check each combination with proper recycling
  for (size_t i = 0; i < n; ++i) {
    double a = alpha(i % alpha.n_elem);
    double b = beta(i % beta.n_elem);
    double g = gamma(i % gamma.n_elem);
    double d = delta(i % delta.n_elem);
    double l = lambda(i % lambda.n_elem);
    
    valid(i) = check_pars(a, b, g, d, l, strict) ? 1 : 0;
  }
  
  return valid;
}

/**
 * check_kkw_pars: Validate parameters for Kw-Kumaraswamy (kkw) distribution
 * 
 * kkw is GKw with gamma = 1: kkw(α, β, δ, λ) = GKw(α, β, 1, δ, λ)
 * 
 * Parameter constraints:
 *   alpha > 0
 *   beta > 0
 *   delta >= 0
 *   lambda > 0
 * 
 * @param alpha Shape parameter (must be > 0)
 * @param beta Shape parameter (must be > 0)
 * @param delta Shape parameter (must be >= 0)
 * @param lambda Shape parameter (must be > 0)
 * @param strict Enable strict bounds checking
 * @return true if parameters are valid, false otherwise
 */
inline bool check_kkw_pars(double alpha,
                           double beta,
                           double delta,
                           double lambda,
                           bool strict = false) {
  // Check for NaN/Inf
  if (std::isnan(alpha) || std::isnan(beta) || std::isnan(delta) || std::isnan(lambda)) {
    return false;
  }
  if (std::isinf(alpha) || std::isinf(beta) || std::isinf(delta) || std::isinf(lambda)) {
    return false;
  }
  
  // Basic constraints
  if (alpha <= 0.0 || beta <= 0.0 || delta < 0.0 || lambda <= 0.0) {
    return false;
  }
  
  // Strict bounds
  if (strict) {
    if (alpha < STRICT_MIN_PARAM || beta < STRICT_MIN_PARAM || lambda < STRICT_MIN_PARAM) {
      return false;
    }
    if (alpha > STRICT_MAX_PARAM || beta > STRICT_MAX_PARAM || 
        delta > STRICT_MAX_PARAM || lambda > STRICT_MAX_PARAM) {
      return false;
    }
    if (delta > 0.0 && delta < STRICT_MIN_PARAM) {
      return false;
    }
  }
  
  return true;
}

/**
 * check_bkw_pars: Validate parameters for Beta-Kumaraswamy (BKw) distribution
 * 
 * BKw is GKw with lambda = 1: BKw(α, β, γ, δ) = GKw(α, β, γ, δ, 1)
 * 
 * Parameter constraints:
 *   alpha > 0
 *   beta > 0
 *   gamma > 0
 *   delta >= 0
 * 
 * @param alpha Shape parameter (must be > 0)
 * @param beta Shape parameter (must be > 0)
 * @param gamma Shape parameter (must be > 0)
 * @param delta Shape parameter (must be >= 0)
 * @param strict Enable strict bounds checking
 * @return true if parameters are valid, false otherwise
 */
inline bool check_bkw_pars(double alpha,
                           double beta,
                           double gamma,
                           double delta,
                           bool strict = false) {
  // Check for NaN/Inf
  if (std::isnan(alpha) || std::isnan(beta) || std::isnan(gamma) || std::isnan(delta)) {
    return false;
  }
  if (std::isinf(alpha) || std::isinf(beta) || std::isinf(gamma) || std::isinf(delta)) {
    return false;
  }
  
  // Basic constraints
  if (alpha <= 0.0 || beta <= 0.0 || gamma <= 0.0 || delta < 0.0) {
    return false;
  }
  
  // Strict bounds
  if (strict) {
    if (alpha < STRICT_MIN_PARAM || beta < STRICT_MIN_PARAM || gamma < STRICT_MIN_PARAM) {
      return false;
    }
    if (alpha > STRICT_MAX_PARAM || beta > STRICT_MAX_PARAM || 
        gamma > STRICT_MAX_PARAM || delta > STRICT_MAX_PARAM) {
      return false;
    }
    if (delta > 0.0 && delta < STRICT_MIN_PARAM) {
      return false;
    }
  }
  
  return true;
}

/**
 * check_ekw_pars: Validate parameters for Exponentiated-Kumaraswamy (EKw) distribution
 * 
 * EKw is GKw with gamma = 1, delta = 0: EKw(α, β, λ) = GKw(α, β, 1, 0, λ)
 * 
 * Parameter constraints:
 *   alpha > 0
 *   beta > 0
 *   lambda > 0
 * 
 * @param alpha Shape parameter (must be > 0)
 * @param beta Shape parameter (must be > 0)
 * @param lambda Shape parameter (must be > 0)
 * @param strict Enable strict bounds checking
 * @return true if parameters are valid, false otherwise
 */
inline bool check_ekw_pars(double alpha, double beta, double lambda, bool strict = false) {
  // Check for NaN/Inf
  if (std::isnan(alpha) || std::isnan(beta) || std::isnan(lambda)) {
    return false;
  }
  if (std::isinf(alpha) || std::isinf(beta) || std::isinf(lambda)) {
    return false;
  }
  
  // Basic constraints
  if (alpha <= 0.0 || beta <= 0.0 || lambda <= 0.0) {
    return false;
  }
  
  // Strict bounds
  if (strict) {
    if (alpha < STRICT_MIN_PARAM || beta < STRICT_MIN_PARAM || lambda < STRICT_MIN_PARAM) {
      return false;
    }
    if (alpha > STRICT_MAX_PARAM || beta > STRICT_MAX_PARAM || lambda > STRICT_MAX_PARAM) {
      return false;
    }
  }
  
  return true;
}

/**
 * check_bp_pars: Validate parameters for Beta-Power (BP) distribution
 * 
 * BP is GKw with alpha = beta = 1: BP(γ, δ, λ) = GKw(1, 1, γ, δ, λ)
 * 
 * Parameter constraints:
 *   gamma > 0
 *   delta >= 0
 *   lambda > 0
 * 
 * @param gamma Shape parameter (must be > 0)
 * @param delta Shape parameter (must be >= 0)
 * @param lambda Shape parameter (must be > 0)
 * @param strict Enable strict bounds checking
 * @return true if parameters are valid, false otherwise
 */
inline bool check_bp_pars(double gamma, double delta, double lambda, bool strict = false) {
  // Check for NaN/Inf
  if (std::isnan(gamma) || std::isnan(delta) || std::isnan(lambda)) {
    return false;
  }
  if (std::isinf(gamma) || std::isinf(delta) || std::isinf(lambda)) {
    return false;
  }
  
  // Basic constraints
  if (gamma <= 0.0 || delta < 0.0 || lambda <= 0.0) {
    return false;
  }
  
  // Strict bounds
  if (strict) {
    if (gamma < STRICT_MIN_PARAM || lambda < STRICT_MIN_PARAM) {
      return false;
    }
    if (gamma > STRICT_MAX_PARAM || delta > STRICT_MAX_PARAM || lambda > STRICT_MAX_PARAM) {
      return false;
    }
    if (delta > 0.0 && delta < STRICT_MIN_PARAM) {
      return false;
    }
  }
  
  return true;
}

/**
 * check_kw_pars: Validate parameters for Kumaraswamy (Kw) distribution
 * 
 * Kw is GKw with gamma = delta = lambda = 1: Kw(α, β) = GKw(α, β, 1, 0, 1)
 * 
 * Parameter constraints:
 *   alpha > 0
 *   beta > 0
 * 
 * @param alpha Shape parameter (must be > 0)
 * @param beta Shape parameter (must be > 0)
 * @param strict Enable strict bounds checking
 * @return true if parameters are valid, false otherwise
 */
inline bool check_kw_pars(double alpha, double beta, bool strict = false) {
  // Check for NaN/Inf
  if (std::isnan(alpha) || std::isnan(beta)) {
    return false;
  }
  if (std::isinf(alpha) || std::isinf(beta)) {
    return false;
  }
  
  // Basic constraints
  if (alpha <= 0.0 || beta <= 0.0) {
    return false;
  }
  
  // Strict bounds
  if (strict) {
    if (alpha < STRICT_MIN_PARAM || beta < STRICT_MIN_PARAM) {
      return false;
    }
    if (alpha > STRICT_MAX_PARAM || beta > STRICT_MAX_PARAM) {
      return false;
    }
  }
  
  return true;
}

/**
 * check_beta_pars: Validate parameters for Beta distribution
 * 
 * Parameter constraints:
 *   gamma > 0  (shape1)
 *   delta > 0  (shape2, note: must be POSITIVE for Beta, unlike GKw where >= 0)
 * 
 * @param gamma Shape parameter 1 (must be > 0)
 * @param delta Shape parameter 2 (must be > 0)
 * @param strict Enable strict bounds checking
 * @return true if parameters are valid, false otherwise
 */
inline bool check_beta_pars(double gamma, double delta, bool strict = false) {
  // Check for NaN/Inf
  if (std::isnan(gamma) || std::isnan(delta)) {
    return false;
  }
  if (std::isinf(gamma) || std::isinf(delta)) {
    return false;
  }
  
  // Basic constraints (note: delta > 0 for Beta, not >= 0)
  if (gamma <= 0.0 || delta <= 0.0) {
    return false;
  }
  
  // Strict bounds
  if (strict) {
    if (gamma < STRICT_MIN_PARAM || delta < STRICT_MIN_PARAM) {
      return false;
    }
    if (gamma > STRICT_MAX_PARAM || delta > STRICT_MAX_PARAM) {
      return false;
    }
  }
  
  return true;
}

#endif // GKWDIST_UTILS_H
