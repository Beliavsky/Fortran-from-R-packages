# gkwdist 1.1.4

## CRAN Fix

* **`test-mle-performance.R`**: Added `skip_on_cran()` to all timing-based benchmark
  tests. These tests compare wall-clock times of analytical vs. numerical gradients and
  are inherently unreliable on shared/loaded CRAN check machines, causing spurious
  `ERROR` results. The tests remain available for local development.

# gkwdist 1.1.3

## Bug Fixes

* **`llgkw()` invalid parameter return** (`gkw.cpp`): Fixed critical error where the
  negative log-likelihood returned `R_NegInf` (−∞) for invalid parameters instead of
  `R_PosInf` (+∞). Gradient-based MLE optimizers interpret −∞ as a global minimum,
  causing them to converge to the invalid boundary rather than the true MLE.

* **`gkwinit.cpp` — delta validation** (`gkwinit.cpp`): Fixed internal `gkw_pdf()`
  rejecting `delta = 0` (a valid GKw parameter value) due to a strict `delta <= 0`
  check that should have been `delta < 0`.

* **`gkwinit.cpp` — EKw/Kw sub-family PDF mapping** (`gkwinit.cpp`): Fixed
  `ekw_pdf()` and `kw_pdf()` passing `delta = 1` instead of the correct `delta = 0`
  when delegating to `gkw_pdf()`. EKw and Kw are GKw sub-families with `delta = 0`,
  not `delta = 1`. This produced wrong starting values for MLE of these families.

* **`hsbkw()` — v^(β−1) computation** (`bkw.cpp`): Fixed the Hessian of the BKw
  negative log-likelihood returning a wrong value for β < 1. The ternary expression
  `(beta > 1.0) ? v_beta/v : 1.0` coincidentally produces the correct result for
  β = 1 but is wrong for all 0 < β < 1. Replaced with the exact formula
  `safe_exp((beta - 1.0) * ln_v)`.

## Numerical Stability

* **`safe_exp()` underflow scaling** (`utils.h`): Fixed a systematic 10× error in
  the moderate-underflow branch. The previous implementation used
  `DBL_MIN_SAFE * exp(x − log(DBL_MIN))` where `DBL_MIN_SAFE = 10 * DBL_MIN`,
  yielding `10 * exp(x)` instead of `exp(x)`. The fix uses
  `DBL_MIN * exp(x − log(DBL_MIN)) = exp(x)` exactly.

* **`dgkw()` silent boundary truncation removed** (`gkw.cpp`): Removed a block that
  silently skipped data points within `SQRT_EPSILON^(1/α)` of 0 or 1, returning
  density 0 for those points without warning. The log-space computation handles
  near-boundary values correctly without this truncation.

* **`llekw()` / `grekw()` — lambda clamping removed** (`ekw.cpp`): Removed the
  arbitrary cap `lambda_factor = min(lambda_factor, 1000)` applied to gradient and
  Hessian terms when λ > 1000. This distorted optimization for large-λ scenarios and
  produced incorrect standard errors.

## Code Quality

* **`gkwinit.cpp`**: Removed `using namespace Rcpp;` at file scope; replaced with
  explicit `Rcpp::` qualifications. Added NA/NaN filtering before moment computation
  to prevent silent corruption when input data contains missing values.

* **`bkw.cpp`**: Removed spurious `try/catch` blocks wrapping `Rcpp::as<arma::vec>()`
  conversions in `grbkw()` and `hsbkw()`. These conversions cannot throw in this
  context and the silent fallback masked type errors.

* **`gkw.cpp` / `ekw.cpp`**: Refactored Hessian accumulation to build only the upper
  triangle inside the observation loop and symmetrize once afterwards with
  `arma::symmatu()`, eliminating O(n × p²) redundant assignments.

* **`utils.h` — `vec_safe_pow()` UB guard**: Added guard preventing undefined
  behaviour when casting large `y_rounded` values (> `INT_MAX`) to `int` for
  odd-exponent sign detection.

* **`utils.h` — `vec_safe_pow()` SIMD fast path**: Added an early-return path
  `arma::exp(y * arma::log(x))` for the common case (y > 0, all x > 0) that
  is fully auto-vectorizable, improving throughput in gradient/Hessian evaluation.

# gkwdist 1.1.2

## Code Cleanup and Testing Enhancement

### C++ Code Cleanup

* **Removed legacy commented code**: Cleaned up all C++ source files (`gkw.cpp`, `bkw.cpp`, `kkw.cpp`, `ekw.cpp`, `kw.cpp`, `bpmc.cpp`, `beta_.cpp`) by removing old commented-out implementations that were kept for reference.
* **Code formatting**: Improved R wrapper formatting with consistent indentation and alignment in `.Call()` invocations and roxygen examples.

### New Test Suites

* **Analytical derivatives validation** (`test-derivatives-validation.R`):
  - 70 comprehensive tests validating gradient (`gr*`) and Hessian (`hs*`) functions
  - Compares analytical derivatives against numerical differentiation via `numDeriv`
  - Covers all 7 subfamilies: GKw, BKw, KKw, EKw, Mc, Kw, Beta
  - Multiple parameter configurations per subfamily for robustness

* **MLE performance benchmarks** (`test-mle-performance.R`):
  - Compares optimization efficiency across three scenarios: numerical-only, analytical gradient, and analytical gradient + Hessian
  - Validates that analytical derivatives provide equivalent or better accuracy
  - Tests convergence rates and computational time across all distribution families

### JOSS Paper

* **Added paper for JOSS submission** (`inst/paper/`):
  - Complete manuscript describing the package's statistical framework
  - Comprehensive bibliography with foundational references
  - Compiled PDF ready for submission

---

# gkwdist 1.1.1

## Major Refactoring Release

This release represents a comprehensive refactoring of the entire package codebase, focusing on numerical stability, code consistency, and maintainability.

### C++ Backend Overhaul

* **Unified utility functions**: Introduced `utils.h` header providing numerically stable implementations of critical functions:

  - `log1mexp()`: Stable computation of log(1 - exp(x)) using Mächler (2012) methodology
  - `log1pexp()`: Overflow-protected computation of log(1 + exp(x))
  - `safe_log()`, `safe_exp()`, `safe_pow()`: Protected arithmetic operations with graceful handling of edge cases
  - Vectorized versions (`vec_safe_log`, `vec_log1mexp`, etc.) for efficient array operations

* **Consistent parameter validation**: All distribution families now use dedicated parameter checkers (`check_pars()`, `check_kw_pars()`, `check_ekw_pars()`, etc.) that properly handle NaN, Inf, and boundary conditions.

* **Complete documentation**: All C++ source files now include comprehensive Doxygen-style documentation headers describing:
  - Mathematical formulas for PDF, CDF, quantile, and random generation
  - Parameter constraints and special cases
  - Numerical stability considerations
  - Relationship to parent GKw distribution

### Bug Fixes

* **Fixed critical bug in `qgkw()`**: Corrected logic error where `lower_tail` transformation was incorrectly applied when `log_p = TRUE`. The probability is now properly converted to linear scale before tail adjustment.

* **Fixed gradient calculation in `grkkw()`**: Resolved issue where `log_z` was not recomputed after clamping `z` to minimum threshold, causing corrupted gradient values near boundaries.

* **Fixed Hessian calculation in `hsmc()`**: Corrected sign errors and formula for the lambda component of the Hessian matrix for the Beta-Power/McDonald distribution.

* **Fixed gradient signs in `grmc()`**: Ensured consistent computation of log-likelihood gradient before negation for optimization.

### Code Quality Improvements

* **Eliminated unused variables**: Removed declared but unused constants (`exp_threshold`) and intermediate variables across all distribution files.

* **Removed redundant calculations**: Streamlined computations, notably in `pgkw()` where logarithm was computed twice for the same quantity.

* **Simplified parameter recycling**: Replaced double-modulo indexing pattern (`idx = i % k; vec[idx % vec.n_elem]`) with direct single-modulo access (`vec[i % vec.n_elem]`) in random generation functions.

* **Standardized function signatures**: All distribution functions now follow consistent patterns for parameter order, validation, and return value handling.

### R Wrapper Layer

* **Complete separation of R and C++ interfaces**: All exported R functions now serve as wrappers around internal C++ implementations (`.dgkw_cpp`, `.pgkw_cpp`, etc.), providing:
  - Enhanced input validation with informative error messages
  - Consistent argument checking across all distribution families
  - Proper NA/NaN propagation
  - Documentation accessible via standard R help system

### Distribution Families

All seven distribution families have been refactored with identical improvements:

| Distribution | Parameters | File |
|--------------|------------|------|
| Generalized Kumaraswamy (GKw) | α, β, γ, δ, λ | `gkw.cpp` |
| Kumaraswamy-Kumaraswamy (KKw) | α, β, δ, λ | `kkw.cpp` |
| Beta-Kumaraswamy (BKw) | α, β, γ, δ | `bkw.cpp` |
| Exponentiated Kumaraswamy (EKw) | α, β, λ | `ekw.cpp` |
| Beta-Power/McDonald (BP/Mc) | γ, δ, λ | `bpmc.cpp` |
| Kumaraswamy (Kw) | α, β | `kw.cpp` |
| Beta (GKw-style) | γ, δ | `beta.cpp` |

Each family includes: density (`d*`), distribution (`p*`), quantile (`q*`), random generation (`r*`), negative log-likelihood (`ll*`), gradient (`gr*`), and Hessian (`hs*`) functions.

### Technical Notes

* Minimum supported R version remains 3.5.0
* C++11 standard required (enabled via Rcpp plugin)
* Depends on RcppArmadillo for efficient linear algebra operations

### Acknowledgments

Special thanks to the thorough code review process that identified subtle numerical issues in edge cases, particularly for extreme parameter values and observations near distribution boundaries.

# gkwdist 1.0.7

# gkwdist 1.0.5

## Documentation Improvements

* **Enhanced Examples for Likelihood Functions**: All `ll*`, `gr*`, and `hs*` functions now include comprehensive examples demonstrating:
  - Maximum likelihood estimation with analytical gradients
  - Univariate profile likelihoods with confidence thresholds
  - 2D likelihood surfaces with confidence regions (90%, 95%, 99%)
  - Confidence ellipses with marginal intervals for parameter pairs
  - Numerical vs analytical derivative verification
  - Likelihood ratio tests and score tests

* **Professional Visualization Standards**: 
  - Consistent color scheme across all examples
  - Grid-adaptive algorithms for computational efficiency
  - Base R only - no external dependencies required

* **Complete Coverage**: Enhanced documentation for all distribution families (Kw, EKw, KKw, GKw) covering 2 to 5 parameters

* **Theoretical References**: Documentation cites foundational work by Carrasco et al. (2010), Jones (2009), Kumaraswamy (1980), and standard inference theory from Casella & Berger (2002)


# gkwdist 1.0.3
* **README.md**: Fix typos and faill link
  - Fix zzz.R file by removing useless texts

# gkwdist 1.0.2

# gkwdist 1.0.1

## Major Improvements

### Enhanced `gkwgetstartvalues()` Function
* **NEW**: Added `family` parameter to support all distribution families
  - Automatically returns correct number of parameters for each family
  - Family-specific initial value strategies for better convergence
  - Supported families: `"gkw"`, `"bkw"`, `"kkw"`, `"ekw"`, `"mc"`, `"kw"`, `"beta"`
  - Case-insensitive family names for user convenience

### Documentation Enhancements
* **README.md**: Complete rewrite with mathematical rigor
  - All LaTeX formulas corrected and verified for proper rendering
  - Eight comprehensive examples using `optim()` with analytical gradients
  - Corrected function signatures: all `ll*()`, `gr*()`, and `hs*()` functions use `(par, data)` signature
  - Added performance benchmarks demonstrating 10-50× speedup with C++ implementation
  - Hierarchical structure diagram for all distribution families
  - Model selection workflow and practical guidelines
  - Removed all references to deprecated `gkwfit()` function

### CRAN Submission Readiness
* **DESCRIPTION**: Fixed to meet CRAN requirements
  - Proper `Authors@R` field formatting
  - Removed unused dependencies (`numDeriv`)
  - Corrected package dependencies (`RcppArmadillo` only in `LinkingTo`)
  - Enhanced description with DOI references
  - Fixed maintainer email formatting

## Bug Fixes

* Fixed function call signatures in all README examples to match actual implementation
* Corrected parameter passing in optimization examples (now consistently use `(par, data)`)
* Fixed LaTeX rendering issues with `\left`/`\right` delimiters in GitHub Markdown

## Testing

* **NEW**: Comprehensive test suite using `testthat`
  - 100+ tests covering all exported functions
  - Tests for all 7 distribution families (GKw, BKw, KKw, EKw, MC, Kw, Beta)
  - PDF, CDF, quantile, and random generation tests
  - Log-likelihood, gradient, and Hessian validation
  - Parameter recovery tests with MLE
  - Edge cases and boundary condition handling
  - Integration tests for PDF-CDF consistency

## Performance

* All functions implemented in C++ for maximum computational efficiency
* Analytical derivatives (gradient and Hessian) provide exact computations
* Optimized numerical stability for extreme parameter values

## Notes

* This is the initial CRAN submission
* Package focuses exclusively on distribution functions (no high-level fitting interface)
* Companion package `gkwreg` provides regression modeling capabilities
* All user-facing functions maintain backward compatibility
* C++ implementation uses RcppArmadillo for linear algebra operations
* Analytical functions use robust log-scale computations to prevent overflow/underflow
* Random generation uses inverse CDF method where closed-form solutions exist

# gkwdist 0.1.0

## New Features

* Initial CRAN release
* Generalized Kumaraswamy distribution (5 parameters)
* Six nested sub-families: Beta, Kumaraswamy, Exponentiated-Kumaraswamy, 
  Kumaraswamy-Kumaraswamy, Beta-Kumaraswamy, and McDonald distributions
* Complete set of distribution functions (d/p/q/r)
* Log-likelihood, gradient, and Hessian functions for all families

## Performance

* Optimized C++ implementation via Rcpp
* Vectorized operations for speed
