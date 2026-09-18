#' @include continuous.R discrete.R configure.R
NULL

#' Create and Use Random Number Generators
#'
#' To create a random number generator (RNG) object use `randompack_rng()`,
#' and to specify the underlying RNG engine use `randompack_rng(engine)` where
#' `engine` is a character string naming the engine (see Available Engines
#' below).
#'
#' Once created, the RNG object provides methods for drawing samples from
#' various distributions (e.g., `$normal()`, `$uniform()`, `$int()`). The object
#' can be configured using configure methods (e.g., `$seed()`, `$randomize()`).
#' Multiple independent RNG objects can be used for parallel random number
#' generation across different processes or threads.
#'
#' @param engine RNG engine
#' @param bitexact Logical; set TRUE to make samples bit-identical across platforms
#' @param full_mantissa Logical; set TRUE to use full mantissas for floating draws
#'
#' @return An RNG object with methods for drawing random variates.
#'
#' @section Available Engines:
#' \tabular{lll}{
#'   \code{x256++simd} \tab\tab xorshift256++ with streams (default) \cr
#'   \code{sfc64simd} \tab\tab sfc64 with streams \cr
#'   \code{x256++} \tab\tab xoshiro256++ (Vigna and Blackman, 2018) \cr
#'   \code{x256**} \tab\tab xoshiro256** (Vigna and Blackman, 2018) \cr
#'   \code{x128+} \tab\tab xorshift128+ (Vigna, 2014) \cr
#'   \code{xoro++} \tab\tab xoroshiro128++ (Vigna and Blackman, 2016) \cr
#'   \code{pcg64} \tab\tab PCG64 DXSM (O'Neill, 2014) \cr
#'   \code{sfc64} \tab\tab sfc64 (Chris Doty-Humphrey, 2013) \cr
#'   \code{squares} \tab\tab squares64 (Widynski, 2021) \cr
#'   \code{philox} \tab\tab Philox-4×64 (Salmon and Moraes, 2011) \cr
#'   \code{cwg128} \tab\tab cwg128 (Działa, 2022) \cr
#'   \code{ranlux++} \tab\tab ranlux++ (Sibidanov, 2017) \cr
#'   \code{chacha20} \tab\tab ChaCha20 (Bernstein, 2008) \cr
#' }
#' 
#' @section Continuous distributions:
#' The RNG object provides methods for generating random variates from common
#' continuous probability distributions. All methods return a numeric vector
#' of length `len`.
#'
#' \describe{
#'   \item{`rng$unif(len)`}{Uniform variates on [0,1).}
#'   \item{`rng$unif(len, a, b)`}{Uniform variates on [a,b) with
#'     `a < b`.}
#'   \item{`rng$normal(len)`}{Standard normal variates (mean 0 and standard
#'     deviation 1).}
#'   \item{`rng$normal(len, mu, sigma)`}{Normal variates with mean `mu`
#'     and standard deviation `sigma` (defaults 0 and 1).}
#'   \item{`rng$skew_normal(len, mu, sigma, alpha)`}{Skew-normal variates with
#'     location `mu` and scale `sigma` (defaults 0 and 1), and shape `alpha`
#'     (required).}
#'   \item{`rng$lognormal(len, mu, sigma)`}{Lognormal variates derived from an
#'     underlying normal distribution (defaults 0 and 1).}
#'   \item{`rng$exp(len)`}{Standard exponential variates (scale 1).}
#'   \item{`rng$exp(len, scale)`}{Exponential variates with scale
#'     `scale`.}
#'   \item{`rng$gamma(len, shape, scale)`}{Gamma variates with given shape and
#'     scale (default 1).}
#'   \item{`rng$chi2(len, nu)`}{Chi-square variates with `nu` degrees of
#'     freedom.}
#'   \item{`rng$beta(len, a, b)`}{Beta variates with shape parameters `a`
#'     and `b`.}
#'   \item{`rng$t(len, nu)`}{Student's t variates with `nu` degrees of
#'     freedom.}
#'   \item{`rng$f(len, nu1, nu2)`}{F variates with `nu1` and `nu2`
#'     degrees of freedom.}
#'   \item{`rng$gumbel(len, mu, beta)`}{Gumbel variates with location `mu`
#'     and scale `beta` (defaults 0 and 1).}
#'   \item{`rng$pareto(len, xm, alpha)`}{Pareto variates with minimum value
#'     `xm` and shape `alpha`.}
#'   \item{`rng$weibull(len, shape, scale)`}{Weibull variates with given shape
#'     and scale (default 1).}
#'   \item{`rng$mvn(n, Sigma, mu = NULL)`}{Multivariate normal variates as an
#'     `n` by `d` matrix, where `d` is the dimension of
#'     `Sigma`.}
#' }
#'
#' @section Discrete distributions:
#' The RNG object provides methods for generating random variates from common
#' discrete distributions and combinatorial constructions.
#'
#' \describe{
#'   \item{`rng$int(len, min, max)`}{Uniform integers on `[min, max]`.}
#'   \item{`rng$perm(n)`}{Random permutation of `1:n`.}
#'   \item{`rng$sample(n, k)`}{Sample `k` elements without replacement
#'     from `1:n`.}
#'   \item{`rng$raw(len)`}{Generate `len` random bytes as a raw vector.}
#' }
#'
#' @section Configuration and Copying:
#' Methods for creating, configuring, and managing RNG state. All state-setting
#' methods accept numeric vectors (double or integer) whose elements must be
#' nonnegative whole numbers not exceeding \eqn{2^{32}-1}. Where applicable,
#' shorter vectors are padded with zeros.
#'
#' \describe{
#'   \item{`rng$seed(seed, spawn_key = integer(0))`}{
#'     Reinitialize the RNG deterministically from `seed` and an optional
#'     numeric vector `spawn_key`.
#'   }
#'   \item{`rng$randomize()`}{Randomize the RNG state from system entropy.}
#'   \item{`rng$jump(p)`}{
#'     Jump an xor-family or `ranlux++` engine ahead by \eqn{2^p} steps. The
#'     `x128+` and `xoro128++` engines support `p = 32, 64, 96`, while
#'     `x256++`, `x256**`, `x256++simd`, and `ranlux++` also support
#'     `p = 128` and `p = 192`.
#'   }
#'   \item{`rng$advance(delta)`}{
#'     Advance the `pcg64` engine by an arbitrary 128-bit delta. `delta` may be
#'     a scalar, interpreted as `[delta, 0]`, or a vector of up to 4 32-bit
#'     words; shorter vectors are zero-padded.
#'   }
#'   \item{`rng$duplicate()`}{Duplicate the RNG, preserving its state.}
#'   \item{`rng$serialize()`}{Serialize the current RNG state as a raw vector.}
#'   \item{`rng$deserialize(raw_state)`}{Restore state from a raw vector created
#'     by `serialize()`.}
#'   \item{`rng$set_state(state)`}{Set the engine state directly (advanced use).}
#'   \item{`rng$pcg64_set_inc(inc)`}{
#'     Set the increment of the PCG64 engine. The increment may have length up
#'     to 4 and shorter vectors are zero-padded.
#'   }
#'   \item{`rng$cwg128_set_weyl(weyl)`}{
#'     Set the Weyl increment of the cwg128 engine. `weyl` may have length up
#'     to 4 and shorter vectors are zero-padded.
#'   }
#'   \item{`rng$sfc64_set_abc(abc)`}{
#'     Set the `a`, `b`, `c` state words of the sfc64 engine. `abc` may have
#'     length up to 6 and shorter vectors are zero-padded.
#'   }
#'   \item{`rng$chacha_set_nonce(nonce)`}{
#'     Set the nonce of the ChaCha20 engine. The nonce may have length up to 3
#'     and shorter vectors are zero-padded.
#'   }
#'   \item{`rng$philox_set_key(key)`}{
#'     Set the key of the Philox engine. The key may have length up to 4 and
#'     shorter vectors are zero-padded.
#'   }
#'   \item{`rng$squares_set_key(key)`}{
#'     Set the key of the Squares engine. The key may have length up to 2 and
#'     shorter vectors are zero-padded.
#'   }
#' }
#'
#' @examples
#' # Create an RNG
#' rng <- randompack_rng()                    # Default engine (xoshiro256++)
#' rng_pcg <- randompack_rng("pcg64")         # Specify engine
#' rng_chacha <- randompack_rng("chacha20")
#'
#' # Continuous distributions
#' x <- rng$unif(5)
#' x <- rng$normal(100)                    # Standard normal
#' x <- rng$normal(100, 1, 2)              # N(1,2)
#' x <- rng$skew_normal(100, mu=0, sigma=1, alpha=2)
#' x <- rng$lognormal(5)
#' x <- rng$beta(5, a=2, b=3)
#' Sigma <- diag(2)
#' x <- rng$mvn(10, Sigma, mu=c(1,2))
#'
#' # Discrete distributions
#' x <- rng$int(5, min=1, max=10)
#' x <- rng$perm(5)
#' x <- rng$sample(10, k=3)
#' x <- rng$raw(4)
#'
#' # Configuration and copying
#' rng$seed(12345)                          # seed for reproducibility
#' rng$randomize()                          # randomize from system entropy
#' rngm <- randompack_rng(full_mantissa = TRUE)  # 53-bit mantissas for doubles
#' rng2 <- rng$duplicate()                  # duplicate with same state
#' identical(rng$unif(3), rng2$unif(3))     # TRUE
#' raw_state <- rng$serialize()             # save state
#' rng3 <- randompack_rng()                 # another default RNG
#' rng3$deserialize(raw_state)              # restore state
#' identical(rng$unif(3), rng3$unif(3))     # TRUE
#' rng_sq <- randompack_rng("squares")      # engine-specific state setting
#' rng_sq$set_state(c(2, 0, 0, 0))          # counter = (2,0)
#' rng_sq$squares_set_key(c(3,4))
#'
#' @seealso \code{\link{randompack_engines}} to list all available engines
#'
#' @export
randompack_rng <- function(engine = "x256++simd", bitexact = FALSE,
                           full_mantissa = FALSE) {
  RandompackRNG$new(engine = engine, bitexact = bitexact,
                    full_mantissa = full_mantissa)
}

#' Available RNG Engines
#'
#' Returns a data frame of supported random number generator engines with 
#' their descriptions.
#'
#' @return A data.frame with columns \code{engine} (short name) and 
#'   \code{description} (full name and citation).
#' 
#' @export
randompack_engines = function() {
  out <- .Call("randompack_engines_R", PACKAGE = "randompack")
  if (is.data.frame(out)) return(out)
  if (is.list(out) && length(out) == 2L) {
    return(data.frame(engine = out[[1]], description = out[[2]],
                      stringsAsFactors = FALSE))
  }
  stop("randompack_engines_R returned an unexpected value")
}

RandompackRNG <- R6::R6Class(
  "RandompackRNG",
  public = c(
    methods_configure,
    methods_discrete,
    methods_continuous,
    list(
      ptr = NULL,
      engine = NULL,
      initialize = function(engine = "", bitexact = FALSE,
                            full_mantissa = FALSE) {
        if (is.null(engine)) engine <- ""
        if (!is.character(engine) || length(engine) != 1L)
          stop("engine must be a single character string")
        if (length(bitexact) != 1L || is.na(bitexact))
          stop("bitexact must be TRUE or FALSE")
        if (length(full_mantissa) != 1L || is.na(full_mantissa))
          stop("full_mantissa must be TRUE or FALSE")
        bitexact <- as.logical(bitexact)
        full_mantissa <- as.logical(full_mantissa)
        self$engine <- engine
        self$ptr <- .Call("randompack_create_R", engine, bitexact, PACKAGE = "randompack")
        if (full_mantissa)
          .Call("randompack_full_mantissa_R", self$ptr, TRUE, PACKAGE = "randompack")
        reg.finalizer(
          self,
          function(e) e$.__enclos_env__$private$finalize(),
          onexit = TRUE
        )
        invisible(self)
      }
    )
  ),
  private = list(
    finalize = function() {
      if (!is.null(self$ptr)) {
        .Call("randompack_free_R", self$ptr, PACKAGE = "randompack")
        self$ptr <- NULL
      }
      invisible(NULL)
    }
  )
)
