#' @importFrom R6 R6Class
#' @useDynLib randompack, .registration = TRUE
#'
#' @section Quick Start:
#' Create an RNG instance with \code{randompack_rng()} and use its methods to
#' generate random variates:
#'
#' \preformatted{
#'   rng <- randompack_rng()                # Default: x256++simd
#'   rng <- randompack_rng("pcg64")         # Specify engine
#'   x <- rng$normal(100, mu=1, sigma=2)    # Generate 100 N(1,2) variates
#' }
#'
#' @section Available Engines:
#' Among the supported underlying random number generators (engines) are
#' xoshiro256++, xoshiro256**, PCG64 DXSM, sfc64, Philox-4x64, and ChaCha20.
#' See \code{\link{randompack_rng}} for a complete list.
#'
#' @seealso
#' \code{\link{randompack_rng}} for creating and using random number generators;
#' \code{\link{randompack_engines}} for a list of available engines
"_PACKAGE"
## usethis namespace: start
## usethis namespace: end
NULL
