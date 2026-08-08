#' @title Non-Negative CG Minimizer
#' 
#' @description Minimize a differentiable function subject to all the variables being non-negative
#' (i.e. >= 0), using a Conjugate-Gradient algorithm based on a modified Polak-Ribiere-Polyak formula (see
#' reference at the bottom for details).
#' 
#' @param evaluate_function function(x, ...) objective evaluation function
#' @param evaluate_gradient function(x, ...) gradient evaluation function
#' @param x0 Starting point. Must be a feasible point (>=0). Be aware that it might be modified in-place.
#' @param tol Tolerance for <gradient, direction>
#' @param maxnfeval Maximum number of function evaluations
#' @param maxiter Maximum number of CG iterations
#' @param decr_lnsrch Number by which to decrease the step size after each unsuccessful line search
#' @param lnsrch_const Acceptance parameter for the line search procedure
#' @param max_ls Maximum number of line search trials per iteration
#' @param extra_nonneg_tol Ensure extra non-negative tolerance by explicitly setting elements that
#' are <=0 to zero at each iteration
#' @param nthreads Number of parallel threads to use (ignored if the package was installed from CRAN)
#' @param verbose Whether to print convergence messages
#' @param ... Extra parameters to pass to the objective and gradient functions
#' @export
#' @examples
#' fr <- function(x) {   ## Rosenbrock Banana function
#'   x1 <- x[1]
#'   x2 <- x[2]
#'   100 * (x2 - x1 * x1)^2 + (1 - x1)^2
#' }
#' grr <- function(x) { ## Gradient of 'fr'
#'   x1 <- x[1]
#'   x2 <- x[2]
#'   c(-400 * x1 * (x2 - x1 * x1) - 2 * (1 - x1),
#'     200 *      (x2 - x1 * x1))
#' }
#' minimize.nonneg.cg(fr, grr, x0 = c(0,2), verbose=TRUE, tol=1e-8)
#' @details The underlying C function can also be called directly from Rcpp with `R_GetCCallable` (see example of such usage
#' in the source code of the 'zoo' package).
#' @references Li, C. (2013). A conjugate gradient type method for the nonnegative constraints optimization problems. Journal of Applied Mathematics, 2013.
minimize.nonneg.cg <- function(evaluate_function, evaluate_gradient, x0,
                               tol=1e-4, maxnfeval=1500, maxiter=200,
                               decr_lnsrch=.5, lnsrch_const=.01, max_ls=20,
                               extra_nonneg_tol=FALSE, nthreads=1, verbose=FALSE, ...){

  if (!(decr_lnsrch >0 & decr_lnsrch <1 )){
    stop("'decr_lnsrch' must be between zero and one.")
  }
  if(! (lnsrch_const > 0 & lnsrch_const < 1) ){
    stop("'lnsrch_const' must be between zero and one.")
  }
  if (tol < 0){
    stop("'tol' must be greater than zero.")
  }
  if (nthreads < 1){nthreads = 1}
  if (min(x0) < 0){stop("'x0' must be a feasible point.")}
  x0 = as.numeric(x0)


  tol = as.numeric(tol)
  decr_lnsrch = as.numeric(decr_lnsrch)
  lnsrch_const = as.numeric(lnsrch_const)
  maxnfeval = as.integer(maxnfeval)
  maxiter = as.integer(maxiter)
  max_ls = as.integer(max_ls)
  extra_nonneg_tol = as.integer(extra_nonneg_tol)
  nthreads = as.integer(nthreads)
  verbose = as.logical(verbose)

  return(
    run_minimizer(evaluate_function, evaluate_gradient, do.call, x0,
                  list(...), tol, maxnfeval, maxiter, decr_lnsrch, lnsrch_const,
                  max_ls, extra_nonneg_tol, nthreads, verbose))
}
