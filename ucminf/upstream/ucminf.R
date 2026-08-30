################################################################################
#                                                                              #
# R function for the ucminf package                                            #
#                                                                              #
################################################################################

#------------------------------------------------------------------------------#
# ucminf: General-Purpose Unconstrained Non-Linear Optimization                #
#------------------------------------------------------------------------------#

#' General-Purpose Unconstrained Non-Linear Optimization
#'
#' An algorithm for general-purpose unconstrained non-linear optimization. The
#' algorithm is of quasi-Newton type with BFGS updating of the inverse
#' Hessian and soft line search with a trust region type monitoring of the
#' input to the line search algorithm. The interface of \sQuote{ucminf} is
#' designed for easy interchange with \sQuote{optim}.
#'
#' @param par Initial estimate of minimum for \code{fn}.
#' @param fn Objective function to be minimized.
#' @param gr Gradient of objective function. If \code{NULL} a finite difference
#' approximation is used.
#' @param ... Optional arguments passed to the objective and gradient functions.
#' @param control A list of control parameters. See \sQuote{Details}.
#' @param hessian Integer value: \describe{\item{0}{No hessian approximation is
#' returned.} \item{1}{Returns a numerical approximation of the Hessian using
#' \sQuote{hessian} in the package \sQuote{numDeriv}.} \item{2}{Returns final
#' approximation of the inverse Hessian based on the series of BFGS updates
#' during optimization.} \item{3}{Same at 2, but will also return the Hessian
#' (the inverse of 2).}} If a \code{TRUE} or \code{FALSE} value is given it will
#' switch between option 1 or 0.
#'
#' @details
#' The algorithm is documented in (Nielsen, 2000) (see References below)
#' together with a comparison to the Fortran subroutine \sQuote{MINF} and the
#' Matlab function \sQuote{fminunc}. The implementation of \sQuote{ucminf} in
#' \R uses the original Fortran version of the algorithm.
#'
#' The interface in R is designed so that it is very easy to switch
#' between using \sQuote{ucminf} and \sQuote{\link[stats]{optim}}. The
#' arguments \code{par}, \code{fn}, \code{gr}, and \code{hessian}
#' are all the same (with a few extra options for \code{hessian} in
#' \sQuote{ucminf}). The difference is that there is no \code{method}
#' argument in \sQuote{ucminf} and that some of the components in the
#' \code{control} argument are different due to differences in the algorithms.
#'
#' The algorithm can be given an initial estimate of the Hessian for the
#' optimization and it is possible to get the final approximation of the
#' Hessian based on the series of BFGS updates. This extra functionality
#' may be useful for optimization in a series of related problems.
#'
#' The functions \code{fn} and \code{gr} can return \code{Inf} or \code{NaN}
#' if the functions cannot be evaluated at the supplied value, but the
#' functions must be computable at the initial value. The functions
#' are not allowed to return \code{NA}. Any names given to \code{par} will be
#' copied to the vectors passed to \code{fn} and \code{gr}. No
#' other attributes of \code{par} are copied over.
#'
#' The \code{control} argument is a list that can supply any of the
#' following components:
#' \describe{\item{\code{trace}}{If trace is positive then detailed tracing
#' information is printed for each iteration.}
#' \item{\code{grtol}}{The algorithm stops when
#' \eqn{||F'(x)||_\infty \leq }{||F'(x)||_inf <=} grtol, that is when the
#' largest absolute value of the gradient is less than grtol. Default value is
#' \code{grtol = 1e-6}. }
#' \item{\code{xtol}}{The algorithm stops when \eqn{||x-x_p||_2 \leq
#' \textrm{xtol}\cdot(\textrm{xtol} + ||x||_2)}{||x-x_p||_2 <= xtol*(xtol +
#' ||x||_2)}, where \eqn{x_p} and \eqn{x} are the previous and current estimate
#' of the minimizer. Thus the algorithm stops when the last relative step length
#' is sufficiently small. Default value is \code{xtol = 1e-12}.}
#' \item{\code{stepmax}}{Initial maximal allowed step length (radius of
#' trust-region). The value is updated during the optimization. Default value
#' is \code{stepmax = 1}.}
#' \item{\code{maxeval}}{The maximum number of function evaluations.
#' A function evaluation is counted as one evaluation of the objective function
#' and of the gradient function. Default value is \code{maxeval = 500}.}
#' \item{\code{grad}}{Either \sQuote{forward} or \sQuote{central}. Controls
#' the type of finite difference approximation to be used for the gradient if no
#' gradient function is given in the input argument \sQuote{gr}. Default value
#' is \code{grad = 'forward'}.}
#' \item{\code{gradstep}}{Vector of length 2. The step length in finite
#' difference approximation for the gradient. Step length is
#' \eqn{|x_i|\cdot\textrm{gradstep[1]+gradstep[2]}}{|x_i|*gradstep[1]+
#' gradstep[2]}. Default value is \code{gradstep = c(1e-6, 1e-8)}.}
#' \item{\code{invhessian.lt}}{A vector with an initial approximation to the
#' lower triangle of the inverse Hessian. If not given, the inverse Hessian is
#' initialized as the identity matrix. If \code{H0} is the initial hessian
#' matrix then the lower triangle of the inverse of \code{H0} can be found as
#' \code{invhessian.lt = solve(H0)[lower.tri(H0,diag=TRUE)]}.}}
#'
#' @return \code{\link{ucminf}} returns a list of class \code{'ucminf'}
#' containing the following elements:
#'
#' \item{par}{Computed minimizer.}
#' \item{value}{Objective function value at computed minimizer.}
#' \item{convergence}{Flag for reason of termination:
#' \describe{
#' \item{1}{Stopped by small gradient (grtol).}
#' \item{2}{Stopped by small step (xtol).}
#' \item{3}{Stopped by function evaluation limit (maxeval).}
#' \item{4}{Stopped by zero step from line search}
#' \item{-2}{Computation did not start: length(par) = 0.}
#' \item{-4}{Computation did not start: stepmax is too small.}
#' \item{-5}{Computation did not start: grtol or xtol <= 0.}
#' \item{-6}{Computation did not start: maxeval <= 0.}
#' \item{-7}{Computation did not start: given Hessian not pos. definite.}}}
#' \item{message}{String with reason of termination.}
#' \item{hessian, invhessian}{Estimate of (inv.) Hessian at computed minimizer.
#' The type of estimate is given by the input argument \sQuote{hessian}.}
#' \item{invhessian.lt}{The lower triangle of the final approximation to the
#' inverse Hessian based on the series of BFGS updates during optimization.}
#' \item{info}{Information about the search:
#' \describe{
#' \item{maxgradient}{\eqn{||F'(x)||_\infty}{||F'(x)||_inf}, the largest element
#' in the absolute value of the gradient at the computed minimizer.}
#' \item{laststep}{Length of last step.}
#' \item{stepmax}{Final maximal allowed step length.}
#' \item{neval}{Number of calls to both objective and gradient function.}}}
#'
#' @author \sQuote{UCMINF} algorithm design and Fortran code by Hans Bruun
#' Nielsen.
#'
#' K Hervé Dakpo took over maintenance of the package in May. 2023.
#'
#' Implementation in \R by Stig B. Mortensen, \email{stigbm@gmail.com}.
#'
#' Modifications by Douglas Bates <bates@stat.wisc.edu>, Nov. 2010, to
#' support nested optimization and correct issues with printing on Windows.
#'
#' @seealso \code{\link[stats]{optim}}, \code{\link[stats]{nlminb}},
#' \code{\link[stats]{nlm}}.
#'
#' @references Nielsen, H. B. (2000) \sQuote{UCMINF - An Algorithm For
#' Unconstrained, Nonlinear Optimization}, Report IMM-REP-2000-19, Department of
#' Mathematical Modelling, Technical University of Denmark.
#' \url{http://www.imm.dtu.dk/documents/ftp/tr00/tr19_00.pdf}
#'
#' The original Fortran source code was found at
#' \code{http://www2.imm.dtu.dk/projects/hbn_software/ucminf.f}.
#' (That URL is no longer available but archived at
#' \url{https://web.archive.org/web/20050418082240/http://www.imm.dtu.dk/~hbn/Software/ucminf.f}
#' -- Dr Nielsen passed away in 2015). The code has been slightly modified in
#' this package to be suitable for use with \R.
#'
#' The general structure of the implementation in \R is based on the
#' package \sQuote{FortranCallsR} by Diethelm Wuertz.
#'
#' @keywords optimize nonlinear
#' @export
#'
#' @examples
#' ## Rosenbrock Banana function
#' fR <- function(x) (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2
#' gR <- function(x) c(-400 * x[1] * (x[2] - x[1] * x[1]) - 2 * (1 - x[1]),
#'                     200 * (x[2] - x[1] * x[1]))
#' #  Find minimum and show trace
#' ucminf(par = c(2,.5), fn = fR, gr = gR, control = list(trace = 1))
ucminf <- function(par, fn, gr = NULL, ..., control = list(),
  hessian = 0) {
  con <- list(trace = 0, grtol = 1e-06, xtol = 1e-12, stepmax = 1,
    maxeval = 500, grad = "forward", gradstep = c(1e-06,
      1e-08), invhessian.lt = NULL)
  stopifnot(names(control) %in% names(con))
  con[names(control)] <- control
  stopifnot(length(con$gradstep) == 2, con$grad %in% c("forward",
    "central"))
  fnstr <- quote(fn(.x, ...))
  grstr <- quote(gr(.x, ...))
  rho = new.env(parent = environment())
  n <- length(par)
  eps <- c(con$grtol, con$xtol)
  if (!is.null(gr)) {
    grad <- 0
  } else {
    grad <- ifelse(con$grad == "forward", 1, 2)
  }
  iw <- n * ceiling(max(n + 1, (n + 11)/2)) + 10
  w <- rep(0, iw)
  trace <- con$trace > 0
  icontr = 1 + trace + 2 * !is.null(con$invhessian.lt)
  if (!is.null(con$invhessian.lt))
    w[(4 * n + 1):(4 * n + n * (n + 1)/2)] <- con$invhessian.lt
  par0 <- rep(0, n)
  for (i in 1:n) par0[i] = par[i]
  xname <- as.double(rep(0, n))
  names(xname) <- names(par)
  assign(".f", fnstr, envir = rho)
  assign(".gr", grstr, envir = rho)
  assign(".n", as.integer(n), envir = rho)
  assign(".x", xname, envir = rho)
  assign(".par", as.double(par0), envir = rho)
  assign(".stepmax", as.double(con$stepmax), envir = rho)
  assign(".eps", as.double(eps), envir = rho)
  assign(".maxfun", as.integer(con$maxeval), envir = rho)
  assign(".w", as.double(w), envir = rho)
  assign(".iw", as.integer(iw), envir = rho)
  assign(".icontr", as.integer(icontr), envir = rho)
  assign(".grad", as.integer(grad), envir = rho)
  assign(".grstep", as.double(con$gradstep), envir = rho)
  #
  .Call(C_mfopt, rho)
  #
  W <- get(".w", envir = rho)
  icontr <- get(".icontr", envir = rho)
  ans = list(par = get(".par", envir = rho), value = W[1],
    convergence = icontr, message = switch(as.character(icontr),
      `1` = "Stopped by small gradient (grtol).", `2` = "Stopped by small step (xtol).",
      `3` = "Stopped by function evaluation limit (maxeval)",
      `4` = "Stopped by zero step from line search", `-2` = "Computation did not start: length(par) = 0.",
      `-4` = "Computation did not start: stepmax is too small.",
      `-5` = "Computation did not start: grtol or xtol <= 0.",
      `-6` = "Computation did not start: maxeval <= 0.",
      `-7` = "Computation did not start: given hessian not pos. definite.",
      `-8` = "Computation did not start: iw too small."))
  if (0 < icontr) {
    if (hessian == 1) {
      if (requireNamespace("numDeriv", quietly = TRUE)) {
        p0 <- ans$par
        names(p0) <- names(par)
        ans$hessian <- numDeriv::hessian(fn, p0, method = "Richardson",
          ...)
      } else {
        warning("Skipped hessian estimation - package 'numDeriv' must be installed for hessian option 1")
      }
    }
    if (hessian == 2 | hessian == 3) {
      logicMat <- (matrix(-(1:n^2), n, n, byrow = TRUE) +
        matrix(1:n^2, n, n)) <= 0
      COV <- matrix(0, n, n)
      COV[logicMat] <- W[(4 * n + 1):(4 * n + n * (n +
        1)/2)]
      COV <- t(COV) + COV - diag(diag(COV))
      ans$invhessian <- COV
    }
    if (hessian == 3)
      ans$hessian <- solve(COV)
    ans$invhessian.lt <- W[(4 * n + 1):(4 * n + n * (n +
      1)/2)]
    ans$info = c(maxgradient = W[2], laststep = W[3], stepmax = get(".stepmax",
      envir = rho), neval = get(".maxfun", envir = rho))
  }
  if (trace) {
    cat(paste(ans$message, "\n"))
    if (!is.null(ans$info))
      print(ans$info)
  }
  nm <- names(par)
  if (!is.null(nm)) {
    names(ans$par) <- nm
    if (!is.null(ans$hessian))
      colnames(ans$hessian) <- rownames(ans$hessian) <- nm
  }
  class(ans) <- "ucminf"
  return(ans)
}
