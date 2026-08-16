#' Solve a nonlinear equation system with multiple roots from multiple initial
#' estimates
#'
#' This function solves a system of nonlinear equations with \code{nleqlsv} for
#' multiple initial estimates of the roots.
#'
#' Each row of \code{x} is a vector of initial estimates for the argument
#' \code{x} of \code{nleqslv}. The function runs \code{nleqslv} for each row of
#' the matrix \code{x}. The first initial value is treated separately and
#' slightly differently from the other initial estimates. It is used to check
#' if all arguments in \code{...} are valid arguments for \code{nleqslv} and
#' the function to be solved. This is done by running \code{nleqslv} with no
#' condition handling. If an error is then detected an error message is issued
#' and the function stops. For the remaining initial estimates \code{nleqslv}
#' is executed silently. Only solutions for which the \code{nleqslv}
#' termination code \code{tcode} equals \code{1} are regarded as valid
#' solutions. The rounded solutions (after removal of duplicates) are used to
#' order the solutions in increasing order. These rounded solutions are not
#' included in the return value of the function.
#'
#' @param x A matrix with each row containing an initial guess of the roots.
#' @param fn A function of \code{x} returning a vector of function values with
#' the same length as the vector \code{x}.
#' @param digits integer passed to \code{\link{round}} for locating and
#' removing duplicate rounded solutions.
#' @param \dots Further arguments to be passed to \code{\link{nleqslv}},
#' \code{fn} and \code{jac}.
#' @return If no solutions are found \code{NULL} is returned. Otherwise a list
#' containing the following components is returned \describe{
#' \item{\code{x}}{a matrix with each row containing a unique solution
#' (unrounded)} \item{\code{xfnorm}}{a vector of the function criterion
#' associated with each row of the solution matrix \code{x}.}
#' \item{\code{fnorm}}{a vector containing the function criterion for every
#' converged result} \item{\code{idxcvg}}{a vector containing the row indices
#' of the matrix of initial estimates for which function value convergence was
#' achieved} \item{\code{idxxtol}}{a vector containing the row indices of the
#' matrix of initial estimates for which x-value convergence was achieved}
#' \item{\code{idxnocvg}}{a vector containing the row indices of the matrix of
#' initial estimates which lead to an \code{nleqslv} termination code > 2}
#' \item{\code{idxfatal}}{a vector containing the row indices of the matrix of
#' initial estimates for which a fatal error occurred that made \code{nleqslv}
#' stop} \item{\code{xstart}}{a matrix of the initial estimates corresponding
#' to the solution matrix} \item{\code{cvgstart}}{a matrix of all initial
#' estimates for which convergence was achieved} }
#' @keywords nonlinear optimize
#' @examples
#'
#' # Dennis Schnabel example 6.5.1 page 149 (two solutions)
#' set.seed(123)
#' dslnex <- function(x) {
#'     y <- numeric(2)
#'     y[1] <- x[1]^2 + x[2]^2 - 2
#'     y[2] <- exp(x[1]-1) + x[2]^3 - 2
#'     y
#' }
#' xstart <- matrix(runif(50, min=-2, max=2),ncol=2)
#' ans <- searchZeros(xstart,dslnex, method="Broyden",global="dbldog")
#' ans
#'
#' # more complicated example
#' # R. Baker Kearfott, Some tests of Generalized Bisection,
#' # ACM Transactions on Methematical Software, Vol. 13, No. 3, 1987, pp 197-220
#'
#' # A high-degree polynomial system (section 4.3 Problem 12)
#' # There are 12 real roots (and 126 complex roots to this system!)
#'
#' hdp <- function(x) {
#'     f <- numeric(length(x))
#'     f[1] <- 5 * x[1]^9 - 6 * x[1]^5 * x[2]^2 + x[1] * x[2]^4 + 2 * x[1] * x[3]
#'     f[2] <- -2 * x[1]^6 * x[2] + 2 * x[1]^2 * x[2]^3 + 2 * x[2] * x[3]
#'     f[3] <- x[1]^2 + x[2]^2 - 0.265625
#'     f
#' }
#'
#'
#' N <- 40 # at least to find all 12 roots
#' set.seed(123)
#' xstart <- matrix(runif(3*N,min=-1,max=1), N, 3)  # N initial guesses, each of length 3
#' ans <- searchZeros(xstart,hdp, method="Broyden",global="dbldog")
#' ans$x
#'
#' @export searchZeros
searchZeros <- function(x, fn, digits=4L, ... )
{
    if( !is.numeric(x) ) stop('argument x should be numeric')
    if( !is.matrix(x)  ) stop('argument x must be a matrix')
    if( any(is.na(x))  ) stop("argument x may not contain NA")
    if( !is.numeric(digits) ) stop('argument digits should be numeric')
    if( is.na(digits)  ) digits <- 4L

    N <- nrow(x)
    if( N < 1 ) stop("Matrix 'x' must have at least 1 row")

    tcode <- numeric(N)
    fnorm <- numeric(N)
    xmat  <- matrix(NA, nrow=N, ncol=ncol(x))

    kerr <- numeric(N)
    kptr <- 0

    # for k==1 check that all arguments are correct --> no try

    for ( k in seq_len(N) ){
        if( k == 1 ) {
            z <- nleqslv(x[k, ], fn=fn, ...)
        } else {
            z <- try(nleqslv(x[k, ], fn=fn, ...), silent=TRUE)
            if( inherits(z, "try-error") ) {
                kptr <- kptr+1
                kerr[kptr] <- k
                next
            }
        }
        tcode[k]  <- z$termcd
        fnorm[k]  <- norm(z$fvec,"2")^2/2 # criterion for global methods
        xmat[k, ] <- z$x
    }

    # locate index of converged trials and store corresponding starting values
    # return NULL if no full convergence obtained
    if(!any(tcode==1)) return(NULL)
    idxcvg <- which(tcode==1)
    xstartcvg <- x[idxcvg,,drop=FALSE]
    # rounded solutions for locating duplicates and remove duplicates
    xsol <- round(xmat[idxcvg,,drop=FALSE], digits)
    notdups <- !duplicated(xsol)
    xsol <- xsol[notdups,,drop=FALSE]
    solstart <- xstartcvg[notdups,,drop=FALSE]
    if( !is.null(colnames(x)) ) {
        colnames(xmat) <- colnames(x)
        colnames(xsol) <- colnames(x)
        colnames(solstart) <- colnames(x)
    }

    # order the rounded solution
    if( nrow(xsol) > 1 ) {
        zidxo <- do.call(order,split(xsol,col(xsol)))
    } else {
        zidxo <- 1
    }

    idxfatal <- if(kptr) kerr[1:kptr] else integer(0)
    idxxtol   <- which(tcode==2)
    idxnocvg <- which(tcode>2)
    # original unrounded solutions with duplicates (above) removed
    xsol <- xmat[idxcvg,,drop=FALSE][notdups,,drop=FALSE]

    # return full precision solutions ordered with rounded ordering
    res <- list(x=xsol[zidxo,,drop=FALSE], xfnorm=fnorm[idxcvg][notdups][zidxo],
                fnorm=fnorm[idxcvg], idxcvg=idxcvg, idxxtol=idxxtol,
                idxnocvg=idxnocvg, idxfatal=idxfatal,
                xstart=solstart[zidxo,,drop=FALSE],cvgstart=xstartcvg)
    res
}
