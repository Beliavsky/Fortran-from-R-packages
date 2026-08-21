#' The Vectorized Multivariate Normal Distribution
#'
#' Vectorized versions of \code{\link[mvtnorm]{rmvnorm}} and \code{\link[mvtnorm]{dmvnorm}}
#' from the \pkg{mvtnorm} package, providing random deviates and density for the multivariate
#' normal distribution with varying vectors of means and varying covariance matrices.
#'
#' \samp{rmvnorm(n, m, s)} is equivalent to \samp{rmultinormal(n, m, as.vector(s))}.
#' \samp{dmvnorm(x, m, s)} is equivalent to \samp{dmultinormal(x, m, as.vector(s))}.
#'
#' If \samp{mean} and/or \samp{sigma} is a matrix, the first random deviate will use
#' the first row of \samp{mean} and/or \samp{sigma}, the second will use the second row, ...
#' recycling being permitted by row.
#'
#' If \samp{mean} is a vector of length \samp{l} or is a matrix with \samp{l} columns,
#' \samp{sigma} should be a vector of length \samp{l x l} or a matrix with \samp{l^2} columns.
#'
#' @note The use of a varying sigma may be very time consuming.
#'
#' @param x vector or matrix of quantiles. If x is a matrix, each row is taken to be a quantile.
#' @param n number of observations. If \samp{length(n) > 1}, the length is taken to be the
#'   number required.
#' @param mean vector or matrix of means. If a matrix, each row is taken to be a mean vector.
#'   Default is a vector of 0 of convenient length.
#' @param sigma covariance vector corresponding to the coercion of the covariance matrix into
#'   a vector (if unique for all \samp{n} or \samp{x}) or array of covariance vectors (if
#'   varying). Default is a diagonal matrix of convenient size.
#' @param method matrix decomposition used to determine the matrix root of sigma. Possible
#'   methods are eigenvalue decomposition (\samp{"eigen"}, default), singular value
#'   decomposition (\samp{"svd"}), and Cholesky decomposition (\samp{"chol"}).
#' @param log logical; if \samp{TRUE}, densities d are given as \samp{log(d)}.
#' @keywords distribution
#' @name multinormal
#' @examples
#' ## mean and sigma as vectors
#' (mean <- c(10, 0))
#' (sigma <- matrix(c(1, 2, 2, 10), ncol = 2))
#' sigma <- as.vector(sigma)
#' (x <- matrix(c(9, 8, 1, -1), ncol = 2))
#' round(rmultinormal(10, mean, sigma))
#' dmultinormal(x, mean, sigma)
#'
#' ## mean as matrix
#' (mean <- matrix(c(10, 0, 0, 10), ncol = 2))
#' round(rmultinormal(10, mean, sigma))
#' @export
rmultinormal <- function(n, mean, sigma , method=c("eigen", "svd", "chol"))
{
  if(length(n) == 0) return(n)
  if(length(n) > 1) n <- length(n)
  if (missing(mean)){
	if(is.vector(sigma))  mean <- rep(0, length = sqrt(length(sigma)))
	else  mean <- rep(0, length = sqrt(ncol(sigma)))}
  if (missing(sigma)){
	if(is.vector(mean))  sigma <- as.vector(diag(length(mean)))
	else  sigma <- as.vector(diag(ncol(mean)))}
  if(is.vector(mean))  mean <- matrix(mean,nrow=1)
  nv <- ncol(mean)
  if(nrow(mean) != n)  mean <- matrix(t(mean), ncol=nv, nrow=n, byrow=TRUE)

  if(is.vector(sigma)) {                       # 'classic' rmvnorm to gain time
    if(length(sigma) != (nv^2)) stop("sigma should be a vector of length:",nv^2)
    return(mean + rmvnorm(n, mean = rep(0, nv), sigma = matrix(sigma,ncol=nv), method=method))
        }

# else Sigma is a matrix

    if ((nv^2) != ncol(sigma)) stop("the length/ncol of sigma should be", nv^2, "\n")
    if (nrow(sigma) != n) sigma <- matrix(t(sigma), ncol = nv^2, nrow = n, byrow = TRUE)
    return(t(sapply(1:n, function(x) rmvnorm(1, mean[x, ], matrix(sigma[x,], ncol = nv), method = method))))
 }

#' @rdname multinormal
#' @export
dmultinormal <- function (x, mean, sigma, log = FALSE)
{
  # 'classic' dmvnorm to gain time
    if (is.vector(x))   x <- matrix(x, nrow = 1)
	nv <- ncol(x)
	if (missing(mean))  mean <- rep(0, length = nv)
    if (missing(sigma)) sigma <- as.vector(diag(nv))

	if(is.vector(mean) & is.vector(sigma)) return(dmvnorm(x=x, mean=mean, sigma = matrix(sigma, ncol=nv), log= log))

	if(is.vector(mean))  mean <- matrix(mean,nrow=1)
	if(is.vector(sigma)) sigma <- matrix(sigma,nrow=1)
    if(ncol(sigma) != (nv^2))  stop("the length/ncol of sigma should be", nv^2, "\n")

	if(nrow(mean) != nrow(x))  mean <- matrix(t(mean), ncol=nv, nrow=nrow(x), byrow=TRUE)
	if(nrow(sigma) != nrow(x)) sigma <- matrix(t(sigma), ncol=nv^2, nrow=nrow(x), byrow=TRUE)

    return(t(sapply(1:nrow(x), function(y) dmvnorm(x=x[y,], mean=mean[y, ], sigma=matrix(sigma[y,], ncol = nv), log = log))))
 }
