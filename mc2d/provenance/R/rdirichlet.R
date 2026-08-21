#' The Dirichlet Distribution
#'
#' Density function and random generation from the Dirichlet distribution.
#'
#' The Dirichlet distribution is the multidimensional generalization of the beta distribution.
#' The original code was adapted to provide a kind of "vectorization" used in multivariate
#' \samp{mcnode} objects.
#'
#' @param x a vector containing a single deviate or a matrix containing one random deviate per row.
#' @param alpha a vector of shape parameters, or a matrix of shape parameters by rows.
#'   Recycling (by row) is permitted.
#' @param n number of random vectors to generate. If \samp{length(n) > 1}, the length is
#'   taken to be the number required.
#' @return \samp{ddirichlet} gives the density. \samp{rdirichlet} returns a matrix with
#' \samp{n} rows, each containing a single Dirichlet random deviate.
#' @note Code is adapted from \samp{MCMCpack}. It originates from Greg's Miscellaneous
#' Functions (gregmisc).
#' @seealso \code{\link{Beta}}
#' @keywords distribution
#' @name dirichlet
#' @examples
#' dat <- c(1, 10, 100, 1000, 1000, 100, 10, 1)
#' (alpha <- matrix(dat, nrow = 4, byrow = TRUE))
#' round(x <- rdirichlet(4, alpha), 2)
#' ddirichlet(x, alpha)
#'
#' ## rdirichlet used with mcstoc
#' mcalpha <- mcdata(dat, type = "V", nsv = 4, nvariates = 2)
#' (x <- mcstoc(rdirichlet, type = "V", alpha = mcalpha, nsv = 4, nvariates = 2))
#' unclass(x)
#' @export
ddirichlet <- function (x, alpha)
{
    if(length(x) == 0) return(numeric(0))

	if (!is.matrix(x)) x <- t(x)
    if (!is.matrix(alpha)) alpha <- t(alpha)

    ncx <- ncol(x)
    nca <- ncol(alpha)
    if(ncx != nca) stop("x and alpha should be of same length or have the same number of columns.")

    nrx <- nrow(x)
    nra <- nrow(alpha)

    if(nra != nrx) alpha <- matrix(t(alpha), ncol = nca, nrow = nrx, byrow = TRUE)

	lgama <- lgamma(alpha)
	logD <- apply(lgama,1,sum) - lgamma(apply(alpha,1,sum))
	s <- apply((alpha - 1) * log(x), 1, sum)
	pd <- exp(s - logD)

	pd[apply(x, 1, function(z) any(z < 0 | z > 1))] <- 0
    pd[apply(x, 1, function(z) all.equal(sum(z), 1) != TRUE)] <- 0
    return(pd)
}

#' @rdname dirichlet
#' @export
rdirichlet <- function (n, alpha)
{
  if(length(n) > 1) n <- length(n)
  if(length(n) == 0 || as.integer(n) == 0) return(numeric(0))
  n <- as.integer(n)
  if(n < 0) stop("integer(n) can not be negative in rdirichlet")

  if(is.vector(alpha)) alpha <- t(alpha)
  l <- dim(alpha)[2]
  x <- matrix(rgamma(l * n, t(alpha)), ncol = l, byrow=TRUE)  # Gere le recycling
  return(x / rowSums(x))
}
