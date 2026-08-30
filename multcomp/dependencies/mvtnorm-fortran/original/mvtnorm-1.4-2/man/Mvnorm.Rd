\name{Mvnorm}
\alias{dmvnorm}
\alias{rmvnorm}
\title{Multivariate Normal Density and Random Deviates}
\description{
  These functions provide the density function and a random number
  generator for the multivariate normal
  distribution with mean equal to \code{mean} and covariance matrix
  \code{sigma}.
}
\usage{
dmvnorm(x, mean = rep(0, p), sigma = diag(p), log = FALSE, checkSymmetry = TRUE)
rmvnorm(n, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)),
           method=c("eigen", "svd", "chol"), pre0.9_9994 = FALSE, 
           checkSymmetry = TRUE, rnorm = stats::rnorm)
}
\arguments{
 \item{x}{vector or matrix of quantiles. When \code{x} is a matrix, each
          row is taken to be a quantile and columns correspond to the number of dimensions, 
          \code{p}.}
 \item{n}{number of observations.}
 \item{mean}{mean vector, default is \code{rep(0, length = ncol(x))}.
             In \code{ldmvnorm} or 
             \code{sldmvnorm}, \code{mean} is a matrix with
             observation-specific means arranged in columns.}
 \item{sigma}{covariance matrix, default is \code{diag(ncol(x))}.}
 \item{log}{logical; if \code{TRUE}, densities d are given as log(d).}
 \item{method}{string specifying the matrix decomposition used to
   determine the matrix root of \code{sigma}.  Possible methods are
   eigenvalue decomposition (\code{"eigen"}, default),
   singular value decomposition (\code{"svd"}), and
   Cholesky decomposition (\code{"chol"}).  The
   Cholesky is typically fastest, not by much though.}
 \item{pre0.9_9994}{logical; if \code{FALSE}, the output produced in mvtnorm
                    versions up to 0.9-9993 is reproduced. In 0.9-9994, the
                    output is organized such that \code{rmvnorm(10,...)} has the
                    same first ten rows as \code{rmvnorm(100, ...)} when called
                    with the same seed.}
 \item{checkSymmetry}{logical; if \code{FALSE}, skip checking whether the
   covariance matrix is symmetric or not. This will speed up the
   computation but may cause unexpected outputs when ill-behaved
   \code{sigma} is provided. The default value is \code{TRUE}.}
 \item{rnorm}{a function with the same interface as
   \code{\link[stats]{rnorm}}. This allows switching to other generators of
   standard normal variables.}
}
\details{
\code{dmvnorm} computes the density function of the multivariate normal
specified by mean and the covariance matrix \code{sigma}. 

\code{rmvnorm} generates multivariate normal variables. 

}
\seealso{\code{\link{pmvnorm}}, \code{\link{rnorm}}, \code{\link{qmvnorm}},
         \code{vignette("lmvnorm_src", package = "mvtnorm")}}
\examples{
dmvnorm(x=c(0,0))
dmvnorm(x=c(0,0), mean=c(1,1))

sigma <- matrix(c(4,2,2,3), ncol=2)
x <- rmvnorm(n=500, mean=c(1,2), sigma=sigma)
colMeans(x)
var(x)
dS <- dmvnorm(x, sigma = sigma)

### alternative interface
C <- t(chol(sigma))
(C <- ltMatrices(C[lower.tri(C, diag = TRUE)], diag = TRUE))
dC <- exp(ldmvnorm(obs = t(x), chol = C, logLik = FALSE))
all.equal(dS, dC)

x <- rmvnorm(n=500, mean=c(1,2), sigma=sigma, method="chol")
colMeans(x)
var(x)

plot(x)
}
\keyword{distribution}
\keyword{multivariate}
