\name{tests}
\alias{maha2_test}
\alias{mardia_test}
\title{Formal Tests of Multivariate Normality}
\description{
  Compute formal tests based on the Mahalanobis distances and
  Mahalanobis angles of multivariate normality (including Mardia's
  kurtosis test and Mardia's skewness test).
}
\usage{
maha2_test(x, type = c("ad.test", "ks.test"), dist = c("chi2", "beta"), ...)
mardia_test(x, type = c("kurtosis", "skewness"), method = c("direct", "chol"))
}
\arguments{
  \item{x}{(n, d)-matrix of data.}
  \item{type}{\code{\link{character}} string indicating the type of test:
    \describe{
      \item{\code{"ad.test"}:}{Anderson-Darling test as computed by the
	underlying \code{\link[ADGofTest]{ad.test}()}.}
      \item{\code{"ks.test"}:}{Kolmogorov-Smirnov test as computed by the
	underlying \code{\link{ks.test}()}.}
      \item{\code{"kurtosis"}:}{Mardia's kurtosis test (based on Mahalanobis
	distances).}
      \item{\code{"skewness"}:}{Mardia's skewness test (based on Mahalanobis
	angles).}
    }
  }
  \item{dist}{distribution to check against.}
  \item{method}{method for computing the Mahalanobis angles.}
  \item{\dots}{additional arguments passed to the underlying
    \code{\link[ADGofTest]{ad.test}()} or \code{\link{ks.test}()}.}
}
\value{
  An \code{htest} object (for \code{maha2_test} the one returned
  by the underlying \code{\link[ADGofTest]{ad.test}()} or
  \code{\link{ks.test}()}).
}
\author{Marius Hofert}
\examples{
set.seed(271)
U <- matrix(runif(3 * 200), ncol = 3)
X <- cbind(qexp(U[,1]), qnorm(U[,2:3]))
maha2_test(X) # at the 'edge' of rejecting
maha2_test(X, type = "ks.test") # at the 'edge', too
mardia_test(X) # clearly rejects at 5%
mardia_test(X, type = "skewness") # clearly rejects at 5%
}
\keyword{htest}
