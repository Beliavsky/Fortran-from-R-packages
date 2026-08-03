\name{muS.sp500}
\alias{muS.sp500}
\docType{data}
\title{Return Expectation and Covariance for "FRAPO"s SP500 data}
\description{
  If \eqn{R_{j,t}} are the basically the scale standardized log returns
  for \eqn{j = 1,2,\dots,476} of 476 stocks from \acronym{S&P 500}, as
  from \code{\link[FRAPO]{SP500}}, then \eqn{mu_j = E[R_{j,*}]} somehow
  averaged over time; actually as predicted by \code{muSigma()} at the
  end of the time period, and
  \eqn{\Sigma_{j,k} = Cov(R_j, R_k)}{Sigma[j,k] = Cov(R[j], R[k])}
  are estimated covariances.

  These are the main \dQuote{inputs} needed for the CLA algorithm, see
  \code{\link{CLA}}.
}
\usage{
data("muS.sp500")
}
\format{
  A list with two components,
  \describe{
    \item{mu}{Named num [1:476] 0.00233 0.0035 0.01209 0.00322 0.00249 ...\cr
            names : chr [1:476] "A" "AA" "AAPL" "ABC" ...}
    \item{covar}{num [1:476, 1:476] 0.001498 0.000531 0.000536 ...}
  }
}
\source{
  It is as simple as this:\preformatted{
    data(SP500, package="FRAPO")
    system.time(muS.sp500 <- muSigmaGarch(SP500))   #   26 sec. (lynne, 2017)
  }
}
\seealso{
  \code{\link{muSigmaGarch}()} which was used to construct it.
}
%% \references{
%% %%  ~~ possibly secondary sources and usages ~~
%% }
\examples{
data(muS.sp500)
str(muS.sp500)
}
\keyword{datasets}
