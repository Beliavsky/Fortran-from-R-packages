\name{muSigmaGarch}
\alias{muSigmaGarch}
\title{Compute (mu, Sigma) for a Set of Assets via GARCH fit}
\description{
  Compute (mu, Sigma) for a set of assets via a GARCH fit to each
  individual asset, using package
  \CRANpkg{fGarch}'s \code{\link[fGarch]{garchFit}()}.
}
\usage{
muSigmaGarch(x, formula = ~garch(1, 1), cond.dist = "std", trace = FALSE,
             \dots)
}
\arguments{
  \item{x}{numeric matrix or data frame (\eqn{T \times d}{T * n}) of log
    returns of \eqn{d} assets, observed on a common set of \eqn{T} time points.
  }
  \item{formula}{optional formula for \code{\link[fGarch]{garchFit}}.}
  \item{cond.dist}{the conditional distribution to be used for the garch process.}
  \item{trace}{logical indicating if some progress of \code{garchFit()}
    should printed to the console.}
  \item{\dots}{optional arguments to \code{\link{cor}}, i.e., \code{use}
    or \code{method}.}
}
%% \details{
%% %%  ~~ If necessary, more details than the description above ~~
%% }
\value{
 a list with components
 \item{mu}{numeric vector of length \eqn{n} of mean returns (\eqn{= E[R_i]}).}
 \item{covar}{covariance matrix (\eqn{n \times n}{n * n}) of the returns.}
}
%% \references{
%% %% ~put references to the literature/web site here ~
%% }
\seealso{
  \code{\link{muS.sp500}} which has been produced via \code{muSigmaGarch}.
  \code{\link{CLA}} which needs \code{(mu, covar)} as crucial input.
}
\examples{
if(requireNamespace("FRAPO")) {
  data(NASDAQ, package = "FRAPO")
  ## 12 randomly picked stocks from NASDAQ data
  iS <- if(FALSE) {  ## created (w/ warning, in new R)  by
    RNGversion("3.5.0"); set.seed(17); iS <- sample(ncol(NASDAQ), 12)
  } else c(341L, 2126L, 1028L, 1704L, 895L, 1181L, 454L, 410L, 1707L, 425L, 950L, 5L)
  X. <- NASDAQ[, iS]
  muSig <- muSigmaGarch(X.)
  stopifnot(identical(names(muSig$mu), names(NASDAQ)[iS]),
            identical(dim(muSig$covar), c(12L,12L)),
     all.equal(unname(muSig$mu),
               c(  7.97, -4.05, -14,     21.5, -5.36, -15.3,
                 -15.9,  11.8,   -1.64, -14,    3.13, 121) / 10000,
               tol = 0.0015))
}% if( <FRAPO> )
}
\keyword{multivariate}
