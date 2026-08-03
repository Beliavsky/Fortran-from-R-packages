\name{findMu}
\alias{findMu} %% keep in sync with ./findSig.Rd
%%                                    ^^^^^^^^^^ ?? TODO ?? one help page also for  findSig() ???
\title{Find  mu(W) and W, given sigma(W) and CLA result}
\description{
  Find \eqn{\mu(W)} and \eqn{W}, given \eqn{\sigma(W}) and
  \code{\link{CLA}} result.
}
\usage{
findMu(Sig0, result, covar, tol.unir = 1e-06, equal.tol = 1e-06)
}
\arguments{
  \item{Sig0}{numeric vector of \eqn{\sigma(W)} values.}
  \item{result}{a \code{\link{list}} with components \code{MS_weight}
    and \code{weights_set} as resulting from \code{\link{CLA}()}.}
  \item{covar}{the same \eqn{n \times n}{n x n} covariance matrix (of
    asset returns) as the argument of \code{\link{CLA}()}.}
  \item{tol.unir}{numeric tolerance passed to \code{\link{uniroot}}.}
  \item{equal.tol}{numeric tolerance to be used in
    \code{\link{all.equal}(.., tolerance = equal.tol)} in the check to
    see if the \eqn{\mu}{mu} of two neighbouring turning points are
    equal.}
}
\value{
  a \code{\link{list}} with components
  \item{Mu}{numeric vector of same length, say \eqn{M}, as \code{Sig0}.}
  \item{weight}{numeric \eqn{n \times M}{n x M} matrix of weights.}
}
\references{
  Master thesis, p.33
}

\seealso{
 \code{\link{findSig}},
 \code{\link{CLA}}, \code{\link{MS}}.
}
\examples{
data(muS.sp500)
## Full data taking too much time for example
if(getRversion() >= "3.6") .Rk <- RNGversion("3.5.0") # for back compatibility & warning
set.seed(2016)
iS <- sample.int(length(muS.sp500$mu), 17)
if(getRversion() >= "3.6") do.call(RNGkind, as.list(.Rk)) # revert
cov17 <- muS.sp500$covar[iS, iS]
CLsp.17 <- CLA(muS.sp500$mu[iS], covar=cov17, lB=0, uB = 1/2)
CLsp.17 # 16 turning points
summary(tpS <- CLsp.17$MS_weights[,"Sig"])
str(s0 <- seq(0.0186, 0.0477, by = 0.0001))
mu.. <- findMu(s0, result=CLsp.17, covar=cov17)
str(mu..)
stopifnot(dim(mu..$weight) == c(17, length(s0)))
plot(s0, mu..$Mu, xlab=quote(sigma), ylab = quote(mu),
     type = "o", cex = 1/4)
points(CLsp.17$MS_weights, col = "tomato", cex = 1.5)
}
\keyword{optimize}
\keyword{arith}
