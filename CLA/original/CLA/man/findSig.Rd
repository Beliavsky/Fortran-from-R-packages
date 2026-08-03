\name{findSig}
\alias{findSig} %% keep in sync with ./findMu.Rd
\title{Find sigma(W) and W, given mu(W) and CLA result}
\description{
  Find \eqn{\sigma(W)} and \eqn{W}, given \eqn{\mu(W}) and
  \code{\link{CLA}} result.
}
\usage{
findSig(Mu0, result, covar, equal.tol)
}
\arguments{
  \item{Mu0}{numeric vector of \eqn{\mu(W)} values.}
  \item{result}{a \code{\link{list}} with components \code{MS_weight}
    and \code{weights_set} as resulting from \code{\link{CLA}()}.}
  \item{covar}{the same \eqn{n \times n}{n x n} covariance matrix (of
    asset returns) as the argument of \code{\link{CLA}()}.}
  \item{equal.tol}{numeric tolerance to be used in
    \code{\link{all.equal}(.., tolerance = equal.tol)} in the check to
    see if the \eqn{\mu}{mu} of two neighbouring turning points are
    equal.}
}
%% \details{
%% %%  ~~ If necessary, more details than the description above ~~
%% }
\value{
  a \code{\link{list}} with components
  \item{Sig}{numeric vector of same length, say \eqn{M}, as \code{Mu0}.}
  \item{weight}{numeric \eqn{n \times M}{n x M} matrix of weights.}
}
\references{
  Master thesis, p.33
}
\seealso{
 \code{\link{findMu}},
 \code{\link{CLA}}, \code{\link{MS}}.
}
\examples{
data(muS.sp500)
## Full data taking too much time for example: Subset of n=21:
if(getRversion() >= "3.6") .Rk <- RNGversion("3.5.0") # for back compatibility & warning
set.seed(2018)
iS <- sample.int(length(muS.sp500$mu), 21)
if(getRversion() >= "3.6") do.call(RNGkind, as.list(.Rk)) # revert
cov21 <- muS.sp500$covar[iS, iS]
CLsp.21 <- CLA(muS.sp500$mu[iS], covar=cov21, lB=0, uB = 1/2)
CLsp.21 # 14 turning points
summary(tpM <- CLsp.21$MS_weights[,"Mu"])
str(m0 <- c(min(tpM),seq(0.00205, 0.00525, by = 0.00005), max(tpM)))
sig. <- findSig(m0, result=CLsp.21, covar=cov21)
str(sig.)
stopifnot(dim(sig.$weight) == c(21, length(m0)))
plot(sig.$Sig, m0, xlab=quote(sigma), ylab = quote(mu),
     type = "o", cex = 1/4)
points(CLsp.21$MS_weights, col = "tomato", cex = 1.5)
title("Efficient Frontier from CLA()")
mtext("findSig() to interpolate between turning points", side=3)
}
\keyword{arith}
\keyword{dplot}
