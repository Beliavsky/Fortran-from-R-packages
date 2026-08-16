\name{estdlaplace2}
\alias{estdlaplace2}
\title{
Sample estimation for the ADSL
}
\description{
The function provides the point estimates for the parameters of the ASDL, resorting to four possible methods: method of moments, maximum likelihood method, method of proportion, modified method of moments. For details, please take a look at the references.
}
\usage{
estdlaplace2(x, method = "M", err = 0.001, parml = c(exp(-1), exp(-1)))
}
\arguments{
  \item{x}{
a vector of observations from the ADSL
}
  \item{method}{
\code{M} for the method of moments, \code{ML} for the maximum likelihood methods, \code{P} for the method of proportion, \code{MM} for the modified method of moments
}
  \item{err}{
a positive tolerance value, as small as possible, used in the definition of lower and upper bounds of the parameters \eqn{p} and \eqn{q} in the minimization algorithm utilized by the method of moments
}
  \item{parml}{
starting values for \eqn{p} and \eqn{q} in the optimization process for the maximum likelihood method
}
}
\value{
a vector with the parameter estimates of \eqn{p} and \eqn{q}.
}
\references{
A. Barbiero, An alternative discrete Laplace distribution, \emph{Statistical Methodology}, 16: 47-67
}
\author{
Alessandro Barbiero, Riccardo Inchingolo
}
\seealso{
\code{\link{dlaplacelike2}}
}
\examples{
p <- 0.4
q <- 0.6
x <- rdlaplace2(n=100, p, q)
est <- matrix(0, 5, 2)
est[1,] <- c(p,q)
est[2,] <- estdlaplace2(x, method="M")
est[3,] <- estdlaplace2(x, method="ML")
est[4,] <- estdlaplace2(x, method="P")
est[5,] <- estdlaplace2(x, method="MM")
dimnames(est)[[1]]<-c("true","M","ML","P","MM")
dimnames(est)[[2]]<-c("p","q")
xlim <- c(min(est[,1])*.98,max(est[,1])*1.02)
ylim <- c(min(est[,2])*.98,max(est[,2])*1.02)
plot(est, pch=19, col=1:5, xlim=xlim, ylim=ylim)
text(est, dimnames(est)[[1]], pos=3, col=1:5, cex= .75)
}
\keyword{distribution,htest}

