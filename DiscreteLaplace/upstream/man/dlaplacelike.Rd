\name{dlaplacelike2}
\alias{dlaplacelike2}
\title{
Log-likelihood function for the ADSL distribution
}
\description{
Log-likelihood function (changed in sign) for the ADSL distribution.
}
\usage{
dlaplacelike2(par, x)
}
\arguments{
  \item{par}{
the vector of parameters \eqn{(p,q)}
}
  \item{x}{
a vector of observations from ADSL
}
}
\value{
The log-likelihood function with changed sign.
}
\references{
A. Barbiero, An alternative discrete Laplace distribution, \emph{Statistical Methodology}, 16: 47-67
}
\author{
Alessandro Barbiero, Riccardo Inchingolo
}
\seealso{
\code{\link{estdlaplace2}}, \code{\link{ddlaplace2}}
}
\examples{
p <- 0.25
q <- 0.7
x <- rdlaplace2(n=100, p, q)
par <- estdlaplace2(x, "ML")
-dlaplacelike2(par, x) # greater than...
-dlaplacelike2(c(p, q), x)
}
\keyword{htest}
