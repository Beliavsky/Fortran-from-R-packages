\name{ioFI2}
\alias{ioFI2}
\title{
Inverse of the observed Fisher Information matrix
}
\description{
Inverse of the observed Fisher Information matrix computed on a random sample of ADSL values.
}
\usage{
ioFI2(x)
}
\arguments{
  \item{x}{
a vector of observations from the ADSL
}
}
\value{
The inverse of the observed Fisher Information matrix.
}
\references{
A. Barbiero, An alternative discrete Laplace distribution, \emph{Statistical Methodology}, 16: 47-67
}
\author{
Alessandro Barbiero, Riccardo Inchingolo
}
\seealso{
\code{\link{iFI2}}
}
\examples{
n <- 100
p <- 0.4
q <- 0.7
x <- rdlaplace2(n, p, q)
M <- ioFI2(x)
par <- estdlaplace2(x, "ML")
se <- diag(sqrt(M))
par # MLEs
se # standard errors
M # compare with  the inverse of Fisher Information matrix
iFI2(par[1], par[2])/n # with MLEs plugged in
iFI2(p, q)/n # or the true values
}
\keyword{htest}
\keyword{distribution}
