\name{iFI2}
\alias{iFI2}
\title{
Inverse of Fisher Information matrix
}
\description{
Inverse of Fisher Information matrix for the ADSL.
}
\usage{
iFI2(p, q)
}
\arguments{
  \item{p}{
first parameter \eqn{p}, in \eqn{(0,1)}, of the ADSL
}
  \item{q}{
second parameter \eqn{q}, in \eqn{(0,1)}, of the ADSL
}
}
\value{
The inverse of Fisher Information matrix. Take a look at the references for more details.
}
\references{
A. Barbiero, An alternative discrete Laplace distribution, \emph{Statistical Methodology}, 16: 47-67
}
\author{
Alessandro Barbiero, Riccardo Inchingolo
}
\seealso{
\code{\link{ioFI2}}
}
\examples{
p <- 0.2
q <- 0.8
iFI2(p, q)
}
\keyword{distribution,htest}
