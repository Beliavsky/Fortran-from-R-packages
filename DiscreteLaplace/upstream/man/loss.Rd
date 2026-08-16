\name{loss}
\alias{loss}
\title{
Loss function for the method of moments
}
\description{
A loss function used for the implementation of the method of moments  (for the ADSL).
}
\usage{
loss(par, x)
}
\arguments{
  \item{par}{
the vector of parameters, \code{p} and \code{q}, of the ADSL
}
  \item{x}{
a vector of sample values from the ADSL
}
}
\value{
The value \eqn{L=[E(X)-m_1(x)]^2+[E(X^2)-m_2(x)]^2}, where \eqn{m_1} and \eqn{m_2} are the first- and second-order sample moments.
}
\references{
A. Barbiero, An alternative discrete Laplace distribution, \emph{Statistical Methodology}, 16: 47-67
}
\author{
Alessandro Barbiero, Riccardo Inchingolo
}
\seealso{
\code{\link{estdlaplace2}}, \code{\link{Edlaplace2}}
}
\examples{
p <- 0.3
q <- 0.7
x <- rdlaplace2(n=100, p, q)
par <- estdlaplace2(x, "M")
loss(par, x) # should be near zero
loss(c(p,q), x) # may be far greater than zero
}
\keyword{distribution}