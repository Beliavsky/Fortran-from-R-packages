\name{lossdw}
\alias{lossdw}
\title{
Loss function
}
\description{
Loss function for the method of moments (type 1 discrete Weibull)
}
\usage{
lossdw(par, x, zero = FALSE, eps = 1e-04, nmax=1000)
}
\arguments{
  \item{par}{
vector of parameters \eqn{q} and \eqn{\beta}
}
  \item{x}{
the vector of sample values
}
  \item{zero}{
\code{TRUE}, if the support contains \eqn{0}; \code{FALSE} otherwise
}
  \item{eps}{
error threshold for the numerical computation of the expected value
}
  \item{nmax}{
maximum value considered for the numerical computation of the expected value
}
}
\details{
The loss function is given by \eqn{L(x;q,\beta)=[m_1-\mathrm{E}(X;q,\beta)]^2+[m_2-\mathrm{E}(X^2;q,\beta)]^2}, where \eqn{\mathrm{E}(\cdot)} denotes the expected value, \eqn{m_1} and \eqn{m_2} are the first and second order sample moments respectively.
}
\value{
the value of the quadratic loss function
}
\author{
Alessandro Barbiero
}
\seealso{
\code{\link{Edweibull}}
}
\examples{
x <- c(1,1,1,1,1,2,2,2,3,4)
lossdw(c(0.5, 1), x)
par <- estdweibull(x, "M") # parameter estimates derived by the method of moments
par
lossdw(par, x) # the loss is zero using these estimates
}
\keyword{distribution}

