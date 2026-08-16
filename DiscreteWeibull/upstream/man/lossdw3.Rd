\name{lossdw3}
\alias{lossdw3}
\title{
Loss function
}
\description{
Loss function for the method of moments (type 3 discrete Weibull)
}
\usage{
lossdw3(par, x, eps = 1e-04)
}
\arguments{
  \item{par}{
vector of parameters \eqn{q} and \eqn{\beta}
}
  \item{x}{
the vector of sample values
}
  \item{eps}{
error threshold for the numerical computation of the expected value
}
}
\details{
The loss function is given by \eqn{L(x;c,\beta)=[m_1-\mathrm{E}(X;c,\beta)]^2+[m_2-\mathrm{E}(X^2;c,\beta)]^2}, where \eqn{\mathrm{E}(\cdot)} denotes the expected value, \eqn{m_1} and \eqn{m_2} are the first and second order sample moments respectively.
}
\value{
the value of the quadratic loss function
}
\author{
Alessandro Barbiero
}
\seealso{
\code{\link{Edweibull3}}
}
\examples{
n <- 25
c <- 1/3
beta <- 2/3
x <- rdweibull3(n, c, beta)
par <- estdweibull3(x, "M")
par
lossdw3(par, x)
}
\keyword{distribution}

