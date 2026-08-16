\name{lossdiw}
\alias{lossdiw}
\title{
Loss function
}
\description{
Quadratic loss function for the method of moments
}
\usage{
lossdiw(x, par, eps = 1e-04, nmax=1000)
}
\arguments{
  \item{x}{
a vector of sample values
}
  \item{par}{
a vector of parameters (\eqn{q} and \eqn{\beta})
}
  \item{eps}{
a tolerance error for the computation of first order moments
}
  \item{nmax}{
a first maximum value for the computation of first order moments
}
}
\value{
the value of the quadratic loss function \eqn{L(x; q, \beta)=(E(X; q, \beta)-m_1)^2+(E(X^2; q, \beta)-m_2)^2} where \eqn{m_1} and \eqn{m_2} are the first and second order sample moments.
}
\seealso{
\code{\link{Ediweibull}}
}
\examples{
n<-100
q<-0.5
beta<-2.5
x<-rdiweibull(n, q, beta)
# loss function computed on the true values
lossdiw(x, c(q, beta))
par<-estdiweibull(x, method="M")
# estimates of the parameters through the method of moments
par
# loss function computed on the estimates derived through
# the method of moments
lossdiw(x, par)
# it should be zero (however, smaller than before...)
}
\keyword{distribution}
