\name{Ediweibull}
\alias{Ediweibull}
\title{
First and second order moments
}
\description{
First and second order moments of the discrete inverse Weibull distribution
}
\usage{
Ediweibull(q, beta, eps = 1e-04, nmax = 1000)
}
\arguments{
  \item{q}{
the value of the \eqn{q} parameter
}
  \item{beta}{
the value of the \eqn{\beta} parameter
}
  \item{eps}{
error threshold for the approximated computation of the moments
}
  \item{nmax}{
a first maximum value of the support considered for the approximated computation of the moments
}
}
\details{
For a discrete inverse Weibull distribution we have \eqn{E(X;q,\beta)=\sum_{x=0}^{+\infty} 1-F(x;q, \beta)} and \eqn{E(X^2;q,\beta)=2\sum_{x=1}^{+\infty} x(1-F(x;q, \beta))+E(X;q, \beta)}.
The expected values are numerically computed considering a truncated support: integer values smaller than or equal to \eqn{\min(nmax;F^{-1}(1-eps;q,\beta))}, where \eqn{F^{-1}} is the inverse of the cumulative distribution function (implemented by the function \code{\link{qdiweibull}}). Increasing the value of \code{nmax} or decreasing the value of \code{eps} improves the approximation, but slows down the calculation speed
}
\value{
a list comprising the (approximate) first and second order moments of the discrete inverse Weibull distribution. Note that the first moment is finite iff \eqn{\beta} is greater than 1; the second order moment is finite iff \eqn{\beta} is greater than 2
}
\references{
Khan M.S., Pasha G.R., Pasha A.H. (2008) Theoretical Analysis of Inverse Weibull Distribution, WSEAS Trabsactions on Mathematics 2(7): 30-38
}
\examples{
# Ex.1
q<-0.75
beta<-1.25
Ediweibull(q, beta)
# Ex.2 
q<-0.5
beta<-2.5
Ediweibull(q, beta)
# Ex.3
q<-0.4
beta<-4
Ediweibull(q, beta)
}
\keyword{distribution}
