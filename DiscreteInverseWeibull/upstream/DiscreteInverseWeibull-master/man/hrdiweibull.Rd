\name{hrdiweibull}
\alias{hrdiweibull}
\title{
Hazard rate function
}
\description{
Hazard rate function for the discrete inverse Weibull distribution
}
\usage{
hrdiweibull(x, q, beta)
}
\arguments{
  \item{x}{
a vector of values
}
  \item{q}{
the value of the \eqn{q} parameter
}
  \item{beta}{
the value of the \eqn{\beta} parameter
}
}
\value{
the hazard rate function computed on the \code{x} values
}
\details{
The hazard rate function is defined as \eqn{r(x)=P(X=x)/P(X\ge x)=(q^{x^{-\beta}}-q^{(x-1)^{-\beta}})/(1-q^{(x-1)^{-\beta}})}
}
\seealso{
\code{\link{ahrdiweibull}}
}
\examples{
q<-0.5
beta<-2.5
x<-1:10
hrdiweibull(x, q, beta)
}
\keyword{distribution}
