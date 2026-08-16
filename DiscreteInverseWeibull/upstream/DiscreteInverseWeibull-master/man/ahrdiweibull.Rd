\name{ahrdiweibull}
\alias{ahrdiweibull}
\title{
Alternative hazard rate function
}
\description{
Alternative hazard rate function for the discrete inverse Weibull distribution
}
\usage{
ahrdiweibull(x, q, beta)
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
the value of the alternative hazard rate function in the \code{x} values
}
\details{
The alternative hazard rate function is defined as \eqn{h(x)=\log(P(X>x-1)/P(X>x))=\log[(1-q^{(x-1)^{-\beta}})/(1-q^{x^{-\beta}})]}
}
\seealso{
\code{\link{hrdiweibull}}
}
\examples{
q<-0.5
beta<-2
x<-1:10
y<-ahrdiweibull(x, q, beta)
y
plot(x,y,ylab="alt.hazard rate")
}
\keyword{distribution}
