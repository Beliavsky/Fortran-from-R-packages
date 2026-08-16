\name{estdiweibull}
\alias{estdiweibull}
\title{
Estimation of parameters
}
\description{
Sample estimation of the parameters of the discrete inverse Weibull distribution
}
\usage{
estdiweibull(x, method="P", control=list())
}
\arguments{
  \item{x}{
a vector of sample values
}
  \item{method}{
the estimation method that will be carried out: \code{"P"} method of proportion, \code{"M"} method of moments, \code{"H"} heuristic-maximum likelihood method, \code{"PP"} graphical method-probability plot
}
  \item{control}{
a list of additional parameters:
\code{eps}, \code{nmax} for the method of moments;
\code{beta1}, \code{z}, \code{r}, \code{Leps} for the heuristic method
}
}
\details{For a description of the methods, have a look at the reference. Note that they may be not applicable to some specific samples. For examples, the method of proportion cannot be applied if there are no 1s in the samples; it cannot be applied for estimating \eqn{\beta} if all the sample values are \eqn{\leq 2}. The method of moments cannot be applied for estimating \eqn{\beta} if all the sample values are \eqn{\leq 2}; besides, it may return unreliable results since the first and second moments can be computed only if \eqn{\beta>2}. The heuristic method cannot be applied for estimating \eqn{\beta} if all the sample values are \eqn{\leq 2}.
}
\value{
a vector containing the two estimates of \eqn{q} and \eqn{\beta}
}
\seealso{
\code{\link{heuristic}}, \code{\link{Ediweibull}}
}
\examples{
n<-100 
q<-0.5 
beta<-2.5
# generation of a sample
x<-rdiweibull(n, q, beta)
# sample estimation through each of the implemented methods
estdiweibull(x, method="P") 
estdiweibull(x, method="M")
estdiweibull(x, method="H")
estdiweibull(x, method="PP")
}
\keyword{distribution, htest}
