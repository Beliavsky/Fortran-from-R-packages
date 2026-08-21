\name{predict.TAR}
\alias{predict.TAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Prediction based on  a fitted TAR model}
\description{
Predictions based on  a fitted TAR model. The errors are assumed to be
normally distributed. 
The predictive distributions are approximated by simulation.
}
\usage{
\method{predict}{TAR}(object, n.ahead = 1, n.sim = 1000,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{a fitted TAR model from the tar function}
  \item{n.ahead}{number of prediction steps ahead}
  \item{n.sim}{simulation size}
  \item{...}{other arguments; not used here but kept for consistency with 
  the generic method}
}
\value{
  \item{fit}{a vector of medians of the 1-step to n.ahead-step predictive
distributions}
  \item{pred.interval}{a matrix whose i-th row consists of the 2.5 and 97.5 precentiles of the i-step predictive distribution}
  \item{pred.matrix}{ a matrix whose j-th column consists of all simulated 
values from the j-step predictive distribution} 
}
\references{ "Time Series Analysis, with Applications in R" by J.D. Cryer and K.S. Chan}
\author{Kung-Sik Chan}
\seealso{ \code{\link{tar}}}
\examples{
data(prey.eq)
prey.tar.1=tar(y=log(prey.eq),p1=4,p2=4,d=3,a=.1,b=.9,print=TRUE)
set.seed(2357125)
pred.prey=predict(prey.tar.1,n.ahead=60,n.sim=1000)
yy=ts(c(log(prey.eq),pred.prey$fit),frequency=1,start=1)
plot(yy,type='n',ylim=range(c(yy,pred.prey$pred.interval)),ylab='Log Prey',
xlab=expression(t))
lines(log(prey.eq))
lines(window(yy, start=end(prey.eq)[1]+1),lty=2)
lines(ts(pred.prey$pred.interval[2,],start=end(prey.eq)[1]+1),lty=2)
lines(ts(pred.prey$pred.interval[1,],start=end(prey.eq)[1]+1),lty=2)
}
\keyword{methods}
