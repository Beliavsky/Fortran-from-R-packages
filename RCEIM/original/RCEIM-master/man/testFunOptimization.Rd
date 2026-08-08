\name{testFunOptimization}
\alias{testFunOptimization}

\title{
1D test problem for RCEIM
}
\description{
An one-dimension problem for testing optimization methods. 

This function was created for demonstrating the RCEIM package. It has the form:

\if{latex}{\deqn{f(x) = \exp(-((x - 2)^2)) + 0.9 * \exp(-((x + 2)^2)) + 0.5 * \sin(8*x) + 0.25 * \cos(2*x)}}
\if{html}{f(x) = exp(-((x - 2)^2)) + 0.9 * exp(-((x + 2)^2)) + 0.5 * sin(8*x) + 0.25 * cos(2*x)}

}
\usage{
testFunOptimization(x)
}
\arguments{
  \item{x}{The point where the function is computed.}
}

\value{
The value of the function at x. 
}
\author{
Alberto Krone-Martins
}

\seealso{
\code{\link{testFunOptimization2d}}
}
\examples{
# Create a graphical representation of the problem with a line plot
dev.new()
xx <- seq(-10,10,by=0.01) 
plot(xx, testFunOptimization(xx), type="l", xlab="x", ylab="Value")
rm(list=c('xx'))
}

\keyword{ misc }
\keyword{ design }
