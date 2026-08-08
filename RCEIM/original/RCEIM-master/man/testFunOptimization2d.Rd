\name{testFunOptimization2d}
\alias{testFunOptimization2d}

\title{
2D test problem for RCEIM
}
\description{
A two-dimensional problem for testing optimization methods. 

This function was created for demonstrating the RCEIM package. It has the form:

\if{latex}{\deqn{ f(x_1,x_2) = \frac{(x_1-4)^2 + (x_2+2)^2}{50} - \frac{(x_1+2)^2 + (x_2+4)^2}{90}}}
\if{latex}{\deqn{ -\exp(-(x_1-2)^2) - 0.9*\exp(-(x_2+2)^2) - 0.5*\sin(8*x_1) - 0.25*\cos(2*x_2)}}
\if{latex}{\deqn{ +0.25*\sin(\frac{x_1*x_2}{2}) + 0.5*\cos(\frac{x_2*x_1}{2.5})}}
\if{html}{
    f(x_1,x_2) = ((x_1-4)^2 + (x_2+2)^2)/50 - ((x_1+2)^2 + (x_2+4)^2)/90 
	-exp(-(x_1-2)^2) - 0.9*exp(-(x_2+2)^2) - 0.5*sin(8*x_1) - 0.25*cos(2*x_2)
	+0.25*sin(x_1*x_2/2) + 0.5*cos(x_2*x_1/2.5)}
}
\usage{
testFunOptimization2d(x)
}
\arguments{
  \item{x}{a vector with the point where the function is computed.}
}
\value{
The value of the function at the requested point (x_1, x_2).
}
\author{
Alberto Krone-Martins
}

\seealso{
\code{\link{testFunOptimization}}
}
\examples{
# Create a graphical representation of the problem with a contour plot
dev.new()
xx <- seq(-10,10,by=0.1)
yy <- seq(-10,10,by=0.1)
zz <- matrix(nrow=length(yy), ncol=length(xx))
for(i in 1:length(xx)){
	for(j in 1:length(yy)){
		zz[i,j] <- testFunOptimization2d( c(xx[i],yy[j]) )
	}
}
image(xx,yy,zz, col=gray((50:100)/100), xlab="x", ylab="y")
contour(xx,yy,zz, add=TRUE)
rm(list=c('xx','yy','zz'))
}
\keyword{ misc }
\keyword{ design }
