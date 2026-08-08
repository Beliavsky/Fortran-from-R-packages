\name{RCEIM-package}
\alias{RCEIM-package}
\alias{RCEIM}
\docType{package}
\title{
An R Cross Entropy Inspired Method for Optimization
RCEIM
}
\description{
RCEIM is a package implementing a stochastic heuristic method for performing multidimensional function optimization. The method is inspired in the Cross-Entropy Method.
}
\details{
\tabular{ll}{
Package: \tab RCEIM\cr
Type: \tab Package\cr
Version: \tab 0.3\cr
Date: \tab 2015-02-19\cr
License: \tab GPL (>= 2)\cr
}
RCEIM implements a simple stochastic heuristic method for optimization in the function \code{\link{ceimOpt}}. This method starts from a random population of solutions, computes the values of the function and selects a fraction of these solutions - the elite members. Afterwards, based on the elite members it creates a gaussian distribution, samples new random solutions from it, and iterates until convergence is reached (this is controlled by an epsilon parameter) or other stopping criteria is fulfilled (such as the maximum number of iterations).

One advantage of this method is that it does not impose strong conditions on the function to be optimized. The function must written as an R function, but it does not need to be continuous, differentiable, neither analytic. Moreover, the method is ready for multicore processing, enabling the optimization of time-consuming functions.

\if{html}{\figure{RCEIM-Problems.jpg}{RCEIM-Problems}}
\if{latex}{\figure{RCEIM-Problems.pdf}{options: width=4.5in}}

Two examples of 1D and 2D functions that can be used as test problems for RCEIM (defined in \code{\link{testFunOptimization}} and \code{\link{testFunOptimization2d}}) are represented in the above figures.

}
\author{
Alberto Krone-Martins, Laurent Gallucio

Maintainer: Alberto Krone-Martins <algol@sim.ul.pt>
}
\keyword{ package }
\keyword{ optimize }
\seealso{
\code{\link{ceimOpt}}
}
\examples{
# Solve a simple optimization problem and shows the results
po <- ceimOpt(OptimFunction=function(x){(x[1]+1)^2+(x[2]+2)^2}, maxIter=100, epsilon=0.3, 
      nParams=2)
plotEliteDistrib(po$EliteMembers)
rm(po)

# A harder problem in 1D
po <- ceimOpt(OptimFunction="testFunOptimization", maxIter=10, epsilon=0.3, 
      nParams=1, verbose=TRUE)
dev.new()
xx <- seq(-10,10,by=0.01)
plot(xx, testFunOptimization(xx), type="l", xlab="x", ylab="Value")
points(po$BestMember[1], po$BestMember[2], col="red")
rm(list=c('xx','po'))

# A harder problem in 2D
po <- ceimOpt(OptimFunction="testFunOptimization2d", maxIter=20, epsilon=0.3, 
      nParams=2, verbose=TRUE)
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
points(po$BestMember[1], po$BestMember[2], col="red", pch=19, cex=0.5)
rm(list=c('xx','yy','zz'))

# Example of multicore processing
slowFunction <- function(x) { 
  p<-runif(50000)
  return((x+3)^2)
} 
system.time(po <- ceimOpt(OptimFunction="slowFunction", maxIter=10,
    Ntot=100, epsilon=0.3, nParams=1, verbose=FALSE, parallel=FALSE))
print(po$BestMember)
system.time(po <- ceimOpt(OptimFunction="slowFunction", maxIter=10,
    Ntot=100, epsilon=0.3, nParams=1, verbose=FALSE, parallel=TRUE))
print(po$BestMember)
rm(po)
}