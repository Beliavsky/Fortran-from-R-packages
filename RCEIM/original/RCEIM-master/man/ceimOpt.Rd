\name{ceimOpt}
\alias{ceimOpt}
\title{
A Cross Entropy Inspired Method for Optimization
}
\description{
This is a Cross-Entropy Inspired Method for optimization.
}
\usage{
ceimOpt(OptimFunction = "testFunOptimization", nParams = 1, minimize = TRUE, 
    Ntot = 1000, N_elite = floor(Ntot/4), N_super = 1, alpha = 1, epsilon = 0.1, 
    q = 2, maxIter = 50, waitGen = maxIter, boundaries = t(matrix(rep(c(-10, 10), 
    nParams), ncol = nParams)),plotConvergence = FALSE, chaosGen = maxIter, 
    handIterative = FALSE, verbose = FALSE, plotResultDistribution = FALSE, 
    parallelVersion = FALSE)
}
\arguments{
  \item{OptimFunction}{
A string with the name of the function that will be optimized.
}
  \item{nParams}{
An integer with the number of parameters of the function that will be optimized.
}
  \item{minimize}{
A boolean indicating if the OptimFunction should be minimized or maximized.
}
  \item{Ntot}{
An integer with the number of individuals per iteration.
}
  \item{N_elite}{
An integer with the number of elite individuals, or in other words, the individuals used to define the individuals of the new iteration.
}
  \item{N_super}{
An integer with the number of super-individuals, or those individuals with the best fitness values, that are directly replicated to the next iteration.
}
  \item{alpha}{
A parameter of the CE method used to control the convergence rate, and to prevent early convergence.
}
  \item{epsilon}{
A convergence control parameter: if the maximum st.dev. of the parameters of the elite individuals divided by its average value is smaller than this number, the method considers that it converged.
}
  \item{q}{
A parameter of the CE method used to control the convergence rate, and to prevent early convergence.  
}
  \item{maxIter}{
The maximum number of iterations that the method will run before stop.
}
  \item{waitGen}{
The number of iterations that the method will wait: after "waitGen" without any improvement in the best individual, the method gives up and return the best individual as an answer.
}
  \item{boundaries}{
A matrix with as many rows as there are parameters and two columns the first column stores the minimum value, while the second, the maximum.
}
  \item{plotConvergence}{
A flag to indicate if the user wants to visually check the convergence of the method.
}
  \item{chaosGen}{
The number of iterations before the method replaces all the solutions, but the super-individuals, by a new random trial.
}
  \item{handIterative}{
A flag to indicate if the user wants to press enter between the each generation.
}
  \item{verbose}{
A flag to indicate if the user wants to receive some convergence and distribution information printed on the screen.
}
  \item{plotResultDistribution}{
A flag to indicate if the user wants to see the resulting distribution of elite members (black curve), the value of the fittest member (red line) and of the average member (blue line), for each parameter.
}
  \item{parallelVersion}{
A flag to indicate if the user wants to use all the cores in his/her computer to compute the fitness functions.
}
}
\details{
This is a simple stochastic heuristic method for optimization. It starts from a random population of points, computes the values of the function and selects a fraction of the points - the elite members. Then, based on these fittest points, it constructs a gaussian distribution, samples new random points from it, and iterates until convergence is reached (this is controlled by the epsilon parameter) or other stopping criteria is fulfilled (such as the maximum number of iterations).

The method does not impose strong conditions on the function to be optimized. The function must written as an R function, but it does not need to be neither continuous, differentiable or analytic. Moreover, the method is ready for multicore processing, enabling the optimization of time-consuming functions.
}
\value{
A list that contains:
\item{BestMember }{The parameters and the fitness value of the best member.}
\item{Convergence }{A boolean indicating if the method reached convergence.}
\item{Criteria }{Stopping criterion.}
\item{Iterations }{The amount of iterations.}
\item{EliteMembers }{The parameters and fitness values of the elite members at the last iteration.}
}
\author{
Alberto Krone-Martins
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

\keyword{ optimize }