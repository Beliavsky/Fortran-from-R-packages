\name{mcga}
\alias{mcga}

\title{
Performs machine coded genetic algorithms on a function subject to be minimized.
}
\description{
Machine coded genetic algorithm (MCGA) is a fast tool for real-valued optimization problems. It uses the byte representation of variables rather than real-values. It performs the classical crossover operations (uniform) on these byte representations. Mutation operator is also similar to classical mutation operator, which is to say, it changes a randomly selected byte value of a chromosome by +1 or -1 with probability 1/2. In MCGAs there is no need for encoding-decoding process and the classical operators are directly applicable on real-values. It is fast and can handle a wide range of a search space with high precision. Using a 256-unary alphabet is the main disadvantage of this algorithm but a moderate size population is convenient for many problems.
}
\usage{
mcga(popsize, chsize, crossprob = 1.0, mutateprob = 0.01, 
	 elitism = 1, minval, maxval, maxiter = 10, evalFunc)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{popsize}{
Number of chromosomes. 
}
  \item{chsize}{
Number of parameters.
}
  \item{crossprob}{
Crossover probability. By default it is 1.0
}
  \item{mutateprob}{
Mutation probability. By default it is 0.01
}
  \item{elitism}{
Number of best chromosomes to be copied directly into next generation. By default it is 1
}
  \item{minval}{
The lower bound of the randomized initial population. This is not a constraint for parameters.
}
  \item{maxval}{
The upper bound of the randomized initial population. This is not a constraint for parameters.
}
  \item{maxiter}{
The maximum number of generations. By default it is 10
}
  \item{evalFunc}{
An R function. By default, each problem is a minimization. 
}
}
\value{
\item{population}{Sorted population resulted after generations}
\item{costs}{Cost values for each chromosomes in the resulted population}
}
\author{
Mehmet Hakan Satman - mhsatman@istanbul.edu.tr
}
\references{
M.H.Satman (2013), Machine Coded Genetic Algorithms for Real Parameter Optimization Problems, Gazi University Journal of Science, Vol 26, No 1, pp. 85-95
}

\examples{
# A sample optimization problem
# Min f(xi) = (x1-7)^2 + (x2-77)^2 + (x3-777)^2 + (x4-7777)^2 + (x5-77777)^2
# The range of xi is unknown. The solution is
# x1 = 7
# x2 = 77
# x3 = 777
# x4 = 7777
# x5 = 77777
# Min f(xi) = 0
require("mcga")
 f<-function(x){
    return ((x[1]-7)^2 + (x[2]-77)^2 +(x[3]-777)^2 +(x[4]-7777)^2 +(x[5]-77777)^2)
 }
 m <- mcga(	popsize=200, 
			chsize=5, 
			minval=0.0, 
			maxval=999999999.9, 
			maxiter=2500, 
			crossprob=1.0, 
			mutateprob=0.01, 
			evalFunc=f)
			
 cat("Best chromosome:\n")
 print(m$population[1,])
 cat("Cost: ",m$costs[1],"\n")
}
