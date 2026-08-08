\name{mcga-package}
\alias{mcga-package}
\docType{package}
\title{
Machine Coded Genetic Algorithms for Real-valued Optimization Problems
}
\description{
Machine coded genetic algorithm (MCGA) is a fast tool for real-valued optimization problems. It uses the byte representation of variables rather than real-values. It performs the classical crossover operations (uniform) on these byte representations. Mutation operator is also similar to classical mutation operator, which is to say, it changes a randomly selected byte value of a chromosome by +1 or -1 with probability 1/2. In MCGAs there is no need for encoding-decoding process and the classical operators are directly applicable on real-values. It is fast and can handle a wide range of a search space with high precision. Using a 256-unary alphabet is the main disadvantage of this algorithm but a moderate size population is convenient for many problems.
}

\author{
Mehmet Hakan Satman

Maintainer: Mehmet Hakan Satman <mhsatman@istanbul.edu.tr>

}
\examples{
\dontrun{
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
}
