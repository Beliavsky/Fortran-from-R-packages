\name{multi_mcga}
\alias{multi_mcga}

\title{
Performs multi objective machine coded genetic algorithms.
}
\description{
Machine coded genetic algorithm (MCGA) is a fast tool for real-valued optimization problems. It uses the byte representation of variables rather than real-values. It performs the classical crossover operations (uniform) on these byte representations. Mutation operator is also similar to classical mutation operator, which is to say, it changes a randomly selected byte value of a chromosome by +1 or -1 with probability 1/2. In MCGAs there is no need for encoding-decoding process and the classical operators are directly applicable on real-values. It is fast and can handle a wide range of a search space with high precision. Using a 256-unary alphabet is the main disadvantage of this algorithm but a moderate size population is convenient for many problems.

This function performs multi objective optimization using the same logic underlying the mcga. Chromosomes are sorted by their 
objective values using a non-dominated sorting algorithm. 
}
\usage{
multi_mcga(popsize, chsize, crossprob = 1.0, mutateprob = 0.01, 
		   elitism = 1, minval, maxval, maxiter = 10, numfunc, evalFunc)
}
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
The maximum number of generations. By default it is 10.
}
  \item{numfunc}{
Number of objective functions.
}
  \item{evalFunc}{
An R function. By default, each problem is a minimization. This function must return a cost vector
with dimension of numfunc. Each element of this vector points to the corresponding function to optimize. 
}
}
\value{
\item{population}{Sorted population resulted after generations}
\item{costs}{Cost values for each chromosomes in the resulted population}
\item{ranks}{Calculated ranks using a non-dominated sorting for each chromosome}
}
\author{
Mehmet Hakan Satman - mhsatman@istanbul.edu.tr
}

\references{
Deb, K. (2000). An efficient constraint handling method for
     genetic algorithms. Computer methods in applied mechanics and
     engineering, 186(2), 311-338.
}


\examples{
\dontrun{
 # We have two objective functions.
 f1<-function(x){
   return(sin(x))
 }

 f2<-function(x){
   return(sin(2*x))
 }

 # This function returns a vector of cost functions for a given x sent from mcga
 f<-function(x){
   return ( c( f1(x), f2(x)) )
 }

 # main loop
 m<-multi_mcga(popsize=200, chsize=1, minval= 0, elitism=2, 
 	      maxval= 2.0 * pi, maxiter=1000, crossprob=1.0, 
	      mutateprob=0.01, evalFunc=f, numfunc=2)

 # Points show best five solutions. 
 curve(f1, 0, 2*pi)
 curve(f2, 0, 2*pi, add=TRUE)

 p <- m$population[1:5,]
 points(p, f1(p))
 points(p, f2(p))
}
}
