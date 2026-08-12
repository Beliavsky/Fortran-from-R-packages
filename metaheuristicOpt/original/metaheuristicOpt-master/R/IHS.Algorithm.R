#############################################################################
#
#  This file is a part of the R package "metaheuristicOpt".
#
#  Author: Iip
#  Co-author: -
#  Supervisors: Lala Septem Riza, Eddy Prasetyo Nugroho
#
#
#  This package is free software: you can redistribute it and/or modify it under
#  the terms of the GNU General Public License as published by the Free Software
#  Foundation, either version 2 of the License, or (at your option) any later version.
#
#  This package is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
#  A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
#############################################################################
#' This is the internal function that implements Improved Harmony Search
#' Algorithm. It is used to solve continuous optimization tasks.
#' Users do not need to call it directly,
#' but just use \code{\link{metaOpt}}.
#'
#' Harmony Search (HS)  was proposed by (Geem et al., 2001)
#' mimicking the improvisation of music players. Furthermore,
#' Improved Harmny Search (HS), proposed by Mahdavi, employs a method for
#' generating new solution vectors that enhances accuracy and convergence
#' rate of harmony search algorithm.
#'
#' In order to find the optimal solution, the algorithm follow the following steps.
#' \itemize{
#' \item Step 1. Initialized the problem and algorithm parameters
#' \item Step 2. Initialize the Harmony Memory, creating the Harmony memory and give
#'       random rumber for each memory.
#' \item Step 3. Improvise new Harmony, Generating new Harmony based on parameter defined by user
#' \item Step 4. Update the Harmony Memory, If new harmony have better fitness than the worst harmony in
#'       Harmony Memory, then replace the worst harmony with new Harmony.
#' \item Step 5. Check termination criteria, if termination criterion is satisfied, return the
#'       best Harmony as the optimal solution for given problem. Otherwise, back to Step 3.
#'}
#'
#' @title Optimization using Harmony Search Algorithm
#'
#' @param FUN an objective function or cost function,
#'
#' @param optimType a string value that represent the type of optimization.
#'        There are two option for this arguments: \code{"MIN"} and \code{"MAX"}.
#'        The default value is \code{"MIN"}, which the function will do minimization.
#'        Otherwise, you can use \code{"MAX"} for maximization problem.
#'        The default value is \code{"MIN"}.
#'
#' @param numVar a positive integer to determine the number variables.
#'
#' @param numPopulation a positive integer to determine the number populations. The default value is 40.
#'
#' @param maxIter a positive integer to determine the maximum number of iterations. The default value is 500.
#'
#' @param PAR a positive integer to determine the value of Pinch Adjusting Ratio. The default value is 0.3.
#'
#' @param HMCR a positive integer to determine the Harmony Memory Consideration Rate. The default value is 0.95.
#'
#' @param bandwith a positive integer to determine the bandwith. The default value is 0.05.
#'
#' @param rangeVar a matrix (\eqn{2 \times n}) containing the range of variables,
#'        where \eqn{n} is the number of variables, and first and second rows
#'        are the lower bound (minimum) and upper bound (maximum) values, respectively.
#'        If all variable have equal upper bound, you can define \code{rangeVar} as
#'        matrix (\eqn{2 \times 1}).
#'
#' @importFrom graphics plot
#' @importFrom stats runif
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @seealso \code{\link{metaOpt}}
#'
#' @examples
#' ##################################
#' ## Optimizing the quartic with noise function
#'
#' # define Quartic with noise function as objective function
#' quartic <- function(x){
#'     dim <- length(x)
#'     result <- sum(c(1:dim)*(x^4))+runif(1)
#'     return(result)
#' }
#'
#' ## Define parameter
#' numVar <- 5
#' rangeVar <- matrix(c(-10,10), nrow=2)
#' PAR <- 0.3
#' HMCR <- 0.95
#' bandwith <- 0.05
#'
#' ## calculate the optimum solution using Harmony Search algorithm
#' resultHS <- HS(quartic, optimType="MIN", numVar, numPopulation=20,
#'                  maxIter=100, rangeVar, PAR, HMCR, bandwith)
#'
#' ## calculate the optimum value using quartic with noise function
#' optimum.value <- quartic(resultHS)
#'
#' @return \code{Vector [v1, v2, ..., vn]} where \code{n} is number variable
#'         and \code{vn} is value of \code{n-th} variable.
#'
#' @references
#' Geem, Zong Woo, Joong Hoon Kim, and G. V. Loganathan (2001). "A new
#' heuristic optimization algorithm: harmony search." Simulation 76.2: pp. 60-68.
#'
#' M. Mahdavi, M. Fesanghary, E. Damangir, An improved harmony search algorithm
#' for solving optimization problems, Applied Mathematics and Computation,
#' Volume 188, Issue 2, 2007, Pages 1567-1579, ISSN 0096-3003,
#' https://doi.org/10.1016/j.amc.2006.11.033
#' @export

HS <- function(FUN, optimType="MIN", numVar, numPopulation=40, maxIter=500, rangeVar, PAR=0.3, HMCR=0.95, bandwith=0.05){
	# calculate the dimension of problem if not specified by user
	dimension <- ncol(rangeVar)

	# parsing rangeVar to lowerBound and upperBound
	lowerBound <- rangeVar[1,]
	upperBound <- rangeVar[2,]

	# if user define the same upper bound and lower bound for each dimension
	if(dimension==1){
		dimension <- numVar
	}

	## convert optimType to numerical form
	## 1 for minimization and -1 for maximization
	if(optimType == "MAX") optimType <- -1 else optimType <- 1

	# generate initial population of harmonyMemory
	harmonyMemory <- generateRandom(numPopulation, dimension, lowerBound, upperBound)

	# find the best position
	bestPos <- engineHS(FUN, optimType, maxIter, lowerBound, upperBound, PAR, HMCR, bandwith, harmonyMemory)

	return(bestPos)
}

## support function for calculating best position with HS algorithm
# @param FUN objective function
# @param optimType type optimization
# @param maxIter maximum number iteration
# @param lowerBound lower bound for each variable
# @param upperBound upper bound for each variable
# @param harmonyMemory a matrix of harmonyMemory

engineHS <- function(FUN, optimType, maxIter, lowerBound, upperBound, PAR, HMCR, bandwith, harmonyMemory){
	# check length lb and ub
	# if user only define one lb and ub, then repeat it until the dimension
	if(length(lowerBound)==1 & length(upperBound)==1){
		lowerBound <- rep(lowerBound,ncol(harmonyMemory))
		upperBound <- rep(upperBound,ncol(harmonyMemory))
	}

	# calculate the harmonyMemory fitness
	harmonyMemoryFitness <- calcFitness(FUN, optimType, harmonyMemory)

	# set only for initialization
	bestPos <- harmonyMemory[which.min(harmonyMemoryFitness),]
	FbestPos <- harmonyMemoryFitness[which.min(harmonyMemoryFitness)]

	# curve to plot
	curve <- c()
	progressbar <- txtProgressBar(min = 0, max = maxIter, style = 3)


	for (t in 1:maxIter){
		# number of new harmonies = number population
		newHarmonies <- matrix(ncol=ncol(harmonyMemory), nrow=nrow(harmonyMemory))

		# looping to number population to create new harmony population
		for (j in 1:nrow(harmonyMemory)){
			newHarmony <- c()
			for (i in 1:ncol(harmonyMemory)){
				# check probability to get value from harmony memory
				if(runif(1) < HMCR){
					# get new position from random harmony memory
					index <- runif(1, 1, nrow(harmonyMemory))
					newHarmony[i] <- harmonyMemory[index,i]

					# check probability to adjusting pich
					if(runif(1) < PAR){
						# adjust pich based on random values
						r <- runif(1)
						if(r < 0.5){
							temp <- newHarmony[i] - r * bandwith
						}else{
							temp <- newHarmony[i] + r * bandwith
						}

						# check bound
						if(temp >= lowerBound[i] & temp <= upperBound[i]){
							newHarmony[i] <- temp
						}
					}
				}else{
					# get new position from random value
					newHarmony[i] <- runif(1, lowerBound[i], upperBound[i])
				}
			}
			newHarmonies[j,] <- newHarmony
		}

		#calc fitness of new harmonies
		newFitness <- calcFitness(FUN, optimType, newHarmonies)

		# merge HM with new harmonies
		doubleMemory <- rbind(harmonyMemory, newHarmonies)
		doubleFitness <- c(harmonyMemoryFitness, newFitness)

		# sort the doubleFitness
		index <- order(doubleFitness)
		doubleFitness <- sort(doubleFitness)
		# sort the double memory
		doubleMemory <- doubleMemory[index,]

		# get the new harmony memory fitness and position
		harmonyMemoryFitness <- doubleFitness[1:nrow(harmonyMemory)]
		harmonyMemory <- doubleMemory[1:nrow(harmonyMemory),]

		# update the best harmony
		if(harmonyMemoryFitness[1] < FbestPos){
			bestPos <- harmonyMemory[1,]
			FbestPos <- harmonyMemoryFitness[1]
		}

		# save the best fitness for iteration t
		curve[t] <- FbestPos

		setTxtProgressBar(progressbar, t)
	}

	close(progressbar)
	curve <- curve*optimType
	# plot(c(1:maxIter), curve, type="l", main="IHS", log="y", xlab="Number Iteration", ylab = "Best Fittness",
		                  # ylim=c(curve[which.min(curve)],curve[which.max(curve)]))
	return(bestPos)
}
