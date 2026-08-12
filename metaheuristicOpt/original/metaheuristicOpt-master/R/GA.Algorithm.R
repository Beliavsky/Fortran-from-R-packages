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
#' This is the internal function that implements Genetic
#' Algorithm. It is used to solve continuous optimization tasks.
#' Users do not need to call it directly,
#' but just use \code{\link{metaOpt}}.
#'
#' Genetic algorithms (GA) were invented by John Holland in the 1960 and
#' were developed by Holland and his students and colleagues at the
#' University of Michigan in the 1960 and the 1970. GA are commonly used
#' to generate high-quality solutions to optimization and search problems
#' by relying on bio-inspired operators such as mutation, crossover and
#' selection.
#'
#' In order to find the optimal solution, the algorithm follow the following steps.
#' \itemize{
#' \item Initialization: Initialize the first population randomly,
#'       calculate the fitness and save the best fitness as bestPopulation.
#' \item Selection: Select set of individual parent for doing crossover. Number of parent
#'       determined by the crossover probability which defined by user.
#'       In this work, we use method called Roulette Whell Selection.
#' \item Crossover: Doing crossover between two parent from Selection step.
#'       This step done by selecting two point randomly and switching the values between them.
#' \item Mutation : All individu in population have a chance to mutate. When mutation occurs, we
#'       generate the random values to replace the old one.
#' \item Calculate the fitness of each individual and update bestPopulation.
#' \item Check termination criteria, if termination criterion is satisfied, return the
#'       bestPopulation as the optimal solution for given problem. Otherwise, back to Selection steps.
#'}
#'
#' @title Optimization using Genetic Algorithm
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
#' @param rangeVar a matrix (\eqn{2 \times n}) containing the range of variables,
#'        where \eqn{n} is the number of variables, and first and second rows
#'        are the lower bound (minimum) and upper bound (maximum) values, respectively.
#'        If all variable have equal upper bound, you can define \code{rangeVar} as
#'        matrix (\eqn{2 \times 1}).
#'
#' @param Pm a positive integer to determine mutation probability. The default value is 0.1.
#'
#' @param Pc a positive integer to determine crossover probability. The default value is 0.8.
#'
#' @importFrom graphics plot
#' @importFrom stats runif
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @seealso \code{\link{metaOpt}}
#'
#' @examples
#' ##################################
#' ## Optimizing the sphere function
#'
#' # define sphere function as objective function
#' sphere <- function(x){
#'     return(sum(x^2))
#' }
#'
#' ## Define parameter
#' Pm <- 0.1
#' Pc <- 0.8
#' numVar <- 5
#' rangeVar <- matrix(c(-10,10), nrow=2)
#'
#' ## calculate the optimum solution using Genetic Algorithm
#' resultGA <- GA(sphere, optimType="MIN", numVar, numPopulation=20,
#'                  maxIter=100, rangeVar, Pm, Pc)
#'
#' ## calculate the optimum value using sphere function
#' optimum.value <- sphere(resultGA)
#'
#' @return \code{Vector [v1, v2, ..., vn]} where \code{n} is number variable
#'         and \code{vn} is value of \code{n-th} variable.
#'
#' @references
#' Holland, J. H. 1975. Adaptation in Natural and Artificial Systems.
#' University of Michigan Press. (Second edition: MIT Press, 1992.)
#'
#' Melanie Mitchell. 1998. An Introduction to Genetic Algorithms.
#' MIT Press, Cambridge, MA, USA.
#' @export

GA <- function(FUN, optimType="MIN", numVar, numPopulation=40, maxIter=500, rangeVar, Pm=0.1, Pc=0.8){
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

	# generate initial population of candidate
	candidate <- generateRandom(numPopulation, dimension, lowerBound, upperBound)

	# find the best position
	bestPos <- engineGA(FUN, optimType, maxIter, lowerBound, upperBound, Pm, Pc, candidate)

	return(bestPos)
}

## support function for calculating best position with HS algorithm
# @param FUN objective function
# @param optimType type optimization
# @param maxIter maximum number iteration
# @param lowerBound lower bound for each variable
# @param upperBound upper bound for each variable
# @param candidate a matrix of candidate solution

engineGA <- function(FUN, optimType, maxIter, lowerBound, upperBound, Pm, Pc, candidate){
	# check length lb and ub
	# if user only define one lb and ub, then repeat it until the dimension
	if(length(lowerBound)==1 & length(upperBound)==1){
		lowerBound <- rep(lowerBound,ncol(candidate))
		upperBound <- rep(upperBound,ncol(candidate))
	}

	# calculate the candidate fitness
	candidateFitness <- calcFitness(FUN, optimType, candidate)

	# set current best
	bestPos <- candidate[which.min(candidateFitness),]
	FbestPos <- candidateFitness[which.min(candidateFitness)]

	# curve to plot
	curve <- c()
	progressbar <- txtProgressBar(min = 0, max = maxIter, style = 3)

	for (t in 1:maxIter){
		# do selection to determine candidate to do crossover phase
		numSelect <- Pc * nrow(candidate)
		# even numSelect
		if (numSelect %% 2 != 0) numSelect <- numSelect + 1
		# select index with machine roulette selection
		index <- c()
		while(length(index) < numSelect){
			temp <- rouletteWhell(candidateFitness)
			index <- c(index,temp)
		}
		# now we have parent for doing crossover
		parent <- candidate[index,]

		# do two point crossover
		offspring <- matrix(ncol=ncol(candidate),nrow=numSelect)
		for (i in seq(from=1, to=numSelect, by=2)) {
			offspring[i,] <- parent[i,]
			offspring[i+1,] <- parent[i+1,]
			# select two point
			p1 <- sample(c(1:ncol(candidate)),1)
			p2 <- sample(c(1:ncol(candidate)),1)
			# swap element
			if(p1 < p2){
				offspring[i,p1:p2] <- parent[i+1,p1:p2]
				offspring[i+1,p1:p2] <- parent[i,p1:p2]
			}else{
				offspring[i,p2:p1] <- parent[i+1,p2:p1]
				offspring[i+1,p2:p1] <- parent[i,p2:p1]
			}
		}

		# do mutation
		newPopulation <- rbind(candidate, offspring)
		for (i in 1:nrow(newPopulation)) {
			# pool number [0,1]
			p <- runif(1)
			if(p < Pm){
				# do mutation
				n <- ceiling(0.1 * ncol(candidate))
				for (j in 1:n) {
					index <- sample(c(1:ncol(candidate)), 1)
					newPopulation[i,index] <- runif(1, lowerBound[index], upperBound[index])
				}
			}
		}

		newFitness <- calcFitness(FUN, optimType, newPopulation)

		# sort new location based on fitness value
		index <- order(newFitness)
		newFitness <- sort(newFitness)
		newPopulation <- newPopulation[index,]

		candidate <- newPopulation[1:nrow(candidate),]
		candidateFitness <- newFitness[1:nrow(candidate)]

		bestPos <- candidate[1,]
		FbestPos <- candidateFitness[1]

		# save the best fitness for iteration t
		curve[t] <- FbestPos

		setTxtProgressBar(progressbar, t)
	}

	close(progressbar)
	curve <- curve*optimType
	# plot(c(1:maxIter), curve, type="l", main="GA", log="y", xlab="Number Iteration", ylab = "Best Fittness",
		                  # ylim=c(curve[which.min(curve)],curve[which.max(curve)]))
	return(bestPos)
}
