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
#' This is the internal function that implements Sine Cosine
#' Algorithm. It is used to solve continuous optimization tasks.
#' Users do not need to call it directly,
#' but just use \code{\link{metaOpt}}.
#'
#' This algorithm was proposed by (Mirjalili, 2016). The SCA creates multiple initial
#' random candidate solutions and requires them to fluctuate outwards or towards the
#' best solution using a mathematical model based on sine and cosine functions. Several
#' random and adaptive variables also are integrated to this algorithm to emphasize
#' exploration and exploitation of the search space in different milestones of optimization.
#'
#' In order to find the optimal solution, the algorithm follow the following steps.
#' \itemize{
#' \item Initialization: Initialize the first population of candidate solution randomly,
#'       calculate the fitness of candidate solution and find the best candidate.
#' \item Update Candidate Position: Update the position with the equation that represent the
#'       behaviour of sine and cosine function.
#' \item Update the best candidate if there are candidate solution with better fitness.
#' \item Check termination criteria, if termination criterion is satisfied, return the
#'       best candidate as the optimal solution for given problem. Otherwise, back to Update Candidate Position steps.
#'}
#'
#' @title Optimization using Sine Cosine Algorithm
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
#' @importFrom graphics plot
#' @importFrom stats runif
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @seealso \code{\link{metaOpt}}
#'
#' @examples
#' ##################################
#' ## Optimizing the step function
#'
#' # define step function as objective function
#' step <- function(x){
#'     result <- sum(abs((x+0.5))^2)
#'     return(result)
#' }
#'
#' ## Define parameter
#' numVar <- 5
#' rangeVar <- matrix(c(-100,100), nrow=2)
#'
#' ## calculate the optimum solution using Sine Cosine Algorithm
#' resultSCA <- SCA(step, optimType="MIN", numVar, numPopulation=20,
#'                  maxIter=100, rangeVar)
#'
#' ## calculate the optimum value using step function
#' optimum.value <- step(resultSCA)
#'
#' @return \code{Vector [v1, v2, ..., vn]} where \code{n} is number variable
#'         and \code{vn} is value of \code{n-th} variable.
#'
#' @references
#' Seyedali Mirjalili, SCA: A Sine Cosine Algorithm for solving optimization problems,
#' Knowledge-Based Systems, Volume 96, 2016, Pages 120-133, ISSN 0950-7051,
#' https://doi.org/10.1016/j.knosys.2015.12.022
#'
#' @export

SCA <- function(FUN, optimType="MIN", numVar, numPopulation=40, maxIter=500, rangeVar){
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
	bestPos <- engine.SCA(FUN, optimType, maxIter, lowerBound, upperBound, candidate)

	return(bestPos)
}

## support function for calculating best position with SCA algorithm
# @param FUN objective function
# @param optimType type optimization
# @param maxIter maximum number iteration
# @param lowerBound lower bound for each variable
# @param upperBound upper bound for each variable
# @param candidate population of candidate

engine.SCA <- function(FUN, optimType, maxIter, lowerBound, upperBound, candidate){
	# calculate the candidate fitness
	candidateFitness <- calcFitness(FUN, optimType, candidate)

	# sort candidate location based on fitness value
	index <- order(candidateFitness)
	candidateFitness <- sort(candidateFitness)
	candidate <- candidate[index,]

	# set the current best position
	bestPos <- candidate[1,]
	FbestPos <- candidateFitness[1]

	# curve to plot
	curve <- c()
	progressbar <- txtProgressBar(min = 0, max = maxIter, style = 3)

	for (t in 1:maxIter){
		# value a in eq (3.4)
		a <- 2

		# value r1 decreased linearly from a to 0
		r1 <- a-t*((a)/maxIter)

		for (i in 1:nrow(candidate)){
			for (j in 1:ncol(candidate)) {
	            # generate random number for each dimension
	            r2 <- (2*pi)*runif(1)
	            r3 <- 2*runif(1)
	            r4 <- runif(1)

	            if(r4 < 0.5){
	            	candidate[i,j] <- candidate[i,j]+(r1*sin(r2)*abs(r3*bestPos[j]-candidate[i,j]))
	            }else{
	            	candidate[i,j] <- candidate[i,j]+(r1*cos(r2)*abs(r3*bestPos[j]-candidate[i,j]))
	            }

			}

			# bring back candidate if it go outside search space
			candidate[i,] <- checkBound(candidate[i,], lowerBound, upperBound)

			fitness <- optimType*FUN(candidate[i,])

			# update bestPos
	        if(fitness<FbestPos){
	            FbestPos <- fitness
	            bestPos <- candidate[i,]
	        }
		}

		# save the best fitness for iteration t
		curve[t] <- FbestPos

		setTxtProgressBar(progressbar, t)
	}

	close(progressbar)
	curve <- curve*optimType
	# plot(c(1:maxIter), curve, type="l", main="SCA", log="y", xlab="Number Iteration", ylab = "Best Fittness",
		                  # ylim=c(curve[which.min(curve)],curve[which.max(curve)]))
	return(bestPos)
}
