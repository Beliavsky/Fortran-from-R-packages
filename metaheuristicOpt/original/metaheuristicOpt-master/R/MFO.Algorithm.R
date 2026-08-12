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
#' This is the internal function that implements Moth Flame Optimization
#' Algorithm. It is used to solve continuous optimization tasks.
#' Users do not need to call it directly,
#' but just use \code{\link{metaOpt}}.
#'
#' This algorithm was proposed (Mirjalili, 2015). The main inspiration of
#' this optimizer is the navigation method of moths in nature called transverse
#' orientation. Moths fly in night by maintaining a fixed angle with respect to
#' the moon, a very effective mechanism for travelling in a straight line
#' for long distances. However, these fancy insects are trapped in a useless/deadly
#' spiral path around artificial lights.
#'
#' In order to find the optimal solution, the algorithm follow the following steps.
#' \itemize{
#' \item Initialization: Initialize the first population of moth randomly,
#'       calculate the fitness of moth and find the best moth as the best flame obtained so far
#'       The flame indicate the best position obtained by motion of moth. So in this step, position of
#'       flame will same with the position of moth.
#' \item Update Moth Position: All moth move around the corresponding flame.
#'       In every iteration, the number flame is decreasing over the iteration.
#'       So at the end of iteration all moth will move around the best solution obtained so far.
#' \item Replace a flame with the position of moth if a moth becomes fitter than flame
#' \item Check termination criteria, if termination criterion is satisfied, return the
#'       best flame as the optimal solution for given problem. Otherwise, back to Update Moth Position steps.
#'}
#'
#' @title Optimization using Moth Flame Optimizer
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
#' ## Optimizing the schewefel's problem 2.22 function
#'
#' # define schewefel's problem 2.22 function as objective function
#' schewefels2.22 <- function(x){
#'    return(sum(abs(x)+prod(abs(x))))
#' }
#'
#' ## Define parameter
#' numVar <- 5
#' rangeVar <- matrix(c(-10,10), nrow=2)
#'
#' ## calculate the optimum solution using Moth Flame Optimizer
#' resultMFO <- MFO(schewefels2.22, optimType="MIN", numVar, numPopulation=20,
#'                  maxIter=100, rangeVar)
#'
#' ## calculate the optimum value using schewefel's problem 2.22 function
#' optimum.value <- schewefels2.22(resultMFO)
#'
#' @return \code{Vector [v1, v2, ..., vn]} where \code{n} is number variable
#'         and \code{vn} is value of \code{n-th} variable.
#'
#' @references
#' Seyedali Mirjalili, Moth-flame optimization algorithm: A novel nature-inspired
#' heuristic paradigm, Knowledge-Based Systems, Volume 89, 2015, Pages 228-249,
#' ISSN 0950-7051, https://doi.org/10.1016/j.knosys.2015.07.006
#'
#' @export

MFO <- function(FUN, optimType="MIN", numVar, numPopulation=40, maxIter=500, rangeVar){
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

	# generate initial population of moth
	moth <- generateRandom(numPopulation, dimension, lowerBound, upperBound)

	# find the best position
	bestPos <- engineMFO(FUN, optimType, maxIter, lowerBound, upperBound, moth)

	return(bestPos)
}

## support function for calculating best position with MFO algorithm
# @param FUN objective function
# @param optimType type optimization
# @param maxIter maximum number iteration
# @param lowerBound lower bound for each variable
# @param upperBound upper bound for each variable
# @param moth population of moth

engineMFO <- function(FUN, optimType, maxIter, lowerBound, upperBound, moth){
	# calculate the moth fitness
	mothFitness <- calcFitness(FUN, optimType, moth)

	# sort moth location based on fitness value
	index <- order(mothFitness)
	flameFitness <- sort(mothFitness)
	flames <- moth[index,]

	# set the current best position
	bestPos <- flames[1,]
	FbestPos <- flameFitness[1]

	# curve to plot
	curve <- c()
	progressbar <- txtProgressBar(min = 0, max = maxIter, style = 3)

	for (t in 1:maxIter){
		# number of flames
		flame_no <- round(nrow(moth)-t*((nrow(moth)-1)/maxIter))

		# store the previous moth position

		# value of a linearly decreased from -1 to -2
		a <- -1+t*((-1)/maxIter)

		for (i in 1:nrow(moth)){
			for (j in 1:ncol(moth)){
                dist2flame <- abs(flames[i,j]-moth[i,j])
                b <- 1
                # r is random number in [-1,1]
                r <- (a-1)*runif(1)+1

				if(i <= flame_no){
	                moth[i,j] <- dist2flame*exp(b*r)*cos(r*2*pi)+flames[i,j]
				}
				if(i > flame_no){
	                moth[i,j] <- dist2flame*exp(b*r)*cos(r*2*pi)+flames[flame_no,j]
				}
			}

			# check bound
			moth[i,] <- checkBound(moth[i,], lowerBound, upperBound)
			# calculate i-th moth fitness
			mothFitness[i] <- optimType*FUN(moth[i,])
		}

		# combine the previous moth with current sorted moth
		mothUnion <- rbind(moth, flames)
		fitnessUnion <- c(mothFitness, flameFitness)

		# then sort the union
		index <- order(fitnessUnion)
		fitnessUnion <- sort(fitnessUnion)
		mothUnion <- mothUnion[index,]

		# get N moth with best fitness
		# N is number of moth in one population
		flames <- mothUnion[1:nrow(moth),]
		flameFitness <- fitnessUnion[1:nrow(moth)]

		# update the best position
		bestPos <- flames[1,]
		FbestPos <- flameFitness[1]

		# save the best fitness for iteration t
		curve[t] <- FbestPos

		setTxtProgressBar(progressbar, t)
	}

	close(progressbar)
	curve <- curve*optimType
	# plot(c(1:maxIter), curve, type="l", main="MFO", log="y", xlab="Number Iteration", ylab = "Best Fittness",
		                  # ylim=c(curve[which.min(curve)],curve[which.max(curve)]))
	return(bestPos)
}
