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
#' This is the internal function that implements Grasshopper
#' Algorithm. It is used to solve continuous optimization tasks.
#' Users do not need to call it directly,
#' but just use \code{\link{metaOpt}}.
#'
#' Grasshopper Optimisation Algorithm (GOA) was proposed by (Mirjalili et al., 2017).
#'The algorithm mathematically models and mimics the
#' behaviour of grasshopper swarms in nature for solving optimisation
#' problems.
#'
#' @title Optimization using Grasshopper Optimisation Algorithm
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
#' ## Optimizing the schewefel's problem 1.2 function
#'
#' # define schewefel's problem 1.2 function as objective function
#' schewefels1.2 <- function(x){
#'   dim <- length(x)
#'   result <- 0
#'     for(i in 1:dim){
#'        result <- result + sum(x[1:i])^2
#'    }
#'   return(result)
#' }
#'
#' ## Define parameter
#' numVar <- 5
#' rangeVar <- matrix(c(-10,10), nrow=2)
#'
#' ## calculate the optimum solution using grasshoper algorithm
#' resultGOA <- GOA(schewefels1.2, optimType="MIN", numVar, numPopulation=20,
#'                  maxIter=100, rangeVar)
#'
#' ## calculate the optimum value using schewefel's problem 1.2 function
#' optimum.value <- schewefels1.2(resultGOA)
#'
#' @return \code{Vector [v1, v2, ..., vn]} where \code{n} is number variable
#'         and \code{vn} is value of \code{n-th} variable.
#'
#' @references
#' Shahrzad Saremi, Seyedali Mirjalili, Andrew Lewis, Grasshopper Optimisation
#' Algorithm: Theory and application, Advances in Engineering Software,
#' Volume 105, March 2017, Pages 30-47, ISSN 0965-9978,
#' https://doi.org/10.1016/j.advengsoft.2017.01.004
#' @export

GOA <- function(FUN, optimType="MIN", numVar, numPopulation=40, maxIter=500, rangeVar){
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

	# generate initial population of grasshopper
	grasshopper <- generateRandom(numPopulation, dimension, lowerBound, upperBound)

	# find the best position
	bestPos <- engineGOA(FUN, optimType, maxIter, lowerBound, upperBound, grasshopper)

	return(bestPos)
}

## support function for calculating best position with HS algorithm
# @param FUN objective function
# @param optimType type optimization
# @param maxIter maximum number iteration
# @param lowerBound lower bound for each variable
# @param upperBound upper bound for each variable
# @param grasshopper a matrix of grasshopper

engineGOA <- function(FUN, optimType, maxIter, lowerBound, upperBound, grasshopper){
	state <- 0
	dimension <- ncol(grasshopper)
	mlb <- mean(lowerBound)
	mub <- mean(upperBound)
	# check length lb and ub
	# if user only define one lb and ub, then repeat it until the dimension
	if(length(lowerBound)==1 & length(upperBound)==1){
		lowerBound <- rep(lowerBound,dimension)
		upperBound <- rep(upperBound,dimension)
	}

	# this algo will calculate distance between two grashopper, so it should run with even number of dimension
	if (dimension %% 2 != 0){
		dimension <- dimension+1
		upperBound <- c(upperBound, mub)
		lowerBound <- c(lowerBound, mlb)
		grasshopper <- generateRandom(nrow(grasshopper), dimension, lowerBound, upperBound)
		# state indicates that one dimension is added
		state <- 1
	}

	# calculate the grasshopper fitness
	grasshopperFitness <- c()
	if(state == 1){
		grasshopperFitness <- calcFitness(FUN, optimType, grasshopper[,1:dimension-1])
	}else{
		grasshopperFitness <- calcFitness(FUN, optimType, grasshopper)
	}

	# sort the fitness
	index <- order(grasshopperFitness)
	grasshopperFitness <- sort(grasshopperFitness)

	# set the current best position
	bestPos <- grasshopper[index[1],]
	FbestPos <- grasshopperFitness[1]

	Cmax <- 1
	Cmin <- 0.00004

	# curve to plot
	curve <- c()
	progressbar <- txtProgressBar(min = 0, max = maxIter, style = 3)

	for (t in 1:maxIter){
		# balancing the c values for exploration and exploitation
		c <- Cmax-t*((Cmax-Cmin)/maxIter)

		newGrasshopper <- grasshopper
		for (i in 1:nrow(grasshopper)) {
			tempG <- t(grasshopper)
			Si_total <- matrix(ncol=1, nrow=dimension)
			for (j in seq(from=1, to=dimension, by=2)) {
				Si <- matrix(c(0,0), ncol=1, nrow=2)
				for (k in 1:nrow(grasshopper)) {
					if(i != k){
						# distance between two point
						nextJ <- j+1;
						R <- distance(tempG[j:nextJ,k], tempG[j:nextJ,i])

						part3 <- (tempG[j:nextJ,k]-tempG[j:nextJ,i])/(R+2.2204e-16)
						xj_xi <- 2+ R%%2

						s_ij <- ((upperBound[j:nextJ] - lowerBound[j:nextJ])*c/2)*S(xj_xi)*part3
						Si <- Si+s_ij
					}
				}
				Si_total[j] <- Si[1]
				Si_total[j+1] <- Si[2]
			}
			newX <- c * t(Si_total) + bestPos
			newGrasshopper[i,] <- newX
		}

		grasshopper <- newGrasshopper

		for (i in 1:nrow(grasshopper)) {
			grasshopper[i,] <- checkBound(grasshopper[i,], lowerBound, upperBound)
			if(state == 1){
				grasshopperFitness[i] <- optimType*FUN(grasshopper[i,1:dimension-1])
			}else{
				grasshopperFitness[i] <- optimType*FUN(grasshopper[i,])
			}

			if(grasshopperFitness[i] < FbestPos){
				bestPos <- grasshopper[i,]
				FbestPos <- grasshopperFitness[i]
			}
		}

		# save the best fitness for iteration t
		curve[t] <- FbestPos

		setTxtProgressBar(progressbar, t)
	}

	if(state==1){
		bestPos <- bestPos[1:dimension-1]
	}

	close(progressbar)
	curve <- curve*optimType
	# plot(c(1:maxIter), curve, type="l", main="GOA", log="y", xlab="Number Iteration", ylab = "Best Fittness",
		                  # ylim=c(curve[which.min(curve)],curve[which.max(curve)]))
	return(bestPos)
}

## support function for calculating distance
# @param a is vector first position
# @param b is vector second position
distance <- function(a, b){
	result <- sqrt(sum((a-b)^2))
	return(result)
}

## support function for calculating S function
# @param R double

S <- function(R){
	f <- 0.5
	l <- 1.5
	o <- f*exp(-R/l)-exp(-R)
	return(o)
}
