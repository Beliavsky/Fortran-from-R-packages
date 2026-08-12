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
#' This is the internal function that implements Whale Optimization
#' Algorithm. It is used to solve continuous optimization tasks.
#' Users do not need to call it directly,
#' but just use \code{\link{metaOpt}}.
#'
#' This algorithm was proposed by (Mirjalili, 2016), which mimics the
#' social behavior of humpback whales. The algorithm is inspired by the
#' bubble-net hunting strategy.
#'
#' In order to find the optimal solution, the algorithm follow the following steps.
#' \itemize{
#' \item Initialization: Initialize the first population of whale randomly,
#'       calculate the fitness of whale and find the best whale position as the
#'       best position obtained so far.
#' \item Update Whale Position: Update the whale position using bubble-net hunting
#'       strategy. The whale position will depend on the best whale position obtained so far.
#'       Otherwise random whale choosen if the specific condition meet.
#' \item Update the best position if there are new whale that have better fitness
#' \item Check termination criteria, if termination criterion is satisfied, return the
#'       best position as the optimal solution for given problem. Otherwise, back to Update Whale Position steps.
#'}
#'
#' @title Optimization using Whale Optimization Algorithm
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
#' ## Optimizing the sphere function
#'
#' # define sphere function as objective function
#' sphere <- function(x){
#'     return(sum(x^2))
#' }
#'
#' ## Define parameter
#' numVar <- 5
#' rangeVar <- matrix(c(-10,10), nrow=2)
#'
#' ## calculate the optimum solution using Ant Lion Optimizer
#' resultWOA <- WOA(sphere, optimType="MIN", numVar, numPopulation=20,
#'                  maxIter=100, rangeVar)
#'
#' ## calculate the optimum value using sphere function
#' optimum.value <- sphere(resultWOA)
#'
#' @return \code{Vector [v1, v2, ..., vn]} where \code{n} is number variable
#'         and \code{vn} is value of \code{n-th} variable.
#'
#' @references
#' Seyedali Mirjalili, Andrew Lewis, The Whale Optimization Algorithm,
#' Advances in Engineering Software, Volume 95, 2016, Pages 51-67,
#' ISSN 0965-9978, https://doi.org/10.1016/j.advengsoft.2016.01.008
#'
#' @export

WOA <- function(FUN, optimType="MIN", numVar, numPopulation=40, maxIter=500, rangeVar){
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

	# generate initial population of whale
	whale <- generateRandom(numPopulation, dimension, lowerBound, upperBound)

	# find the best position
	bestPos <- engineWOA(FUN, optimType, maxIter, lowerBound, upperBound, whale)

	return(bestPos)
}

## support function for calculating best position with SCA algorithm
# @param FUN objective function
# @param optimType type optimization
# @param maxIter maximum number iteration
# @param lowerBound lower bound for each variable
# @param upperBound upper bound for each variable
# @param whale population of whale

engineWOA <- function(FUN, optimType, maxIter, lowerBound, upperBound, whale){
	# calculate the whale fitness
	whaleFitness <- calcFitness(FUN, optimType, whale)

	# sort whale location based on fitness value
	index <- order(whaleFitness)
	whaleFitness <- sort(whaleFitness)
	whale <- whale[index,]

	# set the current best position
	bestPos <- whale[1,]
	FbestPos <- whaleFitness[1]

	# curve to plot
	curve <- c()
	progressbar <- txtProgressBar(min = 0, max = maxIter, style = 3)

	for (t in 1:maxIter){
		# value a decreased linearly from 2 to 0
		a <- 2-t*((2)/maxIter)

		# value a2 decreased linearly from -1 to -2
		a2 <- -1+t*((-1)/maxIter)

		for (i in 1:nrow(whale)){
			# generate random number [0,1]
			r1 <- runif(1)
			r2 <- runif(1)

			# vector A and C
			A <- 2*a*r1-a
            C <- 2*r2

            # parameter b is constant for defining the shape of the logaritmic spiral
            # param l is random number in [-1,1]
            b <- 1
            l <- (a2-1)*runif(1)+1

            # p is random number to define the probability to select
            # shrinking encircling mechanism or spiral model
            p <- runif(1)

			for (j in 1:ncol(whale)) {
	            if(p < 0.5){
	            	if(abs(A) >= 1){
	            		# do exploration phase (search for prey)
	            		## choose random index of whale
	            		rand.index <- floor(nrow(whale)*runif(1)+1)
	            		x.rand <- whale[rand.index,]
	            		D.x.rand <- abs(C*x.rand[j]-whale[i,j])
	            		whale[i,j] <- x.rand[j]-A*D.x.rand
	            	}else if(abs(A) < 1){
	            		# do encircling prey
	            		D.bestPos <- abs(C*bestPos[j]-whale[i,j])
	            		whale[i,j] <- bestPos[j]-A*D.bestPos
	            	}
	            }else if(p >= 0.5){
	            	# distance whale to the prey
	            	distance <- abs(bestPos[j]-whale[i,j])
	            	whale[i,j] <- distance*exp(b*l)*cos(l*2*pi)+bestPos[j]
	            }

			}

			# bring back whale if it go outside search space
			whale[i,] <- checkBound(whale[i,], lowerBound, upperBound)

			fitness <- optimType*FUN(whale[i,])

			# update bestPos
	        if(fitness<FbestPos){
	            FbestPos <- fitness
	            bestPos <- whale[i,]
	        }
		}

		# save the best fitness for iteration t
		curve[t] <- FbestPos

		setTxtProgressBar(progressbar, t)
	}

	close(progressbar)
	curve <- curve*optimType
	# plot(c(1:maxIter), curve, type="l", main="WOA", log="y", xlab="Number Iteration", ylab = "Best Fittness",
		                  # ylim=c(curve[which.min(curve)],curve[which.max(curve)]))
	return(bestPos)
}
