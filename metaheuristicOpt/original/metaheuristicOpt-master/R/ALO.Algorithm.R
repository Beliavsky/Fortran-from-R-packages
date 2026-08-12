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
#' This is the internal function that implements Ant Lion Optimizer
#' Algorithm. It is used to solve continuous optimization tasks.
#' Users do not need to call it directly,
#' but just use \code{\link{metaOpt}}.
#'
#' This algorithm was proposed by (Mirjalili, 2015). The Ant Lion Optimizer (ALO)
#' algorithm mimics the hunting mechanism of antlions in nature. Five main steps
#' of hunting prey such as the random walk of ants, building traps, entrapment of
#' ants in traps, catching preys, and re-building traps are implemented.
#'
#' In order to find the optimal solution, the algorithm follow the following steps.
#' \itemize{
#' \item Initialization: Initialize the first population of ants and antlions randomly,
#'       calculate the fitness of ants and antlions and find the best antlions as the
#'       elite (determined optimum).
#' \item Update Ants Position: Select an antlion using Roulette Whell then update ants
#'       position based on random walk around selected antlion and elite.
#'       Furthermore, calculate the fitness of all ants.
#' \item Replace an antlion with its corresponding ant, if it becomes fitter
#' \item Update elite if an antlion becomes fitter than the elite
#' \item Check termination criteria, if termination criterion is satisfied, return the
#'       elite as the optimal solution for given problem. Otherwise, back to Update Ants Position steps.
#'}
#'
#' @title Optimization using Ant Lion Optimizer
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
#' ## calculate the optimum solution using Ant Lion Optimizer
#' resultALO <- ALO(schewefels2.22, optimType="MIN", numVar,
#' numPopulation=20, maxIter=100, rangeVar)
#'
#' ## calculate the optimum value using schewefel's problem 2.22 function
#' optimum.value <- schewefels2.22(resultALO)
#'
#' @return \code{Vector [v1, v2, ..., vn]} where \code{n} is number variable
#'         and \code{vn} is value of \code{n-th} variable.
#'
#' @references
#' Seyedali Mirjalili, The Ant Lion Optimizer, Advances in Engineering Software,
#' Volume 83, 2015, Pages 80-98, ISSN 0965-9978,
#' https://doi.org/10.1016/j.advengsoft.2015.01.010
#' @export

ALO <- function(FUN, optimType="MIN", numVar, numPopulation=40, maxIter=500, rangeVar){
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

	# generate initial population of antlion and ant
	antlion <- generateRandom(numPopulation, dimension, lowerBound, upperBound)
	ant <- generateRandom(numPopulation, dimension, lowerBound, upperBound)

	# find the best position
	bestPos <- engine.ALO(FUN, optimType, maxIter, lowerBound, upperBound, antlion, ant)

	return(bestPos)
}

## support function for calculating best position with ALO algorithm
# @param FUN objective function
# @param optimType type optimization
# @param maxIter maximum number iteration
# @param lowerBound lower bound for each variable
# @param upperBound upper bound for each variable
# @param antlion population of antlion
# @param ant population of ant

engine.ALO <- function(FUN, optimType, maxIter, lowerBound, upperBound, antlion, ant){
	# calculate the antlion fitness
	antlionFitness <- calcFitness(FUN, optimType, antlion)
	antFitness <- c() # will count later in iteration process

	# sort antlion location based on fitness value
	index <- order(antlionFitness)
	antlionFitness <- sort(antlionFitness)
	antlion <- antlion[index,]

	# set the current best position (bestPos = elite)
	bestPos <- antlion[1,]
	FbestPos <- antlionFitness[1]

	# curve to plot
	curve <- c()
	progressbar <- txtProgressBar(min = 0, max = maxIter, style = 3)

	for (t in 1:maxIter){
		for (i in 1:nrow(ant)){
			# select an antlion by roulette whell selection
			roulette.index <- rouletteWhell(1/antlionFitness)
			# calculate random walk around the selected antlion
			RA <- randomWalk(maxIter, lowerBound, upperBound, antlion[roulette.index,], t)

			# calculate random walk around the elites (best antlion)
			RE <- randomWalk(maxIter, lowerBound, upperBound, bestPos, t)

			ant[i,] <- (RA[t,]+RE[t,])/2
		}

		for (i in 1:nrow(ant)){
			# check boundary and bring back the ant
			ant[i,] <- checkBound(ant[i,], lowerBound, upperBound)

			# check ant fitness
			antFitness[i] <- optimType*FUN(ant[i,])
		}

		# this process show how ant become fitter and antlion goes to ant position
		# to build the new pit
		doublePopulation <- rbind(antlion, ant)
		doubleFitness <- c(antlionFitness, antFitness)

		# sort the doubleFitness
		index <- order(doubleFitness)
		doubleFitness <- sort(doubleFitness)
		# sort the double popu
		doublePopulation <- doublePopulation[index,]

		# get the new antlion fitness and position
		antlionFitness <- doubleFitness[1:nrow(antlion)]
		antlion <- doublePopulation[1:nrow(antlion),]

		# update the best antlion
		if(antlionFitness[1] < FbestPos){
			bestPos <- antlion[1,]
			FbestPos <- antlionFitness[1]
		}
		# save the best fitness for iteration t
		curve[t] <- FbestPos

		setTxtProgressBar(progressbar, t)
	}

	close(progressbar)
	curve <- curve*optimType
	# plot(c(1:maxIter), curve, type="l", main="ALO", log="y", xlab="Number Iteration", ylab = "Best Fittness",
		                  # ylim=c(curve[which.min(curve)],curve[which.max(curve)]))
	return(bestPos)
}

## support function for doing random walk
# @param maxIter maximum number iteration
# @param lowerBound lower bound for each variable
# @param upperBound upper bound for each variable
# @param position the current position of antlion
# @param numIter number of iteration

randomWalk <- function(maxIter, lowerBound, upperBound, position, numIter){
	# check length lb and ub
	# if user only define one lb and ub, then repeat it until the dimension
	if(length(lowerBound)==1 & length(upperBound)==1){
		lowerBound <- rep(lowerBound,length(position))
		upperBound <- rep(upperBound,length(position))
	}

	# I is the ratio I defined by
	# I = 10^w * numIter/maxIter
	I <- 1

	if(numIter > maxIter*0.1){
		I <- 1+100*(numIter/maxIter)
	}

	if(numIter > maxIter*0.5){
		I <- 1+1000*(numIter/maxIter)
	}

	if(numIter > maxIter*0.75){
		I <- 1+10000*(numIter/maxIter)
	}

	if(numIter > maxIter*0.9){
		I <- 1+100000*(numIter/maxIter)
	}

	if(numIter > maxIter*0.95){
		I <- 1+1000000*(numIter/maxIter)
	}

	# decrease boundaries to converge towards antlion
	lowerBound <- lowerBound/I
	upperBound <- upperBound/I

	# move the interval of lb and ub around the antlion
	if(runif(1) < 0.5){
		lowerBound <- lowerBound + position
	}else{
		lowerBound <- -lowerBound + position
	}

	if(runif(1) < 0.5){
		upperBound <- upperBound + position
	}else{
		upperBound <- -upperBound + position
	}

	# create n random walk and normalize according to modified lb and ub
	result <- matrix(ncol=length(position), nrow=maxIter+1)
	for (i in 1:length(position)){
		X <- c(0, cumsum(2*(runif(maxIter)>0.5)-1))

	    # normalize the random walk position using min-max normalization
		a <- min(X)
	    b <- max(X)
	    c <- lowerBound[i]
	    d <- upperBound[i]
	    X_norm <- ((X-a)*(d-c))/(b-a)+c
	    result[,i] <- X_norm
	}
	return(result)
}
