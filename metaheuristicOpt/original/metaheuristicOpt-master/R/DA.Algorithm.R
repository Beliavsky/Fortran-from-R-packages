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
#' This is the internal function that implements Dragonfly
#' Algorithm. It is used to solve continuous optimization tasks.
#' Users do not need to call it directly,
#' but just use \code{\link{metaOpt}}.
#'
#' This algorithm was proposed by (Mirjalili, 2015). The main inspiration of the
#' DA algorithm originates from the static and dynamic swarming behaviours of
#' dragonflies in nature. Two essential phases of optimization, exploration and
#' exploitation, are designed by modelling the social interaction of dragonflies
#' in navigating, searching for foods, and avoiding enemies when swarming
#' dynamically or statistically.
#'
#' In order to find the optimal solution, the algorithm follow the following steps.
#' \itemize{
#' \item Initialization: Initialize the first population of dragonflies randomly,
#'       calculate the fitness of dragonflies and find the best dragonfly as
#'       food source and the worst dragonfly as enemy position.
#' \item Calculating Behaviour Weight that affecting fly direction and distance.
#'       First, find the neighbouring dragonflies for each dragonfly then calculate the behaviour weight.
#'       The behaviour weight consist of separation, alignment, cohesion, attracted toward food sources
#'       and distraction from enemy. The neighbouring dragonfly determined by the neighbouring radius
#'       that increasing linearly for each iteration.
#' \item Update the position each dragonfly using behaviour weight and the delta (same as velocity in PSO).
#' \item Calculate the fitness and update food and enemy position
#' \item Check termination criteria, if termination criterion is satisfied, return the
#'       food position as the optimal solution for given problem. Otherwise, back to Calculating Behaviour Weight steps.
#'}
#'
#' @title Optimization using Dragonfly Algorithm
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
#' ## calculate the optimum solution using dragonfly algorithm
#' resultDA <- DA(schewefels1.2, optimType="MIN", numVar, numPopulation=20,
#'                  maxIter=100, rangeVar)
#'
#' ## calculate the optimum value using schewefel's problem 1.2 function
#' optimum.value <- schewefels1.2(resultDA)
#'
#' @return \code{Vector [v1, v2, ..., vn]} where \code{n} is number variable
#'         and \code{vn} is value of \code{n-th} variable.
#'
#' @references
#' Seyedali Mirjalili. 2015. Dragonfly algorithm: a new meta-heuristic optimization
#' technique for solving single-objective, discrete, and multi-objective problems.
#' Neural Comput. Appl. 27, 4 (May 2015), 1053-1073.
#' DOI=https://doi.org/10.1007/s00521-015-1920-1
#'
#' @export

DA <- function(FUN, optimType="MIN", numVar, numPopulation=40, maxIter=500, rangeVar){
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

	# generate initial population of dragonfly
	dragonfly <- generateRandom(numPopulation, dimension, lowerBound, upperBound)

	# find the best position
	bestPos <- engineDA(FUN, optimType, maxIter, lowerBound, upperBound, dragonfly)

	return(bestPos)
}

## support function for calculating best position with Dragonfly algorithm
# @param FUN objective function
# @param optimType type optimization
# @param maxIter maximum number iteration
# @param lowerBound lower bound for each variable
# @param upperBound upper bound for each variable
# @param dragonfly population of dragonfly

engineDA <- function(FUN, optimType, maxIter, lowerBound, upperBound, dragonfly){
	# check length lb and ub
	# if user only define one lb and ub, then repeat it until the dimension
	if(length(lowerBound)==1 & length(upperBound)==1){
		lowerBound <- rep(lowerBound,ncol(dragonfly))
		upperBound <- rep(upperBound,ncol(dragonfly))
	}

	# boundary of delta
	deltaMax <- (upperBound-lowerBound)/20
	delta <- generateRandom(nrow(dragonfly), ncol(dragonfly), -deltaMax, deltaMax)

	# calculate the dragonfly fitness
	dragonflyFitness <- calcFitness(FUN, optimType, dragonfly)

	# sort dragonfly location based on fitness value
	index <- order(dragonflyFitness)
	dragonflyFitness <- sort(dragonflyFitness)
	dragonfly <- dragonfly[index,]

	# set the current food position
	food <- dragonfly[1,]
	Ffood <- dragonflyFitness[1]

	# set current enemy position
	enemy <- dragonfly[nrow(dragonfly),]
	Fenemy <- dragonflyFitness[length(dragonflyFitness)]

	# curve to plot
	curve <- c()
	progressbar <- txtProgressBar(min = 0, max = maxIter, style = 3)

	# code is translated from MATLAB version
	# you can found the original
	for (t in 1:maxIter){
		# define neighbour distance
		# increasing for each iteration
		r <- (upperBound-lowerBound)/4+((upperBound-lowerBound)*(t/maxIter)*2)

		# w is inertia weight
		# decreased lineraly from 0.9 to 0.4
		w <- 0.9-t*((0.9-0.4)/maxIter)

		# my c decreasing from 0.1 to -0.1
		# this parameter is for define the weight of each behaviour in each iteration
		my_c <- 0.1-t*((0.1)/(maxIter/2))
		if(my_c<0){
			my_c <- 0
		}

		s <- 2*runif(1)*my_c # Seperation weight
	    a <- 2*runif(1)*my_c # Alignment weight
	    c <- 2*runif(1)*my_c # Cohesion weight
	    f <- 2*runif(1)      # Food attraction weight
	    e <- my_c        	 # Enemy distraction weight

		for (i in 1:nrow(dragonfly)){
			# first find the neighbours of i-th dragonfly
			index <- 1
			neighboursNo <- 0

			# save the neighbours delta and position
			neighboursDelta <- matrix(ncol=ncol(delta))
			neighboursDragonfly <- matrix(ncol=ncol(dragonfly))

			for (j in 1:nrow(dragonfly)){
				distance_ij <- euclideanDistance(dragonfly[i,], dragonfly[j,])
				if(all(distance_ij<=r) & all(distance_ij!=0)){
					neighboursNo <- neighboursNo+1
					if(index==1){
						neighboursDelta[index,] <- delta[j,]
						neighboursDragonfly[index,] <- dragonfly[j,]
					}else{
						neighboursDelta <- rbind(neighboursDelta, delta[j,])
						neighboursDragonfly <- rbind(neighboursDragonfly, dragonfly[j,])
					}
					index <- index+1
				}
			}

			# then count the behaviour of dragonfly
			# Separation
			S <- c(rep(0, ncol(dragonfly)))
			if(neighboursNo > 1){
				for (j in 1:neighboursNo){
					S <- S + (neighboursDragonfly[j,]-dragonfly[i,])
				}
				S <- -S
			}

			# Alignment
			A <- delta[i,]
			if(neighboursNo > 1){
				A <- colSums(neighboursDelta)/neighboursNo
			}else{

			}

			# Cohesion
			C_temp <- dragonfly[i,]
			if(neighboursNo > 1){
				C_temp <- colSums(neighboursDragonfly)/neighboursNo
			}
			C <- C_temp - dragonfly[i,]

			# attracted toward food
			dist2food <- euclideanDistance(dragonfly[i,], food)
			F <- c(rep(0, ncol(dragonfly)))
			if(all(dist2food <= r)){
				F <- food - dragonfly[i,]
			}

			# distracted from enemy
			dist2enemy <- euclideanDistance(dragonfly[i,], enemy)
			E <- c(rep(0, ncol(dragonfly)))
			if(all(dist2enemy <= r)){
				E <- enemy + dragonfly[i,]
			}

			if(any(dist2food>r)){
				if(neighboursNo > 1){
					for (j in 1:ncol(dragonfly)){
						delta[i,j] <- w*delta[i,j]+runif(1)*A[j]+runif(1)*C[j]+runif(1)*S[j]
						if(delta[i,j] > deltaMax[j]){
							delta[i,j] <- deltaMax[j]
						}
						if(delta[i,j] < -deltaMax[j]){
							delta[i,j] <- -deltaMax[j]
						}
						dragonfly[i,j] <- dragonfly[i,j]+delta[i,j]
					}
				}else{
					# using levy flight if there are no neighbouring dragonfly
					dragonfly[i,] <- dragonfly[i,]+t(Levy(ncol(dragonfly)))*dragonfly[i,]
					delta[i,] <- c(rep(0, ncol(delta)))
				}
			}else{
				for (j in 1:ncol(dragonfly)){
					delta[i,j] <- (a*A[j]+c*C[j]+s*S[j]+f*F[j]+e*E[j]) + w*delta[i,j]
					if(delta[i,j] > deltaMax[j]){
						delta[i,j] <- deltaMax[j]
					}
					if(delta[i,j] < -deltaMax[j]){
						delta[i,j] <- -deltaMax[j]
					}
					dragonfly[i,j] <- dragonfly[i,j]+delta[i,j]
				}
			}

			# bring back dragonfly if it go outside search space
			dragonfly[i,] <- checkBound(dragonfly[i,], lowerBound, upperBound)
		}

		for (i in 1:nrow(dragonfly)){
			fitness <- optimType*FUN(dragonfly[i,])
			# update food position
	        if(fitness<Ffood){
	            Ffood <- fitness
	            food <- dragonfly[i,]
	        }
	        # update enemy position
	        if(fitness>Fenemy){
	            Fenemy <- fitness
	            enemy <- dragonfly[i,]
	        }
	    }

		# save the best fitness for iteration t
		curve[t] <- Ffood

		setTxtProgressBar(progressbar, t)
	}

	close(progressbar)
	curve <- curve*optimType
	# plot(c(1:maxIter), curve, type="l", main="DA", log="y", xlab="Number Iteration", ylab = "Best Fittness",
		                  # ylim=c(curve[which.min(curve)],curve[which.max(curve)]))
	return(food)
}

# this function is for calculating the euclidean distance between two grashoper
# @param a is the vector position of frist dragonfly
# @param b is the vector position of second dragonfly
# @return Vector containing distance between two dragonfly

euclideanDistance <- function(a,b){
	result <- c()
	for (i in 1:length(a)) {
		result[i] <- sqrt((a[i]-b[i])^2)
	}
	return(result)
}

# this function is for creating the Lefy flight
# @param dim is a integer value that indicates the dimension of levy flight
# @return double value

Levy <- function(dim){
	beta <- 3/2

	sigma <- (gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta)
	u <- runif(dim)*sigma
	v <- runif(dim)
	step <- u/abs(v)^(1/beta)

	result <- 0.01*step
	return(result)
}
