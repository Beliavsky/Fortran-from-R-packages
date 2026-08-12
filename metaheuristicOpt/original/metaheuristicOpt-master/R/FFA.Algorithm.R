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
#' This is the internal function that implements Firefly
#' Algorithm. It is used to solve continuous optimization tasks.
#' Users do not need to call it directly,
#' but just use \code{\link{metaOpt}}.
#'
#' This algorithm was proposed by (Yang, 2009).
#' The firefly algorithm (FFA) mimics the behavior of fireflies, which use
#' a kind of flashing light to communicate with other members of their species.
#' Since the intensity of the light of a single firefly diminishes with
#' increasing distance, the FFA is implicitly able to detect local solutions
#' on its way to the best solution for a given objective function.
#'
#' In order to find the optimal solution, the algorithm follow the following steps.
#' \itemize{
#' \item Initialization: Initialize the first population of fireflies randomly,
#'       calculate the fitness of fireflies and assumes fitness values as
#'       Light Intensity.
#' \item Update the firefly position based on the attractiveness. The firefly that have higher light
#'       intensity will tend to attract other fireflies. The attracted firefly will move based on
#'       the parameter that given by user.
#' \item Calculate the fitness and update the best firefly position.
#' \item Check termination criteria, if termination criterion is satisfied, return the
#'       best position as the optimal solution for given problem. Otherwise, back to Update firefly position steps.
#'}
#'
#' @title Optimization using Firefly Algorithm
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
#' @param B0 a positive integer to determine the attractiveness firefly at r=0. The default value is 1.
#'
#' @param gamma a positive integer to determine light absorption coefficient. The default value is 1.
#'
#' @param alphaFFA a positive integer to determine randomization parameter. The default value is 0.2.
#'
#' @importFrom graphics plot
#' @importFrom stats runif
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @seealso \code{\link{metaOpt}}
#'
#' @examples
#' ##################################
#' ## Optimizing the quartic with noise function
#' # define Quartic with noise function as objective function
#' quartic <- function(x){
#'     dim <- length(x)
#'     result <- sum(c(1:dim)*(x^4))+runif(1)
#'     return(result)
#' }
#'
#' ## Define parameter
#' B0 <- 1
#' gamma <- 1
#' alphaFFA <- 0.2
#' numVar <- 5
#' rangeVar <- matrix(c(-1.28,1.28), nrow=2)
#'
#' ## calculate the optimum solution using Firefly Algorithm
#' resultFFA <- FFA(quartic, optimType="MIN", numVar, numPopulation=20,
#'                  maxIter=100, rangeVar, B0, gamma, alphaFFA)
#'
#' ## calculate the optimum value using sphere function
#' optimum.value <- quartic(resultFFA)
#'
#' @return \code{Vector [v1, v2, ..., vn]} where \code{n} is number variable
#'         and \code{vn} is value of \code{n-th} variable.
#'
#' @references
#' X.-S. Yang, Firefly algorithms for multimodal optimization, in:
#' Stochastic Algorithms: Foundations and Applications, SAGA 2009,
#' Lecture Notes in Computer Sciences, Vol. 5792, pp. 169-178 (2009).
#' @export

FFA <- function(FUN, optimType="MIN", numVar, numPopulation=40, maxIter=500, rangeVar, B0=1, gamma=1, alphaFFA=0.2){
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

	# generate initial population of particle
	fireflies <- generateRandom(numPopulation, dimension, lowerBound, upperBound)

	# find the best particle position
	bestFirefly <- engineFFA(FUN, optimType, maxIter, lowerBound, upperBound, B0, gamma, alphaFFA, fireflies)
	# pilihan kedua parameter nya diganti jadi list

	return(bestFirefly)
}

## support function for calculating best position with PSO algorithm
# @param FUN objective function
# @param optimType type optimization
# @param maxIter maximum number iteration
# @param lowerBound lower bound for each variable
# @param upperBound upper bound for each variable
# @param B0 a positive integer to determine the attractiveness at r=0.
# @param gamma a positive integer to determine light absorption coefficient.
# @param alpha a positive integer to determine randomization parameter.
# @param firefly population of candidate solution

engineFFA <- function(FUN, optimType, maxIter, lowerBound, upperBound, B0, gamma, alpha, fireflies){
	curve <- c()
	# calculate the fitness and sort
	Light <- calcFitness(FUN, optimType, fireflies)
	## save the index order
	index <- order(Light)
	## sort Light
	Light <- sort(Light)
	## sort firefly in Light intensity order
	fireflies <- fireflies[index,]
	# end sort

	Best <- fireflies[1,]

	progressbar <- txtProgressBar(min = 0, max = maxIter, style = 3)
	for (t in 1:maxIter){
		for (i in 1:nrow(fireflies)) {
			for (j in 1:nrow(fireflies)) {
				if (Light[j] < Light[i]){
					# move firefly i towards j in d-dimension
					## calculate distance between i and j
					r <- sqrt(sum((fireflies[i,] - fireflies[j,])^2))

					## calculate new position for firefly i
					# newPos <- fireflies[i,] + B0*exp(-gamma*r*r) * (fireflies[j,] - fireflies[i,]) + alpha * (runif(ncol(fireflies))-0.5)

					# new pos based on Xin-She Yang program
					beta0 <- 1;
					beta <- (beta0-B0)*exp(-gamma*r^2)+B0;
					randomization <- alpha*(runif(ncol(fireflies))-0.5);
					fireflies[i,] <- fireflies[i,]*(1-beta) + fireflies[j,]*beta + randomization;

					# bring back the firefly if it goes outside search range
					fireflies[i,] <- checkBound(fireflies[i,], lowerBound, upperBound)
				}
			}
		}
		# calculate the fitness and sort
		Light <- calcFitness(FUN, optimType, fireflies)
		## save the index of best Light
		bestIndex <- which.max(Light)

		# update Best firefly
		Best <- fireflies[bestIndex,]
		# save best fitness for plot
		curve[t] <- Light[bestIndex]

		# next progress bar
		setTxtProgressBar(progressbar, t)
	}
	close(progressbar)
	curve <- curve*optimType
	# plot(c(1:maxIter), curve, type="l", main="FFA", log="y", xlab="Number Iteration", ylab = "Best Fittness",
		                  # ylim=c(curve[which.min(curve)],curve[which.max(curve)]))
	return(Best)
}
