#############################################################################
#
#  This file is a part of the R package "metaheuristicOpt".
#
#  Author: Muhammad Bima Adi Prabowo
#  Co-author: -
#  Supervisors: Lala Septem Riza, Enjun Junaeti
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
#' This is the internal function that implements Artificial Bee Colony
#' Algorithm. It is used to solve continuous optimization tasks.
#' Users do not need to call it directly,
#' but just use \code{\link{metaOpt}}.
#'
#' This algorithm was proposed by (Karaboga & Akay, 2009).
#' It inspired by type of bee. They are three types of bee employeed,
#' onlooker and scout. Employed bee work by finding food source.
#' Onlooker bee work by finding better food source other than foods
#' that Employed bee found. Scout bee work by removing abandoned food
#' source. Each candidate solution in ABC algorithm represent as bee
#' and they will move in 3 phases employed, onlooker and scout.
#'
#' In order to find the optimal solution, the algorithm follow the following steps.
#' \itemize{
#' \item initialize population randomly.
#' \item Employed bee phase (Perform local search and greedy algorithm for each candidate solution).
#' \item Onlooker bee phase (Perform local search and greedy algorithm for some candidate solutions).
#' \item Scout bee phase (Remove abandoned candidate solutions).
#' \item If a termination criterion (a maximum number of iterations or a sufficiently good fitness) is met,
#'       exit the loop, else back to employed bee phase.
#' }
#'
#' @title Optimization using Artificial Bee Colony Algorithm
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
#' @param cycleLimit a positive integer to determine number of times allowed for
#'        candidate solution to not move. The default value is as.integer(numVar * numPopulation).
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
#' ## calculate the optimum solution using artificial bee colony algorithm
#' resultABC <- ABC(sphere, optimType="MIN", numVar, numPopulation=20,
#'                  maxIter=100, rangeVar)
#'
#' ## calculate the optimum value using sphere function
#' optimum.value <- sphere(resultABC)
#'
#' @return \code{Vector [v1, v2, ..., vn]} where \code{n} is number variable
#'         and \code{vn} is value of \code{n-th} variable.
#' @references
#' Karaboga, D., & Akay, B. (2009). A comparative study of artificial bee colony algorithm.
#' Applied mathematics and computation, 214(1), 108-132.
#'
#' @export
# Artificial Bee Colony Algorithm (ABC)

ABC <- function(FUN, optimType="MIN", numVar, numPopulation=40, maxIter=500, rangeVar, cycleLimit=as.integer(numVar*numPopulation)){
  # Validation
  if(numPopulation < 1){
    stop("numPopulation must greater than 0")
  }

  if(maxIter < 0){
    stop("maxIter must greater than or equal to 0")
  }

  if(cycleLimit < 0){
    stop("cycleLimit must greater than 0")
  }else if(!is.integer(cycleLimit)){
    stop("cycleLimit must be integer (as.integer())")
  }
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

  # if user only define one lb and ub, then repeat it until the dimension
  if(length(lowerBound)==1 & length(upperBound)==1){
    lowerBound <- rep(lowerBound, dimension)
    upperBound <- rep(upperBound, dimension)
  }

  # generate candidate solution
  candidateSolution <- generateRandomABC(numPopulation, numVar, min(lowerBound), max(upperBound))
  bestPos <- engineABC(FUN, optimType, maxIter, lowerBound, upperBound, candidateSolution, cycleLimit)

  return(bestPos)
}

engineABC <- function(FUN, optimType, maxIter, lowerBound, upperBound, candidateSolution, cycleLimit){
  numVar <- ncol(candidateSolution)
  numPopulation <- nrow(candidateSolution)
  limit <- rep(cycleLimit, numPopulation)

  progressbar <- txtProgressBar(min = 0, max = maxIter, style = 3)
  for(t in 1:maxIter){
    # Employeed Bee Phase
    i <- 1:numPopulation
    j <- sample(1:numVar, numPopulation, replace = TRUE)
    k <- sample(1:numPopulation, numPopulation, replace = TRUE)
    while(!all((k == i) == FALSE)){
      k <- sample(1:numPopulation, numPopulation, replace = TRUE)
    }
    Xij <- matrix(rep(0, numVar*numPopulation), ncol = numVar)
    Xkj <- Xij
    Xij[cbind(i, j)] <- candidateSolution[cbind(i, j)]
    Xkj[cbind(i, j)] <- candidateSolution[cbind(k, j)]
    randomMatrix <- matrix(runif(numPopulation*numVar, min = -1, max = 1), ncol = numVar)
    new <- candidateSolution + randomMatrix *(Xij - Xkj)
    newFitness <- calcFitness(FUN, optimType, new)
    CSFitness <- calcFitness(FUN, optimType, candidateSolution)
    applyGreedy <- newFitness <= CSFitness
    candidateSolution[applyGreedy,] <- new[applyGreedy,]
    limit[!applyGreedy] <- limit[!applyGreedy] - 1

    # Onlooker Bee Phase
    randomVector <- runif(numPopulation)
    CSFitness <- calcFitness(FUN, optimType, candidateSolution)
    # prevent division by 0
    if(sum(CSFitness) == 0){
      prob <- rep(0, numPopulation)
    }else{
      prob <- CSFitness/sum(CSFitness)
    }
    isOnlooker <- randomVector < prob
    if(!all(isOnlooker == F)){
      onLooker <- 1:numPopulation
      onLooker <- onLooker[isOnlooker]
      i <- onLooker
      j <- sample(1:numVar, length(i), replace = TRUE)
      k <- sample(1:numPopulation, length(i), replace = TRUE)
      while(!all((k == i) == FALSE)){
        k <- sample(1:numPopulation, length(i), replace = TRUE)
      }
      Xij <- matrix(rep(0, numVar*length(i)), ncol = numVar)
      Xkj <- Xij
      Xij[cbind(1:length(i), j)] <- candidateSolution[cbind(i, j)]
      Xkj[cbind(1:length(i), j)] <- candidateSolution[cbind(k, j)]
      randomMatrix <- matrix(runif(length(onLooker)*numVar, min = -1, max = 1), ncol = numVar)
      new <- candidateSolution[i,] + randomMatrix *(Xij - Xkj)
      newFitness <- calcFitness(FUN, optimType, new)
      CSFitness <- calcFitness(FUN, optimType, matrix(candidateSolution[i,], ncol = numVar))
      applyGreedy <- newFitness <= CSFitness
      allApplyGreedy <- rep(F, numPopulation)
      allApplyGreedy[i] <- applyGreedy
      candidateSolution[allApplyGreedy,] <- new[applyGreedy,]
      limit[i[!applyGreedy]] <- limit[i[!applyGreedy]] - 1
    }
    # Scout Bee Phase
    CSFitness <- calcFitness(FUN, optimType, matrix(candidateSolution[i,], ncol = numVar))
    best <- order(CSFitness)[1]
    limit[best] <- cycleLimit
    isAbandoned <- rep(T, numPopulation)
    isAbandoned[limit > 0] <- F
    if(!all(isAbandoned == F)){
      candidateSolution[isAbandoned] <- generateRandomABC(length(isAbandoned[isAbandoned == T]), numVar, min(lowerBound), max(upperBound))
    }

    # check bound
    for(i in 1:nrow(candidateSolution)){
      candidateSolution[i,] <- checkBound(candidateSolution[i,], lowerBound, upperBound)
    }

    setTxtProgressBar(progressbar, t)
  }
  close(progressbar)
  return(calcBest(FUN, -1*optimType, candidateSolution))
}

generateRandomABC <- function(numPopulation, dimension, lowerBound, upperBound){
  matrixLowerBound <- matrix(rep(lowerBound, numPopulation * dimension), ncol = dimension, byrow = TRUE)
  matrixUpperBound <- matrix(rep(upperBound, numPopulation * dimension), ncol = dimension, byrow = TRUE)
  matrixRandom <- matrix(runif(numPopulation * dimension), ncol = dimension, byrow = TRUE)
  result <- matrixLowerBound + matrixRandom * (matrixUpperBound - matrixLowerBound)
  return(result)
}
