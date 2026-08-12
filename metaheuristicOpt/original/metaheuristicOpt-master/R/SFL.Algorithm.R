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
#' This is the internal function that implements Shuffled Frog Leaping
#' Algorithm. It is used to solve continuous optimization tasks.
#' Users do not need to call it directly,
#' but just use \code{\link{metaOpt}}.
#'
#' This algorithm was proposed by (Eusuff, Lansey & Pasha, 2006).
#' The main inspiration for SFL algorithm originates from how swarm of frogs
#' finding foods.
#'
#' In order to find the optimal solution, the algorithm follow the following steps.
#' \itemize{
#' \item initialize population randomly.
#' \item separate population into "numMemeplex" memeplexes.
#' \item update worst candidate solution using best candidate solution on
#'       each memeplex as much as "frogLeaping Iteration".
#' \item Shuffled back each memeplexes into population.
#' \item Sort population based on fitness.
#' \item If a termination criterion (a maximum number of iterations or a sufficiently good fitness) is met,
#'       exit the loop, else back to separate population into memeplexes.
#' }
#'
#' @title Optimization using Shuffled Frog Leaping Algorithm
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
#' @param numMemeplex a positive integer (as.integer()) between 0 and numVar to
#'        determine number of memeplexes (see details). The default value is \code{as.integer(numPopulation/3)}.
#'
#' @param frogLeapingIteration a positive integer (as.integer()) to determine number
#'        of iterations for each memeplex. The default value is \code{as.integer(10)}.
#'
#' @importFrom graphics plot
#' @importFrom stats runif
#' @importFrom utils setTxtProgressBar txtProgressBar head tail
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
#' rangeVar <- matrix(c(-1.28, 1.28), nrow=2)
#'
#' ## calculate the optimum solution shuffled frog leaping algorithm
#' resultSFL <- SFL(quartic, optimType="MIN", numVar, numPopulation=20,
#'                  maxIter=100, rangeVar)
#'
#' ## calculate the optimum value using quartic with noise function
#' optimum.value <- quartic(resultSFL)
#'
#' @return \code{Vector [v1, v2, ..., vn]} where \code{n} is number variable
#'         and \code{vn} is value of \code{n-th} variable.
#' @references
#' Eusuff, M., Lansey, K., & Pasha, F. (2006). Shuffled frog-leaping algorithm:
#' a memetic meta-heuristic for discrete optimization. Engineering Optimization,
#' 38(2), 129–154.
#'
#' @export
# Shuffled Frog Leaping -----

SFL <- function(FUN, optimType="MIN", numVar, numPopulation=40, maxIter=500, rangeVar,
                numMemeplex=as.integer(numPopulation/3), frogLeapingIteration=as.integer(10)){
  if(numPopulation < 1){
    stop("numPopulation must greater than 0")
  }

  if(!is.integer(numMemeplex)){
    stop("numMemeplex must be integer (as.integer())")
  }

  if(numMemeplex > numPopulation){
    stop("numMemeplex must less than or equal to numPopulation")
  }
  if(numMemeplex < 1){
    stop("numMemeplex must greater than 0")
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
  # initialize candidate solution
  candidateSolution <- generateRandom(numPopulation, dimension, lowerBound, upperBound)
  bestPos <- engineSFL(FUN, optimType, maxIter, lowerBound, upperBound, candidateSolution,
                       numMemeplex, frogLeapingIteration)
  return(bestPos)
}

engineSFL <- function(FUN, optimType, maxIter, lowerBound, upperBound, candidateSolution,
                      numMemeplex, frogLeapingIteration){
  numVar <- ncol(candidateSolution)
  numPopulation <- nrow(candidateSolution)
  CS <- candidateSolution
  gbest <- calcBest(FUN, -1*optimType, CS)
  memeplexId <- rep(1:numMemeplex, length.out = numPopulation)
  memeplexId <- memeplexId[order(memeplexId)]

  # set index
  index <- 1:numPopulation
  bestId <- c()
  worstId <- c()
  for(i in 1:numMemeplex){
    bestId <- c(bestId, head(index[memeplexId == i], n=1))
    worstId <- c(worstId, tail(index[memeplexId == i], n=1))
  }

  progressbar <- txtProgressBar(min = 0, max = maxIter, style = 3)
  for(t in 1:maxIter){
    # reshuffle
    memeplexId <- rep(1:numMemeplex, length.out = numPopulation)
    CSFitness <- calcFitness(FUN, 1, CS)
    CS <- matrix(CS[order(CSFitness),],ncol = numVar)
    CS <- matrix(CS[order(memeplexId),],ncol = numVar)
    memeplexId <- memeplexId[order(memeplexId)]
    for(i in 1:frogLeapingIteration){
      # update worst using best
      isItMove <- rep(F, length(worstId))
      randomMatrix <- matrix(runif(length(worstId) * numVar), ncol = numVar)
      new <- matrix(CS[worstId,], ncol = numVar) + randomMatrix * (CS[bestId,] - CS[worstId,])
      worstFitness <- calcFitness(FUN, optimType, matrix(CS[worstId,], ncol = numVar))
      newFitness <- calcFitness(FUN, optimType, new)
      prob <- newFitness <= worstFitness
      isItMove <- prob
      if(!all(prob == FALSE)){
        CS[worstId[prob],] <- new[prob,]
      }
      # if best fail to update worst
      if(!all(isItMove == TRUE)){
        randomMatrix <- matrix(runif(length(worstId) * numVar), ncol = numVar)
        new <- matrix(CS[worstId,], ncol = numVar) + randomMatrix * t(gbest - t(CS[worstId,]))
        newFitness <- calcFitness(FUN, optimType, new)
        prob <- newFitness <= worstFitness
        prob <- prob & !isItMove
        if(!all(prob == FALSE)){
          CS[worstId[prob],] <- new[prob,]
          isItMove <- isItMove | prob
        }
      }
      # if gbest fail to update worst
      if(!all(isItMove == TRUE)){
        CS[worstId[!isItMove],] <- generateRandom(length(which(isItMove == F)),numVar, lowerBound, upperBound)
      }
      # reorder
      CSFitness <- calcFitness(FUN, 1, CS)
      CS <- matrix(CS[order(CSFitness),], ncol = numVar)
      memeplexId <- memeplexId[order(CSFitness)]
      CS <- matrix(CS[order(memeplexId),], ncol = numVar)
      memeplexId <- memeplexId[order(memeplexId)]
      gbest <- calcBest(FUN, -1*optimType, rbind(CS, gbest))
    }
    setTxtProgressBar(progressbar, t)
  }
  close(progressbar)
  return(gbest)
}
