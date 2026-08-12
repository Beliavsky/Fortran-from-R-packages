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
#' This is the internal function that implements Black-Hole based Optimization
#' Algorithm. It is used to solve continuous optimization tasks.
#' Users do not need to call it directly,
#' but just use \code{\link{metaOpt}}.
#'
#' This algorithm was proposed by (Hatamlou, 2013).
#' The main inspiration for BHO algorithm originates from black hole
#' that swallow all nearest star. Black hole represent candidate solution
#' with best fitness and other candidate solutions as star, so all star
#' search new best candidate solution while moving towards black-hole.
#' if star reaches better fitness than black hole, exchange its position.
#' star that too close to black hole (pass event horizon) wiil be replace
#' by new random candidate solution.
#'
#' In order to find the optimal solution, the algorithm follow the following steps.
#' \itemize{
#' \item initialize population randomly.
#' \item select best candidate solution as black hole other as stars.
#' \item change each star location to moving toward black hole.
#' \item If a star reaches a location with lower cost than the black hole, exchange their locations.
#' \item If a star crosses the event horizon of the black hole, replace it
#'      with a new star in a random location in the search space.
#' \item If a termination criterion (a maximum number of iterations or a sufficiently good fitness) is met,
#'       exit the loop.
#' }
#'
#' @title Optimization using Black Hole Optimization Algorithm
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
#' ## calculate the optimum solution using black hole optimization
#' resultBHO <- BHO(step, optimType="MIN", numVar, numPopulation=20,
#'                  maxIter=100, rangeVar)
#'
#' ## calculate the optimum value using step function
#' optimum.value <- step(resultBHO)
#'
#' @return \code{Vector [v1, v2, ..., vn]} where \code{n} is number variable
#'         and \code{vn} is value of \code{n-th} variable.
#' @references
#' Hatamlou, A. (2013). Black hole: A new heuristic optimization approach for data clustering.
#' Information Sciences, 222(December), 175–184. https://doi.org/10.1016/j.ins.2012.08.023
#'
#' @export

# Black Hole-based Optimization (BHO)

BHO <- function(FUN, optimType="MIN", numVar, numPopulation=40, maxIter=500, rangeVar){
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
  candidateSolution <- generateRandom(numPopulation, dimension, lowerBound, upperBound)
  bestPos <- engineBHO(FUN, optimType, maxIter, lowerBound, upperBound, candidateSolution)

  return(bestPos)
}

engineBHO <- function(FUN, optimType, maxIter, lowerBound, upperBound, candidateSolution){
  numVar <- ncol(candidateSolution)
  numPopulation <- nrow(candidateSolution)

  progressbar <- txtProgressBar(min = 0, max = maxIter, style = 3)
  for(t in 1:maxIter){
    blackhole <- calcBest(FUN, -1*optimType, candidateSolution)
    # change location each star candidate solution (cs)
    randomMatrix <- matrix(runif(numPopulation * numVar), ncol = numVar)
    candidateSolution <- candidateSolution + randomMatrix * t(blackhole - t(candidateSolution))

    # if a star reaches better fitness exchange it with blackhole
    blackhole <- calcBest(FUN, -1*optimType, candidateSolution)

    # if a star cross event horizon generate new candidate solution for that star
    bhFitness <- calcFitness(FUN, optimType, matrix(blackhole, ncol = numVar, byrow = T))
    fitness <- calcFitness(FUN, optimType, candidateSolution)
    eventHorizon <- bhFitness / sum(fitness)
    if(is.nan(eventHorizon)) eventHorizon <- 0
    isCrossEventHorizon <- abs(fitness - bhFitness) < eventHorizon
    isCrossEventHorizon[order(fitness)[1]] <- F # blackhole exception
    if(!all(isCrossEventHorizon == FALSE)){
      candidateSolution[isCrossEventHorizon, ] <- generateRandom(length(which(isCrossEventHorizon == TRUE)), numVar, lowerBound, upperBound)
    }
    setTxtProgressBar(progressbar, t)
  }
  close(progressbar)
  return(calcBest(FUN, -1*optimType, candidateSolution[, 1:numVar]))
}
