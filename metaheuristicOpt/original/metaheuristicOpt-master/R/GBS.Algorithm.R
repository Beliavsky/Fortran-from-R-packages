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
#' This is the internal function that implements Gravitational Based Search
#' Algorithm. It is used to solve continuous optimization tasks.
#' Users do not need to call it directly,
#' but just use \code{\link{metaOpt}}.
#'
#' This algorithm was proposed by (Rashedi, 2009).
#' GBS use newton law of universal gravitation and second law of motion
#' to optimize. Every candidate solution in population consider having mass and it move
#' using newton law of universal gravitation and second law of motion.
#'
#' In order to find the optimal solution, the algorithm follow the following steps.
#' \itemize{
#' \item initialize population randomly.
#' \item calculate gravitational mass of every candidate solution in population.
#' \item calculate total force of every candidate solution in population using
#'       newton law of universal gravitation.
#' \item calculate acceleration of every candidate solution in population using
#'       newton second law of motion.
#' \item update velocity of every candidate solution in population based on its acceleration.
#' \item move every candidate solution in population based on its velocity.
#' \item If a termination criterion (a maximum number of iterations or a sufficiently good fitness) is met,
#'       exit the loop, else back to calculate gravitational mass.
#' }
#'
#' @title Optimization using Gravitational Based Search Algorithm.
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
#' @param gravitationalConst a numeric to determine gravitational constant while
#'        calculating total force. The default value is \code{max(rangeVar)}.
#'
#' @param kbest a positive numeric between 0 and 1 to determine fraction of population
#'        with best fitness which will affect every candidate solution in population.
#'        The default value is 0.1.
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
#' ## calculate the optimum solution using Gravitational Based Search
#' resultGBS <- GBS(schewefels2.22, optimType="MIN", numVar, numPopulation=20,
#'                  maxIter=100, rangeVar)
#'
#' ## calculate the optimum value using schewefel's problem 2.22 function
#' optimum.value <- schewefels2.22(resultGBS)
#'
#' @return \code{Vector [v1, v2, ..., vn]} where \code{n} is number variable
#'         and \code{vn} is value of \code{n-th} variable.
#' @references
#' Rashedi, E., Nezamabadi-Pour, H., & Saryazdi, S. (2009).
#' GSA: a gravitational search algorithm. Information sciences, 179(13), 2232-2248.
#'
#' @export
# Gravitational Based Search Algorithm (GBS)

GBS <- function(FUN, optimType="MIN", numVar, numPopulation=40, maxIter=500, rangeVar,
                gravitationalConst=max(rangeVar), kbest=0.1){
  # Validation
  if(numPopulation < 1){
    stop("numPopulation must greater than 0")
  }

  if(maxIter < 0){
    stop("maxIter must greater than or equal to 0")
  }

  if(kbest <= 0 | kbest > 1){
    stop("kbest must between 0 and 1")
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

  # generate initial population
  candidateSolution <- generateRandom(numPopulation, dimension, lowerBound, upperBound)
  bestPos <- engineGBS(FUN, optimType, maxIter, lowerBound, upperBound, candidateSolution,
                       gravitationalConst, kbest)
  return(bestPos)
}

engineGBS <- function(FUN, optimType, maxIter, lowerBound, upperBound, candidateSolution,
                      gravitationalConst=100, kbest=0.5){
  numVar <- ncol(candidateSolution)
  numPopulation <- nrow(candidateSolution)
  # give every candidate solution initial velocity
  velocity <- matrix(rep(0, numPopulation * numVar), ncol = numVar)

  gbest <- calcBest(FUN, -1*optimType, candidateSolution)
  progressbar <- txtProgressBar(min = 0, max = maxIter, style = 3)
  for(t in 1:maxIter){
    # update gravitational constant
    eulerConst <- 0.5772156649
    gravitationalConstT <- gravitationalConst / exp(0.01 * t)
    # get best and worst candidate solution in this current t
    CSFitness <- calcFitness(FUN, optimType, candidateSolution)
    best <- candidateSolution[order(CSFitness)[1],]
    worst <- candidateSolution[order(CSFitness)[numPopulation],]
    bestFitness <- calcFitness(FUN, optimType, matrix(best, ncol = numVar))
    worstFitness <- calcFitness(FUN, optimType, matrix(worst, ncol = numVar))

    # calculate gravitaional mass for each candidate solution
    totalGM <- sum((CSFitness - worstFitness) / (bestFitness - worstFitness))
    if(is.nan(totalGM)){
      print("out of bound")
      break
    }
    GM <- (CSFitness - worstFitness) / (bestFitness - worstFitness)
    GM <- GM / totalGM

    # calculate total force
    epsilon <- 8.854e-12
    k <- round(numPopulation * kbest)
    candidateSolution <- as.matrix(candidateSolution[as.numeric(order(CSFitness)), ])
    velocity <- as.matrix(velocity[as.numeric(order(CSFitness)), ])
    totalForce <- c()
    for(i in 1:numPopulation){
      totalI <- rep(0, numVar)
      for(j in 1:k){
        distance <- as.numeric(dist(rbind(candidateSolution[i,],candidateSolution[j,])))
        Force <- runif(1) * gravitationalConstT * GM[j] * (candidateSolution[j,] - candidateSolution[i, ])/distance
        if(all(is.nan(Force))){
          Force <- as.numeric(rep(0, numVar))
        }
        totalI <- totalI + Force
      }
      totalForce <- rbind(totalForce, totalI)
    }

    randomMatrix <- matrix(runif(numPopulation * numVar), ncol = numVar)
    velocity <- randomMatrix * velocity + totalForce
    candidateSolution <- candidateSolution + velocity
    gbest <- calcBest(FUN, -1*optimType, rbind(candidateSolution, gbest))

    setTxtProgressBar(progressbar, t)
  }
  close(progressbar)
  #checkBound
  gbest <- checkBound(gbest, lowerBound, upperBound)
  return(gbest)
}
