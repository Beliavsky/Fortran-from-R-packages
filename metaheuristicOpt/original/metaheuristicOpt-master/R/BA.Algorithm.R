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
#' This is the internal function that implements Bat
#' Algorithm. It is used to solve continuous optimization tasks.
#' Users do not need to call it directly,
#' but just use \code{\link{metaOpt}}.
#'
#' This algorithm was proposed by (Yang, 2011). It was inspired
#' by echolocation of bats. Candidate solutions in bat algorithm
#' are represented by bat. They have flying speed, pulse rate,
#' loudness and pulse frequency and they move based on them.
#'
#' In order to find the optimal solution, the algorithm follow the following steps.
#' \itemize{
#' \item initialize population randomly.
#' \item move every candidate solutions based on velocity and pulse frequnecy.
#' \item move some candidate solutions near globak best randomly.
#' \item If a candidate solution have better fitness than global best replace it with new
#'       random candidate solution then increase its pulse rate and decrease its loudness.
#' \item If a termination criterion (a maximum number of iterations or a sufficiently good fitness) is met,
#'       exit the loop, else back to move every candidate solutions.
#' }
#'
#' @title Optimization using Bat Algorithm
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
#' @param maxFrequency a numeric to determine maximum frequency. The default value is 0.1.
#'
#' @param minFrequency a numeric to determine minimum frequency. The default value is -0.1.
#'
#' @param gama a numeric greater than equal to 1. It use to increase pulse rate. The default value is 1.
#'
#' @param alphaBA a numeric between 0 and 1. It use to decrease loudness. The default value is 0.1.
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
#' ## calculate the optimum solution using bat algorithm
#' resultBA <- BA(schewefels1.2, optimType="MIN", numVar,
#' numPopulation=20, maxIter=100, rangeVar)
#'
#' ## calculate the optimum value using schewefel's problem 1.2 function
#' optimum.value <- schewefels1.2(resultBA)
#'
#' @return \code{Vector [v1, v2, ..., vn]} where \code{n} is number variable
#'         and \code{vn} is value of \code{n-th} variable.
#' @references
#' Yang, X. S., (2011), Bat Algorithm for Multiobjective Optimization,
#' Int. J. Bio-Inspired Computation, Vol. 3, No. 5, pp.267-274.
#'
#' @export
# Bat Algorithm (BA)

BA <- function(FUN, optimType="MIN", numVar, numPopulation=40, maxIter=500, rangeVar,
               maxFrequency=0.1, minFrequency=-0.1, gama=1, alphaBA=0.1){
  # Validation
  if(numPopulation < 1){
    stop("numPopulation must greater than 0")
  }

  if(maxIter < 0){
    stop("maxIter must greater than or equal to 0")
  }

  if(alphaBA < 0 | alphaBA > 1){
    stop("alpha must between 0 and 1")
  }

  if(gama < 1){
    stop("gama must greater or equal to 1")
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
  bestPos <- engineBA(FUN, optimType, maxIter, lowerBound, upperBound, candidateSolution,
                      maxFrequency, minFrequency, gama, alphaBA)

  return(bestPos)
}

engineBA <- function(FUN, optimType, maxIter, lowerBound, upperBound, candidateSolution,
                     maxFrequency, minFrequency, gama, alpha){
  numVar <- ncol(candidateSolution)
  numPopulation <- nrow(candidateSolution)

  pf <- matrix(runif(numPopulation * numVar), ncol = numVar)
  velocity <- matrix(rep(0, numPopulation * numVar), ncol = numVar)
  A <- runif(numPopulation)
  pr <- runif(numPopulation)

  best <- c()
  progressbar <- txtProgressBar(min = 0, max = maxIter, style = 3)
  for(t in 1:maxIter){
    best <- calcBest(FUN, -1*optimType, rbind(candidateSolution, best))
    bestFitness <- calcFitness(FUN, optimType, matrix(best, ncol = numVar))

    # update pf, velocity and candidate solution
    randomMatrix <- matrix(runif(numPopulation * numVar), ncol = numVar)
    pf <- minFrequency + randomMatrix *(maxFrequency - minFrequency)
    velocity <- velocity + t(t(candidateSolution) - best) * pf
    candidateSolution <- candidateSolution + velocity
    candidateSolution <- checkBound(candidateSolution, lowerBound, upperBound)

    # pr comparison
    best <- calcBest(FUN, -1*optimType, rbind(candidateSolution, best))
    prob <- runif(numPopulation) > pr
    if(!all(prob == FALSE)){
      candidateSolution[prob, ] <- outer(A[prob], best, FUN = "*")
    }

    # loudness comparison
    flyingRandomly <- generateRandom(numPopulation, numVar, lowerBound, upperBound)
    frFitness <- calcFitness(FUN, optimType, flyingRandomly)
    best <- calcBest(FUN, -1*optimType, rbind(candidateSolution, best))
    bestFitness <- calcFitness(FUN, optimType, matrix(best, ncol = numVar))
    prob <- frFitness < bestFitness & runif(numPopulation) < A
    if(!all(prob == FALSE)){
      candidateSolution[prob, ] <- flyingRandomly[prob,]
      pr[prob] <- pr[prob]*(1 - exp(-1*gama))
      A[prob] <- alpha*A[prob]
    }

    setTxtProgressBar(progressbar, t)
  }
  close(progressbar)
  result <- calcBest(FUN, -1*optimType, rbind(candidateSolution, best))
  result <- checkBound(result, lowerBound, upperBound)
  return(result)
}
