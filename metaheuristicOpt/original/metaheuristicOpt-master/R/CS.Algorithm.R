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
#' This is the internal function that implements cuckoo search
#' Algorithm. It is used to solve continuous optimization tasks.
#' Users do not need to call it directly,
#' but just use \code{\link{metaOpt}}.
#'
#' This algorithm was proposed by (Yang & Deb, 2009). This
#' algorithhm was inspired by behaviour of cuckoo birds which
#' place its egg on other bird nest. While cuckoo birds putting
#' the eggs in the nests of other birds they are two possible
#' outcome. First the owner of the nest will stay on the nest.
#' Second the owner of the nest will abandon the nest.
#'
#' In order to find the optimal solution, the algorithm follow the following steps.
#' \itemize{
#' \item initialize population randomly.
#' \item create a mutant vector.
#' \item select a candidate solution in population randomly then compare it with mutant vector.
#'       if mutant vector have better fitness than candidate solution replace candidate solution with
#'       mutant vector.
#' \item replace fraction of population with worst fitness with new random candidate solutions.
#' \item If a termination criterion (a maximum number of iterations or a sufficiently good fitness) is met,
#'       exit the loop, else back to create a mutant vector.
#' }
#'
#' @title Optimization using Cuckoo Search algorithm
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
#' @param abandonedFraction a positive numeric between 0 and 1 to determine fraction
#'        of population to be replaced. The default value is 0.5.
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
#' ## calculate the optimum solution cuckoo search
#' resultCS <- CS(sphere, optimType="MIN", numVar, numPopulation=20,
#'                  maxIter=100, rangeVar)
#'
#' ## calculate the optimum value using sphere function
#' optimum.value <- sphere(resultCS)
#'
#' @return \code{Vector [v1, v2, ..., vn]} where \code{n} is number variable
#'         and \code{vn} is value of \code{n-th} variable.
#' @references
#' Yang, X. S., & Deb, S. (2009, December). Cuckoo search via Lévy flights. In 2009 World Congress on
#' Nature & Biologically Inspired Computing (NaBIC) (pp. 210-214). IEEE.
#'
#' @export
# Cuckoo Search (CS)

CS <- function(FUN, optimType="MIN", numVar, numPopulation=40, maxIter=500, rangeVar, abandonedFraction=0.5){
  # Validation
  if(numPopulation < 1){
    stop("numPopulation must greater than 0")
  }

  if(maxIter < 0){
    stop("maxIter must greater than or equal to 0")
  }

  if(abandonedFraction < 0 | abandonedFraction > 1){
    stop("abandonedFraction must between 0 and 1")
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

  #generate candidate solutions
  candidateSolution <- generateRandom(numPopulation, dimension, lowerBound, upperBound)
  bestPos <- engineCS(FUN, optimType, maxIter, lowerBound, upperBound, candidateSolution, abandonedFraction)

  return(bestPos)
}

engineCS <- function(FUN, optimType, maxIter, lowerBound, upperBound, candidateSolution, abandonedFraction){
  numPopulation <- nrow(candidateSolution)
  numVar <- ncol(candidateSolution)
  best <- c()

  progressbar <- txtProgressBar(min = 0, max = maxIter, style = 3)
  for(t in 1:maxIter){
    best <- calcBest(FUN, -1*optimType, rbind(candidateSolution, best))
    # create a cuckoo egg
    choosen <- sample.int(numPopulation, 1)
    choosen <- candidateSolution[choosen,]
    cuckooegg <- choosen + levyFlight(choosen, best)

    # compare with candidate solution
    otherChoosen <- sample.int(numPopulation, 1)
    otherFitness <- calcFitness(FUN, optimType, matrix(candidateSolution[otherChoosen, ], ncol = numVar))
    cuckooFitness <- calcFitness(FUN, optimType, matrix(cuckooegg, ncol = numVar))
    if(cuckooFitness < otherFitness){
      candidateSolution[otherChoosen, ] <- cuckooegg
    }

    # remove fraction
    CSFitness <- calcFitness(FUN, optimType, candidateSolution)
    candidateSolution <- as.matrix(candidateSolution[order(CSFitness), ])
    numCSToRemove <- round(numPopulation * abandonedFraction)
    removeIndex <- as.numeric((numPopulation - numCSToRemove + 1):numPopulation)
    candidateSolution[removeIndex,] <- generateRandom(numCSToRemove, numVar, lowerBound, upperBound)

    # check bound
    for(i in 1:nrow(candidateSolution)){
      candidateSolution[i,] <- checkBound(candidateSolution[i,], lowerBound, upperBound)
    }
    setTxtProgressBar(progressbar, t)
  }
  close(progressbar)
  best <- calcBest(FUN, -1*optimType, rbind(candidateSolution, best))
  return(unname(best))
}

levyFlight <- function(CS, best){
  n <- length(CS)
  beta <- 3/2
  sigma <- (gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta)
  u <- rnorm(n) * sigma
  v <- rnorm(n)
  step <- u/abs(v)^(1/beta)
  stepsize <- 0.01 * step * (CS - best)
  stepsize <- stepsize * rnorm(n)
  return(stepsize)
}
