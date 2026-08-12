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
#' This is the internal function that implements Differential Evolution
#' Algorithm. It is used to solve continuous optimization tasks.
#' Users do not need to call it directly,
#' but just use \code{\link{metaOpt}}.
#'
#' This Differential Evolution algorithm based on jurnal by (Das & Suganthan, 2011).
#' Differential Evolution algorithm use genetic operator for optimization such as
#' mutation, crossover and selection.
#'
#' In order to find the optimal solution, the algorithm follow the following steps.
#' \itemize{
#' \item initialize population randomly.
#' \item create some mutation vectors as new candidate solutions (mutation operator).
#' \item perform crossover operator.
#' \item perform selection operator.
#' \item If a termination criterion (a maximum number of iterations or a sufficiently good fitness) is met,
#'       exit the loop, else back to create some mutation vector.
#' }
#'
#' @title Optimization using Differential Evolution Algorithm
#'
#' @param FUN an objective function or cost function
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
#' @param scalingVector a positive numeric between 0 and 1 to determine scalingVector
#'        for mutation operator. The default value is 0.8.
#'
#' @param crossOverRate a positive numeric between 0 and 1 to determine crossOver probability.
#'        The default value is 0.5.
#'
#' @param strategy characters to determine mutation method. They are six methods to choose:
#' \itemize{
#' \item "classical".
#' \item "best 1"
#' \item "target to best"
#' \item "best 2"
#' \item "rand 2"
#' \item "rand 2 dir"
#' }
#' details of the mutation methods are on the references. The default value is "best 1".
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
#' ## calculate the optimum solution using differential evolution
#' resultDE <- DE(step, optimType="MIN", numVar, numPopulation=20,
#'                  maxIter=100, rangeVar)
#'
#' ## calculate the optimum value using step function
#' optimum.value <- step(resultDE)
#'
#' @return \code{Vector [v1, v2, ..., vn]} where \code{n} is number variable
#'         and \code{vn} is value of \code{n-th} variable.
#' @references
#' Das, S., & Suganthan, P. N. (2011). Differential evolution: A survey of the state-of-the-art.
#' IEEE transactions on evolutionary computation, 15(1), 4-31.
#'
#' @export
# Differential Evolution (DE)

DE <- function(FUN, optimType="MIN", numVar, numPopulation=40, maxIter=500, rangeVar,
               scalingVector=0.8, crossOverRate=0.5, strategy="best 1"){
  # Validation
  if(numPopulation < 1){
    stop("numPopulation must greater than 0")
  }

  if(maxIter < 0){
    stop("maxIter must greater than or equal to 0")
  }
  # check parameter scalingVector
  if(scalingVector < 0 || scalingVector > 1){
    stop("parameter scalingVector must between 0 and 1")
  }

  # check parameter crossOverRate
  if(crossOverRate < 0 || crossOverRate > 1){
    stop("parameter crossOverRate must between 0 and 1")
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
  candidateSolution <- generateRandomDE(numPopulation, dimension, lowerBound, upperBound)
  bestPos <- engineDE(FUN, optimType, maxIter, lowerBound, upperBound, candidateSolution,
                      scalingVector, crossOverRate, strategy)

  return(bestPos)
}

engineDE <- function(FUN, optimType, maxIter, lowerBound, upperBound, candidateSolutions,
                     scalingVector, crossOverRate, strategy){
  numVar <- ncol(candidateSolutions)
  numPopulation <- nrow(candidateSolutions)

  index <- 1:numPopulation
  progressbar <- txtProgressBar(min = 0, max = maxIter, style = 3)
  for(t in 1:maxIter){
    # create mutant vector
    indexMutant1 <- sample(index, numPopulation)
    isExclusive <- FALSE
    while(isExclusive == FALSE){
      indexMutant2 <- sample(index, numPopulation)
      if(all((indexMutant2 == indexMutant1) == FALSE)){
        isExclusive <- TRUE
      }
    }
    mutant1 <- candidateSolutions[indexMutant1,]
    mutant2 <- candidateSolutions[indexMutant2,]
    # create mutant 3
    if(strategy == "best 2" | strategy == "rand 2 dir" | strategy == "rand 2"){
      isExclusive <- FALSE
      while(isExclusive == FALSE){
        indexMutant3 <- sample(index, numPopulation)
        if((all((indexMutant3 == indexMutant1) == FALSE))&(all((indexMutant3 == indexMutant2) == FALSE))){
          mutant3 <- candidateSolutions[indexMutant3,]
          isExclusive <- TRUE
        }
      }
    }
    # create mutant 4
    if(strategy == "best 2" | strategy == "rand 2"){
      isExclusive <- FALSE
      while(isExclusive == FALSE){
        indexMutant4 <- sample(index, numPopulation)
        if((all((indexMutant4 == indexMutant1) == FALSE))&(all((indexMutant4 == indexMutant2) == FALSE))&(all((indexMutant4 == indexMutant3) == FALSE))){
          mutant4 <- candidateSolutions[indexMutant4,]
          isExclusive <- TRUE
        }
      }
    }
    # create mutant 5
    if(strategy == "rand 2"){
      isExclusive <- FALSE
      while(isExclusive == FALSE){
        indexMutant5 <- sample(index, numPopulation)
        if((all((indexMutant5 == indexMutant1) == FALSE))&(all((indexMutant5 == indexMutant2) == FALSE))&(all((indexMutant5 == indexMutant3) == FALSE))&(all((indexMutant5 == indexMutant4) == FALSE))){
          mutant5 <- candidateSolutions[indexMutant5,]
          isExclusive <- TRUE
        }
      }
    }

    # apply based on strategy
    if(strategy == "clasical"){
      mutant <- candidateSolutions + scalingVector * (mutant1 - mutant2)
    }else if(strategy == "best 1"){
      best <- calcBest(FUN, -1*optimType, candidateSolutions)
      mutant <- t(best + t(scalingVector * (mutant1 - mutant2)))
    }else if(strategy == "target to best"){
      best <- calcBest(FUN, -1*optimType, candidateSolutions)
      mutant <- candidateSolutions + scalingVector * t(best - t(candidateSolutions)) +
        scalingVector * (mutant1 - mutant2)
    }else if(strategy == "best 2"){
      best <- calcBest(FUN, -1*optimType, candidateSolutions)
      mutant <- t(best + t(scalingVector * (mutant1 - mutant2) + scalingVector * (mutant3 - mutant4)))
    }else if(strategy == "rand 2"){
      mutant <- mutant1 + scalingVector * (mutant2 - mutant3) + scalingVector * (mutant4 - mutant5)
    }else if(strategy == "rand 2 dir"){
      mutant <- mutant1 + scalingVector/2*(mutant1 - mutant2 - mutant3)
    }else{
      stop(paste("there is no strategy name", strategy))
    }

    # crossover(CO) step
    COProb <- matrix(runif(numPopulation * numVar), ncol = numVar)
    COLogical <- matrix(rep(F, numPopulation * numVar), ncol = numVar)
    COLogical[COProb <= crossOverRate] <- T

    COLogical <- !COLogical
    mutant[COLogical] <- candidateSolutions[COLogical]

    # Selection step
    mutantFitness <- calcFitness(FUN, optimType, mutant)
    CSFitness <- calcFitness(FUN, optimType, candidateSolutions)
    isMutantBetterThanCS <- mutantFitness < CSFitness
    candidateSolutions[isMutantBetterThanCS,] <- mutant[isMutantBetterThanCS,]

    # check bound
    for(i in 1:nrow(candidateSolutions)){
      candidateSolutions[i,] <- checkBound(candidateSolutions[i,], lowerBound, upperBound)
    }

    setTxtProgressBar(progressbar, t)
  }
  close(progressbar)
  return(calcBest(FUN, -1*optimType, candidateSolutions))
}

generateRandomDE <- function(numPopulation, dimension, lowerBound, upperBound){
  matrixLowerBound <- matrix(rep(lowerBound, numPopulation), ncol = dimension, byrow = TRUE)
  matrixUpperBound <- matrix(rep(upperBound, numPopulation), ncol = dimension, byrow = TRUE)
  matrixRandom <- apply(matrix(ncol = dimension, nrow = numPopulation), c(1, 2), function(x){
    runif(1)
  })
  result <- matrixLowerBound + matrixRandom * (matrixUpperBound - matrixLowerBound)
  return(result)
}
