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
#' This is the internal function that implements Cat Swarm Optimization
#' Algorithm. It is used to solve continuous optimization tasks.
#' Users do not need to call it directly,
#' but just use \code{\link{metaOpt}}.
#'
#' This algorithm was proposed by (Chu, Tsai & Pan, 2006).
#' This algorithm was inspired by behaviours of felyne. Behaviours of
#' felyne can be devided into two seeking mode (when flyne rest)
#' and tracing mode (when felyne chase its prey). candidate solutions divided
#' into seeking and tracing mode. candidate solution in seeking mode move using
#' local search while candidate solution in tracing mode move using genetic operator.
#'
#' In order to find the optimal solution, the algorithm follow the following steps.
#' \itemize{
#' \item initialize population randomly.
#' \item flaging (tracing or seeking) every candidate solution in population based on mixtureRatio randomly.
#' \item candidate solutions in seeking mode move using local search
#' \item candidate solutions in tracing mode move using genetic operator
#' \item If a termination criterion (a maximum number of iterations or a sufficiently good fitness) is met,
#'       exit the loop, else back to flaging candidate solutions.
#' }
#'
#' @title Optimization using Cat Swarm Optimization Algorithm
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
#' @param mixtureRatio a positive numeric between 0 and 1 to determine flaging proportion.
#'        higher mixtureRatio increase number of candidate solutions in seeking mode
#'        and vice versa. The default value is 0.5.
#'
#' @param tracingConstant a positive numeric between 0 and 1 to determine tracingConstant. The default value is 0.1.
#'
#' @param maximumVelocity a positive numeric to determine maximumVelocity while candidate solutions
#'        in tracing mode performing local search. The default value is 1.
#'
#' @param smp a positive integer to determine number of duplication in genetic operator. The default value is \code{as.integer(20)}.
#'
#' @param srd a positive numeric between 0 and 100 to determine mutation length in genetic operator. The default value is 20.
#'
#' @param cdc a positive integer between 0 and numVar to determine number of variabel in
#'        candidate solutions in seeking mode to be mutated during mutation step in
#'        genetic operator. The default value is \code{as.integer(numVar)}.
#'
#' @param spc a logical. if spc is TRUE smp = smp else smp = smp - 1.  The default value is TRUE.
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
#' resultCSO <- CSO(schewefels2.22, optimType="MIN", numVar, numPopulation=20,
#'                  maxIter=100, rangeVar)
#'
#' ## calculate the optimum value using schewefel's problem 2.22 function
#' optimum.value <- schewefels2.22(resultCSO)
#'
#' @return \code{Vector [v1, v2, ..., vn]} where \code{n} is number variable
#'         and \code{vn} is value of \code{n-th} variable.
#' @references
#' Chu, S. C., Tsai, P. W., & Pan, J. S. (2006, August). Cat swarm optimization.
#' In Pacific Rim international conference on artificial intelligence (pp. 854-858).
#' Springer, Berlin, Heidelberg.
#'
#' @export
# Cat Swarm Optimization (CSO)

CSO <- function(FUN, optimType="MIN", numVar, numPopulation=40, maxIter=500, rangeVar,
                mixtureRatio=0.5, tracingConstant=0.1, maximumVelocity=1, smp=as.integer(20), srd=20, cdc=as.integer(numVar), spc=TRUE){
  # Validation
  if(numPopulation < 1){
    stop("numPopulation must greater than 0")
  }

  if(maxIter < 0){
    stop("maxIter must greater than or equal to 0")
  }

  if(mixtureRatio < 0 | mixtureRatio > 1){
    stop("mixtureRatio must between 0 and 1")
  }

  if(tracingConstant < 0 | tracingConstant > 1){
    stop("mixtureRatio must between 0 and 1")
  }

  if(!is.integer(smp)){
    stop("smp must be integer (as.integer())")
  }

  if(cdc > numVar){
    stop("cdc must less than or equal to numVar")
  }else if(!is.integer(cdc)){
    stop("smp must be integer (as.integer())")
  }

  if(!is.logical(spc)){
    stop("spc must be logical")
  }

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
  bestPos <- engineCSO(FUN, optimType, maxIter, lowerBound, upperBound, candidateSolution,
                       mixtureRatio, tracingConstant, maximumVelocity, smp, srd, cdc, spc)
  return(bestPos)
}

engineCSO <- function(FUN, optimType, maxIter, lowerBound, upperBound, candidateSolution,
                      mixtureRatio, tracingConstant, maximumVelocity, smp, srd, cdc, spc){
  numVar <- ncol(candidateSolution)
  numPopulation <- nrow(candidateSolution)
  # generate candidate solutions
  fitness <- calcFitness(FUN, optimType, candidateSolution)
  velocity <- apply(matrix(rep(NA, numPopulation*numVar), ncol = numVar), c(1, 2), function(x){
    runif(1, min = 0, max = maximumVelocity)
  })
  candidateSolutions <- data.frame(candidateSolution, velocity, fitness)
  bestCandidate <- candidateSolutions[order(candidateSolutions$fitness)[1],]

  progressbar <- txtProgressBar(min = 0, max = maxIter, style = 3)
  for(t in 1:maxIter){
    # Give each candidate solutions flag
    candidateSolutions$flag <- flagingCSO(mixtureRatio, numPopulation)
    # Determine/update best and worst candidate solution
    bestCandidateinThisIteration <- candidateSolutions[order(candidateSolutions$fitness)[1],]
    if(bestCandidate$fitness > bestCandidateinThisIteration$fitness) bestCandidate <- bestCandidateinThisIteration
    worstCandidate <- candidateSolutions[order(candidateSolutions$fitness)[numPopulation],]
    indexVariable <- 1:numVar # index column variable in data frame
    indexVelocity <- (numVar+1):(numVar+numVar)# index column velocity in data frame

    # if flag == "tracing" ----
    tracing <- candidateSolutions[candidateSolutions$flag == "tracing",]
    tracingVariable <- as.matrix(tracing[,indexVariable])
    tracingVelocity <- as.matrix(tracing[,indexVelocity])
    # Update velocity using
    randomMatrix <- apply(tracingVelocity, c(1, 2), function(x){
      runif(1)
    })
    bestCandidateVariable <- as.numeric(bestCandidate[,1:numVar])
    bestCandidateVariable <- matrix(rep(bestCandidateVariable, nrow(tracing)), ncol = numVar, byrow = TRUE)
    tracingVelocity <- tracingVelocity + randomMatrix * tracingConstant * (bestCandidateVariable - tracingVariable)
    # check velocity
    tracingVelocity[tracingVelocity > maximumVelocity] <- maximumVelocity
    # update position
    tracingVariable <- tracingVariable + tracingVelocity
    # update candidate solution
    candidateSolutions[candidateSolutions$flag == "tracing", indexVariable] <- tracingVariable
    candidateSolutions[candidateSolutions$flag == "tracing", indexVelocity] <- tracingVelocity


    # if flag == "seeking" ----
    seeking <- candidateSolutions[candidateSolutions$flag == "seeking",]
    seekingVariable <- seeking[,indexVariable]
    if(numVar == 1){
      x <- seekingVariable
      copyId <- 1:length(seekingVariable)
      seekingVariable <- data.frame(x, copyId)
    }else{
      seekingVariable$copyId <- 1:nrow(seekingVariable)
    }

    # make smp copies
    copies <- data.frame()
    if(spc == TRUE){
      for(i in 1:smp){
        copies <- rbind(copies, seekingVariable)
      }
    }else{
      for(i in 1:(smp-1)){
        copies <- rbind(copies, seekingVariable)
      }
    }

    # modified copies
    if(cdc != 0){
      modified <- apply(as.matrix(copies[,indexVariable]), c(1), function(x){
        pickedVariables <- sample(1:numVar, cdc)
        posOrNeg <- sapply(1:cdc, function(y){
          sample(c(1, -1), 1)
        })
        x[pickedVariables] <- x[pickedVariables]*posOrNeg*srd/100
        return(x)
      })
      if(numVar == 1) copies[,indexVariable] <- modified else copies[,indexVariable] <- t(modified)
    }

    # calculate probabilty of all candidate solution (flag == "seeking")
    copies <- rbind(seekingVariable, copies)
    copies$probability <- probabilityCSO(as.matrix(copies[,indexVariable]), bestCandidate, worstCandidate, FUN, optimType)

    # chose one candidate solution for each copy based on probability
    for(i in 1:nrow(seekingVariable)){
      numCopies <- nrow(copies[copies$copyId == i,])
      probCopies <- copies[copies$copyId == i, "probability"]
      if(all(probCopies == 0)){
        choosenCopy <- sample(1:numCopies, 1, replace = TRUE)
      }else{
        choosenCopy <- sample(1:numCopies, 1, prob = probCopies, replace = TRUE)
      }
      choosenCopy <- copies[copies$copyId == i, ][choosenCopy, indexVariable]
      seekingVariable[seekingVariable$copyId == i, indexVariable] <- choosenCopy
    }

    # update candidate solution
    candidateSolutions[candidateSolutions$flag == "seeking", indexVariable] <- seekingVariable[,indexVariable]

    # update candidate solution fitness
    candidateSolutions$fitness <- calcFitness(FUN, optimType, as.matrix(candidateSolutions[,indexVariable]))

    # check Bound
    tempCS <- as.matrix(candidateSolutions[, indexVariable])
    for(i in 1:nrow(tempCS)){
      tempCS[i,] <- checkBound(tempCS[i,], lowerBound, upperBound)
    }
    candidateSolutions[,indexVariable] <- tempCS
    candidateSolutions$fitness <- calcFitness(FUN, optimType, as.matrix(candidateSolutions[,indexVariable]))

    setTxtProgressBar(progressbar, t)
  }
  close(progressbar)
  return(as.numeric(bestCandidate[,indexVariable]))
}

flagingCSO <- function(mixtureRatio, numPopulation){
  numSeeking <- mixtureRatio * numPopulation
  numTracing <- (1 - mixtureRatio) * numPopulation
  seeking <- rep("seeking", ceiling(numSeeking))
  tracing <- rep("tracing", ceiling(numTracing))
  result <- sample(c(seeking, tracing), replace = FALSE)
  if(length(result) != numPopulation){
    result <- result[1:numPopulation]
  }
  return(result)
}

probabilityCSO <- function(input, best, worst, FUN, optimType){
  inputFitness <- calcFitness(FUN, optimType, input)
  bestFitness <- best$fitness
  worstFitness <- worst$fitness
  result <- abs(inputFitness - bestFitness)/(worstFitness - bestFitness)
  return(result)
}

