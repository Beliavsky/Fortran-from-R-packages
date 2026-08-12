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
#' This is the internal function that implements Clonal Selection
#' Algorithm. It is used to solve continuous optimization tasks.
#' Users do not need to call it directly,
#' but just use \code{\link{metaOpt}}.
#'
#' This algorithm was proposed by (Castro & Zuben, 2002). The Clonal Selection Algorithm (CLONALG)
#' mimics maturation proses of imumune system. CLONALG consist 5 step initialize, selection, clonal,
#' hypermutation and maturation.
#'
#' In order to find the optimal solution, the algorithm follow the following steps.
#' \itemize{
#' \item initialize population randomly.
#' \item select top selectionSize candidate solutions from population with best fitness.
#' \item clone each selected candidate solutions.
#' \item hypermutation each variable in cloned candidate solutions.
#' \item maturation combine each hypermutated candidate solution with population.
#'       Select top n candidate solution from population as new population.
#' \item If a termination criterion (a maximum number of iterations or a sufficiently good fitness) is met,
#'       exit the loop.
#' }
#'
#' @title Optimization using Clonal Selection Algorithm
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
#' @param selectionSize a positive integer between 0 and numVar
#'        to determine selection size (see details). The default value is \code{as.integer(numPopulation/4)}.
#'
#' @param multipicationFactor a positive numeric between 0 and 1 to determine number of clones. The default value is 0.5.
#'
#' @param hypermutationRate a positive numeric between 0 and 1 to determine probabilty of variable in
#'        clone candidate solutions to be mutated, close to 1 probability is high and vice versa.
#'        The default value is 0.1.
#'
#' @importFrom stats runif
#' @importFrom utils setTxtProgressBar txtProgressBar
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
#' ## calculate the optimum solution clonal selection algorithm
#' resultCLONALG <- CLONALG(quartic, optimType="MIN", numVar, numPopulation=20,
#'                  maxIter=100, rangeVar)
#'
#' ## calculate the optimum value using quartic with noise function
#' optimum.value <- quartic(resultCLONALG)
#'
#' @return \code{Vector [v1, v2, ..., vn]} where \code{n} is number variable
#'         and \code{vn} is value of \code{n-th} variable.
#' @references
#' Castro, L. & Zuben, F. J. V. (2002).
#' Learning and optimization using the clonal selection principle.
#' IEEE Transactions on Evolutionary Computation, Special Issue on
#' Artificial. Immune Systems, 6(3), 239–251. https://doi.org/10.1109/TEVC.2002.1011539
#'
#' @export
# Clonal Selection Algorithm (CLONALG)

CLONALG <- function(FUN, optimType="MIN", numVar, numPopulation=40, maxIter=500, rangeVar,
                    selectionSize=as.integer(numPopulation/4), multipicationFactor=0.5, hypermutationRate=0.1){
  # Validation
  if(numPopulation < 1){
    stop("numPopulation must greater than 0")
  }

  if(maxIter < 0){
    stop("maxIter must greater than or equal to 0")
  }

  if(!is.integer(selectionSize)){
    stop("selectionSize must be integer")
  }else if(selectionSize > numPopulation | selectionSize <= 0){
    stop("selectionSize must less than or equal numPopulation")
  }

  if(multipicationFactor < 0 | multipicationFactor > 1) stop("multipicationFactor must between 0 and 1")
  if(hypermutationRate < 0 | hypermutationRate > 1) stop("hypermutationRate must between 0 and 1")

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
  candidateSolutions <- generateRandom(numPopulation, dimension, lowerBound, upperBound)
  bestPos <- engineCLONALG(FUN, optimType, maxIter, rangeVar, lowerBound, upperBound,
                           selectionSize, multipicationFactor, hypermutationRate, candidateSolutions)
  return(bestPos)
}

engineCLONALG <- function(FUN, optimType, maxIter, rangeVar, lowerBound, upperBound,
                          selectionSize, multipicationFactor, hypermutationRate,
                          candidateSolution){
  numVar <- ncol(candidateSolution)
  numPopulation <- nrow(candidateSolution)
  # evaluate candidate solution
  fitness <- calcFitness(FUN, optimType, candidateSolution)
  candidateSolutions <- data.frame(candidateSolution, fitness)

  progressbar <- txtProgressBar(min = 0, max = maxIter, style = 3)
  for(t in 1:maxIter){
    # Select top "selectionSize" with best fitness from candidateSolutions as topSelections
    candidateSolutions <- candidateSolutions[order(candidateSolutions$fitness), ]
    topSelections <- data.matrix(candidateSolutions[1:selectionSize, 1:numVar])

    # clone topSelections as clone
    # create empty matrix
    clone <- matrix()[-1:-1]
    for(i in 1:selectionSize){
      numClone <- round(multipicationFactor*numPopulation/i, digits = 0)
      clone <- rbind(clone, matrix(data=rep(topSelections[i,], numClone), ncol = numVar, byrow = TRUE))
    }

    # hypermutate clone
    probMatrix <- apply(matrix(NA, ncol = ncol(clone), nrow= nrow(clone), byrow = TRUE), c(1, 2), function(x){
      runif(1)
    })
    rangeVarMatrix <- rbind(lowerBound[col(clone)[probMatrix <= hypermutationRate]],
                            upperBound[col(clone)[probMatrix <= hypermutationRate]])

    clone[probMatrix <= hypermutationRate] <- apply(rangeVarMatrix, c(2), function(x){
      runif(1, min = x[1], max = x[2])
    })

    # maturation step
    candidateSolution <- clone
    fitness <- calcFitness(FUN, optimType, candidateSolution)
    candidateSolutions <- rbind(candidateSolutions, data.frame(candidateSolution, fitness))
    candidateSolutions <- candidateSolutions[order(candidateSolutions$fitness), ]
    candidateSolutions <- candidateSolutions[1:numPopulation, ]

    setTxtProgressBar(progressbar, t)
  }
  close(progressbar)
  return(unname(calcBest(FUN, optimType, as.matrix(candidateSolutions[, 1:numVar]))))
}
