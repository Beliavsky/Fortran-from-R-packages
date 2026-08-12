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
#' This is the internal function that implements Krill-Herd
#' Algorithm. It is used to solve continuous optimization tasks.
#' Users do not need to call it directly,
#' but just use \code{\link{metaOpt}}.
#'
#' This algorithm was proposed by (Gandomi & Alavi, 2012).
#' It was inspired by behaviours of swarm of krill. Every
#' krill move based on motion induced (such as obstacle, predators),
#' foraging speed (food source) and physical difussion (swarm density).
#' In KH algorithm candidate solution represented by krill. KH algorithm
#' also use genetic operator mutation and crossover.
#'
#' In order to find the optimal solution, the algorithm follow the following steps.
#' \itemize{
#' \item initialize population randomly.
#' \item calculate total motion based on motion induced, foraging speed and physical
#'       difussion for each candidate solutions and move it based on total motion.
#' \item perform genetic operator crossover and mutation
#' \item If a termination criterion (a maximum number of iterations or a sufficiently good fitness) is met,
#'       exit the loop, else back to calculate total motion.
#' }
#'
#' @title Optimization using Krill-Herd Algorithm
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
#' @param maxMotionInduced a positive numeric between 0 and 1 to determine
#'        maximum motion induced. The default value is 0.01.
#'
#' @param inertiaWeightOfMotionInduced a positive numeric between 0 and 1 to determine
#'        how much motion induced affect krill (candidate solution) movement. the
#'        greater the value the greater the affect of motion induced on krill movement.
#'        The default value is 0.01.
#'
#' @param epsilon a positive numeric between 0 and 1 to determine epsilon constant. The default value is 1e-05.
#'
#' @param foragingSpeed a positive numeric between 0 and 1 to determine foraging speed. The default value is 0.02
#'
#' @param inertiaWeightOfForagingSpeed a positive numeric between 0 and 1 to determine
#'        how much foraging speed affect krill (candidate solution) movement. the
#'        greater the value the greater the affect of foraging speed on krill movement.
#'        The default value is 0.01.
#'
#' @param maxDifussionSpeed a positive numeric between 0 and 1 to determine maximum
#'        difussion speed. The default value is 0.01.
#'
#' @param constantSpace a numeric between 0 and 1 to determine how much range affect
#'        krill movement. The default value is 1.
#'
#' @param mu a numeric between 0 and 1 to determine constant number for mutation operator.
#'        The default value is 0.1.
#'
#' @importFrom graphics plot
#' @importFrom stats runif dist rnorm
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
#' ## calculate the optimum solution
#' resultKH <- KH(sphere, optimType="MIN", numVar, numPopulation=20,
#'                  maxIter=100, rangeVar)
#'
#' ## calculate the optimum value using sphere function
#' optimum.value <- sphere(resultKH)
#'
#' @return \code{Vector [v1, v2, ..., vn]} where \code{n} is number variable
#'         and \code{vn} is value of \code{n-th} variable.
#'
#' @references
#' Gandomi, A. H., & Alavi, A. H. (2012). Krill herd: a new bio-inspired optimization algorithm.
#' Communications in nonlinear science and numerical simulation, 17(12), 4831-4845.
#'
#' @export

# Krill-Heard Algorithm(KH)

KH <- function(FUN, optimType="MIN", numVar, numPopulation=40, maxIter=500, rangeVar,
               maxMotionInduced=0.01, inertiaWeightOfMotionInduced=0.01, epsilon=1e-05, foragingSpeed=0.02,
               inertiaWeightOfForagingSpeed=0.01, maxDifussionSpeed=0.01, constantSpace=1, mu=0.1){
  # Validation
  if(numPopulation < 1){
    stop("numPopulation must greater than 0")
  }

  if(maxIter < 0){
    stop("maxIter must greater than or equal to 0")
  }

  if(inertiaWeightOfMotionInduced < 0 | inertiaWeightOfMotionInduced > 1){
    stop("inertiaWeightOfMotionInduced must between 0 and 1")
  }

  if(inertiaWeightOfForagingSpeed < 0 | inertiaWeightOfForagingSpeed > 1){
    stop("inertiaWeightOfForagingSpeed must between 0 and 1")
  }

  if(maxMotionInduced < 0 | maxMotionInduced > 1){
    stop("maxMotionInduced must between 0 and 1")
  }

  if(foragingSpeed < 0 | foragingSpeed > 1){
    stop("foragingSpeed must between 0 and 1")
  }

  if(maxDifussionSpeed < 0 | maxDifussionSpeed > 1){
    stop("maxDifussionSpeed must between 0 and 1")
  }

  if(constantSpace < 0 | constantSpace > 2){
    stop("constantSpace must between 0 and 2")
  }

  if(mu < 0 | mu > 1){
    stop("mu must between 0 and 1")
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
  candidateSolution <- generateRandom(numPopulation, dimension, lowerBound, upperBound)
  bestPos <- engineKH(FUN, optimType, maxIter, lowerBound, upperBound, candidateSolution,
                      maxMotionInduced, inertiaWeightOfMotionInduced, epsilon, foragingSpeed,
                      inertiaWeightOfForagingSpeed, maxDifussionSpeed, constantSpace, mu)
  return(bestPos)
}

engineKH <- function(FUN, optimType, maxIter, lowerBound, upperBound, candidateSolution,
                     maxMotionInduced, inertiaWeightOfMotionInduced, epsilon, foragingSpeed,
                     inertiaWeightOfForagingSpeed, maxDifussionSpeed, constantSpace, mu){
  numVar <- ncol(candidateSolution)
  numPopulation <- nrow(candidateSolution)

  N <- matrix(rep(0, numPopulation * numVar), ncol = numVar)
  f <- matrix(rep(0, numPopulation * numVar), ncol = numVar)

  gbest <- calcBest(FUN, -1*optimType, candidateSolution)

  progressbar <- txtProgressBar(min = 0, max = maxIter, style = 3)
  for(t in 1:maxIter){
    CSFitness <- calcFitness(FUN, optimType, candidateSolution)
    best <- candidateSolution[order(CSFitness)[1], ]
    worst <- candidateSolution[order(CSFitness)[numPopulation], ]
    bestFitness <- calcFitness(FUN, optimType, matrix(best, ncol = numVar))
    worstFitness <- calcFitness(FUN, optimType, matrix(worst, ncol = numVar))
    gbest <- calcBest(FUN, -1*optimType, rbind(candidateSolution, gbest))


    # motion iduced
    sensingDistance <- 1/5*ncol(candidateSolution)*colSums(as.matrix(dist(candidateSolution)))
    isIncludeSD <- as.matrix(dist(candidateSolution, diag = T, upper = T)) < sensingDistance

    alpha <- c()
    for(index in 1:numPopulation){
      X <- apply(as.matrix(candidateSolution[isIncludeSD[index,],]), c(1), function(x, y){
        xijKH(y, x, epsilon)
      }, y=candidateSolution[index,])

      K <- sapply(CSFitness[isIncludeSD[index,]], function(x, y){
        kijKH(y,x, bestFitness, worstFitness)
      }, y=CSFitness[index])

      if(numVar == 1){
        alphaLocal <- sum(X * K)
      }else{
        alphaLocal <- colSums(t(X) * K)
      }

      Cbest <- 2*(runif(1)+t/maxIter)
      X <- xijKH(candidateSolution[index,], best, epsilon)
      K <- kijKH(CSFitness[index], bestFitness, bestFitness, worstFitness)
      alphaTarget <- Cbest*K*X
      alpha <- rbind(alpha, alphaLocal + alphaTarget)
    }

    N <- maxMotionInduced * alpha + inertiaWeightOfMotionInduced * N

    # foraging motion
    if(numVar == 1){
      Xfood  <- sum(candidateSolution * 1 / CSFitness) / sum(1/CSFitness)
      if(is.nan(Xfood)){
        Xfood <- 0
      }
    }else{
      Xfood  <- colSums(candidateSolution * 1 / CSFitness) / sum(1/CSFitness)
    }
    xFoodFitness <- calcFitness(FUN, optimType, matrix(Xfood, ncol = numVar))
    Cfood <- 2*(1-t/maxIter)
    Kifood <- sapply(CSFitness, function(x){
      kijKH(x, xFoodFitness, bestFitness, worstFitness)
    })
    Xifood <- apply(candidateSolution, c(1), function(x){
      xijKH(x, Xfood, epsilon)
    })

    Kibest <- sapply(CSFitness, function(x){
      kijKH(x, bestFitness, bestFitness, worstFitness)
    })
    Xibest <- apply(candidateSolution, c(1), function(x){
      xijKH(x, best, epsilon)
    })

    if(numVar == 1){
      betaFood <- Cfood*Kifood*Xifood
      betaBest <- Xibest * Kibest
    }else{
      betaFood <- t(Cfood*Kifood*Xifood)
      betaBest <- t(Xibest) * Kibest
    }
    beta <- betaFood + betaBest

    f <- foragingSpeed*beta + inertiaWeightOfForagingSpeed*f

    # physical difussion
    D <- maxDifussionSpeed * (1 - t/maxIter)*runif(1, min = -1, max = 1)

    # Motion calculation
    TotalMotion <- N + f + D
    deltaT <- constantSpace*sum(upperBound - lowerBound)
    candidateSolution <- candidateSolution + deltaT * TotalMotion

    # implement genetic operator ----

    # CrossOver
    CSFitness <- calcFitness(FUN, optimType, candidateSolution)
    best <- candidateSolution[order(CSFitness)[1], ]
    worst <- candidateSolution[order(CSFitness)[numPopulation], ]
    bestFitness <- calcFitness(FUN, optimType, matrix(best, ncol = numVar))
    worstFitness <- calcFitness(FUN, optimType, matrix(worst, ncol = numVar))
    gbest <- calcBest(FUN, -1*optimType, rbind(candidateSolution, gbest))
    Kibest <- sapply(CSFitness, function(x){
      kijKH(x, bestFitness, bestFitness, worstFitness)
    })
    randomMatrix <- matrix(runif(numVar * numPopulation), ncol = numVar)
    Cr <- 0.2 * Kibest
    prob <- randomMatrix < Cr
    if(!all(prob == FALSE)){
      Xrm <- sapply(col(candidateSolution)[prob], function(x){
        choosen <- sample.int(numPopulation, 1)
        return(candidateSolution[choosen, x])
      })
      candidateSolution[prob] <- Xrm
    }

    # Mutation
    CSFitness <- calcFitness(FUN, optimType, candidateSolution)
    best <- candidateSolution[order(CSFitness)[1], ]
    worst <- candidateSolution[order(CSFitness)[numPopulation], ]
    bestFitness <- calcFitness(FUN, optimType, matrix(best, ncol = numVar))
    worstFitness <- calcFitness(FUN, optimType, matrix(worst, ncol = numVar))
    gbest <- calcBest(FUN, -1*optimType, rbind(candidateSolution, gbest))
    Kibest <- sapply(CSFitness, function(x){
      kijKH(x, bestFitness, bestFitness, worstFitness)
    })
    randomMatrix <- matrix(runif(numVar * numPopulation), ncol = numVar)
    Mu <- 0.05 * Kibest
    prob <- randomMatrix < Mu
    if(!all(prob == FALSE)){
      Xgbest <- sapply(col(candidateSolution)[prob], function(x){
        P <- sample.int(numPopulation, 1)
        Q <- sample.int(numPopulation, 1)
        return(gbest[x] + mu * (candidateSolution[P, x] - candidateSolution[Q, x]))
      })
      candidateSolution[prob] <- Xgbest
    }

    setTxtProgressBar(progressbar, t)
  }
  close(progressbar)
  gbest <- calcBest(FUN, -1*optimType, rbind(candidateSolution, gbest))
  gbest <- checkBound(gbest, lowerBound, upperBound)
  return(unname(gbest))
}


xijKH <- function(i, j, epsilon){
  (j - i)/(dist(rbind(j, i)) + epsilon)
}

kijKH <- function(i, j, best, worst){
  if(worst == best){
    0
  }else{
    (i - j)/(worst - best)
  }
}
