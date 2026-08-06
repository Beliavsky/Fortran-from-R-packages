# ------------------------------------------------------------------------------
# ----- Step2_ScenarioGeneration.R ---------------------------------------------
# ------------------------------------------------------------------------------

#' Generate Index Scenerio
#' 
#' @description  Simulate a 3D array, numScen-by-numStep-by-numIndex, of Black-Scholes return
#' factors for numIndex indices in each of numStep time steps and each of
#' numScen scenarios. Covariances among indices are specified in covMatrix.
#' Stepsize is given is dT and interpolated discount factors are given in vDF.
#' Random seed is optional for reproducibility.
#'
#' @param covMatrix A numIndex-by-numIndex matrix of doubles of covariances
#'   among numIndex indices.
#' @param numScen An integer of number of scenario (sample paths) to be
#'   simulated.
#' @param numStep An integer of number of periods to be simulated.
#' @param indexNames A vector of strings containing index names.
#' @param dT A double of stepsize in years; dT = 1 / 12 would be monthly.
#' @param forwardCurve A vector of doubles of discount rates at each time step.
#' @param seed An integer of the deterministic seed for random sampling.
#' @return Outputs a 3D array (numScen-by-numStep-by-numIndex) of index
#'   scenarios
#' @examples
#' genIndexScen(mCov, 100, 360, indexNames, 1 / 12, cForwardCurve, 1)
#' @export
#' @importFrom stats rnorm
genIndexScen <- function(covMatrix, numScen, numStep,
                           indexNames, dT = 1 / 12, forwardCurve, seed){
    if (!missing(seed)) {
        set.seed(seed)  # fix random seed if it is given
    }

    if (numScen < 1 || numStep < 1) {
      stop("Please input a numScen and numStep greater than zero")
    }

    numIndex <- length(indexNames) # no. of indices
    # Cholesky decomposed covariance matrix among indices
    cd <- chol(covMatrix)
    stdFactor <- rowSums(cd ^ 2) / 2       # std.dev factor in BS growth factor
    dtSqrt <- sqrt(dT)

    # construct matrices for faster computations
    stdFactorMat <- matrix(rep(stdFactor, each = numStep), nrow = numStep)
    fcMat <- matrix(rep(forwardCurve, times = numIndex), ncol = numIndex)
    mmu <- (fcMat - stdFactorMat) * dT

    # 3D array of scenarios, numScen-by numStep-by-numIndex
    indexScen <- array(0, dim = c(numScen, numStep, numIndex))

    for (i in 1:numScen) {
        # numIndex by numStep matrix of temporary standard normal rv's
        temp <- matrix(rnorm(numIndex * numStep), nrow = numStep,
                       ncol = numIndex)

        # Black-Scholes return factors, vectorized implementation for speed
        tempMat <- exp(mmu + dtSqrt * (temp %*% cd))
        indexScen[i, , ] <- tempMat
    }
    return (indexScen)
}
#----- random covariance matrix
#' @importFrom stats runif
dim <- 5
covMatrix <- matrix(runif(dim ^ 2), nrow = dim)
covMatrix <- t(covMatrix) %*% covMatrix
while (min(eigen(covMatrix)$values) < 1e-6) {
    covMatrix <- matrix(runif(dim ^ 2), nrow = dim)
    covMatrix <- t(covMatrix) %*% covMatrix
}

# ------------------------------------------------------------------------------
#' @importFrom stats runif
rFundMap <- function(indexNames, numFund){
    # Simulate a random fund map that maps the indices to numFund funds.
    # The first numIndex rows of the fund map is an identity matrix.
    # Inputs:
    # -- indexNames: vector of strings, names of indices
    # -- numFund: integer, number of funds
    # Outputs:
    # -- fundMap: (numIndex + numFund)-by-numIndex of doubles

    numIndex <- length(indexNames)

    # placeholder of fundMap
    fundMap <- mat.or.vec(nr = numFund, nc = numIndex)
    colnames(fundMap) <- indexNames

    # randomly assign the number of indices that a fund will be mapped with
    indexMap <- ceiling(runif(numFund, min = 1, max = numIndex))
    # randomly pick indexMap from the pool of indices and assign random weights
    for (i in 1:numFund) {
        # select index
        select <- sample(indexNames, indexMap[i], replace = FALSE)
        # generate weight
        fundMap[i, select] <- sample(c(1:10), length(select), replace = TRUE)
        fundMap[i, ] <- fundMap[i, ] / sum(fundMap[i, ])
    }
    # combine identity matrix to fund map

    return (fundMap = rbind(diag(numIndex), fundMap))
}

# ------------------------------------------------------------------------------
#' 
#' Generate Fund Scenerio
#' 
#' @description Calculate numScen-by-numStep-by-numFund fund scenarios based on given index
#' scenarios indexScen and fund map fundMap that maps indices to funds.
#'
#' @param fundMap A numFund-by-numIndex matrix of doubles,
#'   mapping indices to funds.
#' @param indexScen A numScen-by-numStep-by-numIndex array of doubles,
#'   index scenarios.
#' @return Outputs a numScen-by-numStep-by-numFund array of doubles of
#'   fund scenarios.
#' @examples
#' genFundScen(fundMap, indexScen)
#' @export
genFundScen <- function(fundMap, indexScen){
    # Calculate numScen-by-numIndex-by-numStep fund scenarios based on
    # given index scenarios indexScen and fund map fundMap that maps indices
    # to funds

    # extract dimensions from given input
    scenDim <- dim(indexScen)
    numFund <- nrow(fundMap)
    if (length(scenDim) == 3) {
        if (dim(fundMap)[2] != dim(indexScen)[3]) {
          stop("Funds from indexScen should align with funds from fundMap")
        }

        numScen <- scenDim[1]
        numStep <- scenDim[2]

        # placeholder for the fund scenarios, numScen-by-numStep-by-numFund
        fundScen <- array(0, dim = c(numScen, numStep, numFund))

        for (j in 1:numScen) {
            # the j-th scenario for all indices
            scenJ <- indexScen[j, , ]
            # convert an index scen. to fund scen.
            fundScen[j, , ] <- t(fundMap %*% t(scenJ))
        }
    } else if (length(scenDim) == 2) {
        if (dim(fundMap)[2] != dim(indexScen)[2]) {
          stop("Funds from indexScen should align with funds from fundMap")
        }
        # convert an index scen. to fund scen.
        fundScen <- t(fundMap %*% t(indexScen))
    }

    return (fundScen)
}
