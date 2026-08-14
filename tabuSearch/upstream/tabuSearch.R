#' R based tabu search algorithm for binary configurations
#'
#' @param size the length of the binary configuration.
#' @param iters the number of iterations in the preliminary search of the algorithm.
#' @param objFunc a user supplied method that evaluates the objective function for a given binary string. The objective function is required to take as an argument a vector of zeros and ones.
#' @param config a starting configuration.
#' @param neigh a number of neighbour configurations to check at each iteration. The default is all, which is the length of the string. If \code{neigh} < \code{size}, the neighbours are chosen at random.
#' @param listSize tabu list size.
#' @param nRestarts the maximum number of restarts in the intensification stage of the algorithm.
#' @param repeatAll the number of times to repeat the search.
#' @param verbose if true, the name of the current stage of the algorithm is printed e.g. preliminary stage, intensification stage, diversification stage.
#' @returns A tabu object which is a list giving:
#' * `configKeep` a matrix of configurations chosen at each iteration of the algorithm.
#' * `eUtilityKeep` value of the objective function for the above configurations.
#' * `iters` the number of iterations in the preliminary search of the algorithm.
#' * `neigh` a number of neighbour configurations checked at each iteration.
#' * `listSize` tabu list size.
#' * `repeatAll` the number of times the search was repeated.
#' @export
#' @description A tabu search algorithm for optimizing binary strings. It takes a user defined objective function and reports the best binary configuration found throughout the search i.e. the one with the highest objective function value. The algorithm can be used for variable selection. The results can be plotted and summarized using \code{plot.tabu} and  \code{summary.tabu}.
#' @seealso [plot.tabu()] [summary.tabu()]
#' @references Glover, F., 1977. Heuristics for integer programming using surrogate constraints. Decision Sciences 8, 156-166.
#' @references Glover, F., 1986. Future paths for integer programming and links to artificial intelligence. Computers and Operations Research 13, 533-549.
#' @references Fouskakis, D., Draper, D., 2002. Stochastic optimization: a review. International Statistical Review 70, 315-349.
#' @author K. Domijan
#' @examples
#' # A simple example
#' evaluateSimple <- function(th)return(1)
#' result <- tabuSearch(size = 20, iters = 100, objFunc = evaluateSimple)
#' \dontrun{
#' # simulate 10-d data: 150 samples from 3 bivariate normals and 8 noise variables.
#' # Variable selection should recover the first two variables
#' library(MASS)
#' NF <- 10
#' G <- 3
#' NTR <- 50
#' NTE <- 50
#'  muA <- c(1,3)
#'  SigmaA <- matrix(c(0.2, 0.04, 0.04, 0.2), 2, 2)
#'  muB <- c(1.2,1)
#'    SigmaB <- matrix(c(0.1, -0.06, 0.004, 0.2), 2, 2)
#'    muC <- c(3,2)
#' SigmaC <- matrix(c(0.2, 0.004, 0.004, 0.2), 2, 2)
#'
#' train <- rbind(mvrnorm(NTR, muA, SigmaA), mvrnorm(NTR, muB, SigmaB), mvrnorm(NTR, muC, SigmaC))
#' test <- rbind(mvrnorm(NTE, muA, SigmaA), mvrnorm(NTE, muB, SigmaB), mvrnorm(NTE, muC, SigmaC))
#'
#' train <- cbind(train, matrix(runif(G * NTR * (NF - 2), 0, 4), nrow = G * NTR, ncol = (NF-2)))
#' test <- cbind(test, matrix(runif(G * NTE * (NF - 2), 0, 4), nrow = G * NTE, ncol = (NF-2)))
#'
#' wtr <-  as.factor(rep(1:G, each = NTR))
#' wte <-  as.factor(rep(1:G, each = NTE))
#' pairs(train, col = wtr)
#'
#'
#' library(e1071)
#'
#' evaluate <- function(th){
#'   if (sum(th) == 0)return(0)
#'  model <- svm(train[ ,th==1], wtr , gamma = 0.1)
#'  pred <- predict(model, test[ ,th==1])
#'  csRate <- sum(pred == wte)/NTE
#'  penalty <- (NF - sum(th))/NF
#'  return(csRate + penalty)
#' }
#' res <- tabuSearch(size = NF, iters = 50, objFunc = evaluate, config = matrix(1,1,NF),
#'                   listSize = 5, nRestarts = 4)
#' plot(res)
#' plot(res, "tracePlot")
#' summary(res, verbose = TRUE)
#' }





tabuSearch <-
function(size = 10, iters = 100, objFunc = NULL, config = NULL, neigh = size, listSize = 9, nRestarts = 10,
 repeatAll = 1, verbose = FALSE){
    if (size < 2){
        stop("error: config too short!")
    }
    if (iters < 2){
        stop("error: not enough iterations!")
    }
    if (listSize >= size){
        stop("error: listSize too big!")
    }
    if (neigh > size){
        stop("error: too many neighbours!")
    }
    if (is.null(objFunc)) {
        stop("A evaluation function must be provided. See the objFunc parameter.")
    }
    if (is.null(config)) {
        config <- matrix(0, 1, size)
        config[sample(1:size, sample(1:size, 1))] <- 1
    }
    else if (size != length(config)){
        stop("Length of the starting configuration != size")
    }
    if (repeatAll < 1){
        stop("error: repeatAll must be > 0")
    }

    iter <- 1
    configKeep <- matrix(, repeatAll * iters * (nRestarts + 3), size)
    eUtilityKeep <- vector(, repeatAll * iters * (nRestarts + 3))

    for(j in 1:repeatAll){

        if (j > 1) {
            config <- matrix(0, 1, size)
            config[sample(1:size, sample(1:size, 1))] <- 1
        }

        tabuList <- matrix(0, 1, size)
        listOrder <- matrix(0, 1, listSize) #remembers order of tabu moves

        eUtility <- objFunc(config)
        aspiration <- eUtility    #highest utility found so far


        preliminarySearch<- function(){

            configKeep[iter, ] <- config
            eUtilityKeep[iter] <- eUtility
            iter <- iter + 1

            for (i in 2:iters){
                neighboursEUtility <- matrix(0, 1, size)
                configTemp <- t(matrix(config, size, neigh))
                randomNeighbours <- sample(size, neigh) #pick random neighbours
                diag(configTemp[, randomNeighbours]) <- abs(diag(configTemp[, randomNeighbours]) - 1)#flip
                neighboursEUtility[randomNeighbours] <- apply(configTemp, 1, objFunc)
                maxNontaboo <- max(neighboursEUtility[tabuList == 0])
                maxTaboo <- max(neighboursEUtility[tabuList == 1], -Inf)

                #find new move
                move <- ifelse(maxTaboo > maxNontaboo & maxTaboo > aspiration,
                ifelse(length(which(neighboursEUtility == maxTaboo)) == 1,
                       which(neighboursEUtility == maxTaboo), sample(which(neighboursEUtility == maxTaboo), 1)),
                ifelse(length(which(neighboursEUtility == maxNontaboo & tabuList == 0)) == 1,
                       which(neighboursEUtility == maxNontaboo & tabuList == 0), sample(which(neighboursEUtility == maxNontaboo & tabuList == 0), 1)))

                #if new utility is lower than old, add move to tabu list, else adjust aspiration, if necessary
                if (eUtility >= neighboursEUtility[move]){
                    tabuList[move] <- 1
                    if(sum(tabuList) > listSize){ #if tabu list is full
                        tabuList[listOrder[1]] <- 0
                        listOrder[1:listSize] <- c(listOrder[2:listSize], 0)
                        }
                    listOrder[min(which(listOrder == 0))] <- move
                    }
                else if(neighboursEUtility[move] > aspiration) aspiration <- neighboursEUtility[move]

                #make the new move
                eUtility <- neighboursEUtility[move]
                config[move] <- abs(config[move]-1)
                configKeep[iter,] <- config
                eUtilityKeep[iter] <- eUtility
                iter <- iter + 1
            }

            result = list(aspiration = aspiration,  configKeep = configKeep, eUtilityKeep = eUtilityKeep, iter = iter)
            return(result)
        }


        #PRELIMINARY SEARCH
        if (verbose) cat("Preliminary search stage...\n")
        result <- preliminarySearch()
        aspiration <- result$aspiration
        configKeep <- result$configKeep
        eUtilityKeep <- result$eUtilityKeep
        iter <- result$iter

        #INTENSIFICATION
    	  tempo <- 0
    	  restarts <- 0
        while(tempo < aspiration & restarts < nRestarts){
            if (verbose) cat("Intensification stage...\n")
            eUtility <- max(eUtilityKeep)
            tempo <- aspiration
            config <- configKeep[max(which(eUtilityKeep == max(eUtilityKeep))), ]
            result <- preliminarySearch()
            aspiration <- result$aspiration
            configKeep <- result$configKeep
            eUtilityKeep <- result$eUtilityKeep
            iter <- result$iter
            restarts <- restarts + 1
        }


        #DIVERSIFICATION
        if (verbose) cat("Diversification stage...\n")
        config <- matrix(0, 1, size)
        config[sample(1:size, sample(1:size, 1))] <- 1
        eUtility <- objFunc(config)

        frequent <- apply(configKeep, 2, function(x)sum(diff(x) != 0))   #create new tabu list from most frequent moves
        tabuList <- as.numeric(rank(frequent, ties.method =  "random") > (size - listSize))
        listOrder <- sample(which(rank(frequent, ties.method =  "random") > (size - listSize)), listSize)

        result <- preliminarySearch()
        iter <- result$iter
        configKeep <- result$configKeep
        eUtilityKeep <- result$eUtilityKeep
    }

    endResult <- list(type = "binary configuration", configKeep = configKeep[1:(iter - 1), ],
    eUtilityKeep = eUtilityKeep[1:(iter - 1)], iters = iters, neigh = neigh,
    listSize = listSize,  repeatAll = repeatAll)
    class(endResult) = "tabu"
    return(endResult)
}
