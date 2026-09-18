##  roxygen2::roxygenise("C:/users/valen/onedrive/myrepo/r/tclust", load_code=roxygen2:::load_installed)

#'
#' TKMEANS method for robust K-means clustering
#' 
#' @name tkmeans
#' @aliases print.tkmeans
#' @description This function searches for \code{k} (or less) spherical clusters 
#'  in a data matrix \code{x}, whereas the \code{ceiling(alpha n)} most outlying 
#'  observations are trimmed.
#' @param x A matrix or data.frame of dimension n x p, containing the observations (row-wise). 
#' @param k The number of clusters initially searched for.
#' @param alpha The proportion of observations to be trimmed.
#' @param nstart The number of random initializations to be performed.
#' @param niter1 The number of concentration steps to be performed for the nstart initializations.
#' @param niter2 The maximum number of concentration steps to be performed for the 
#'  \code{nkeep} solutions kept for further iteration. The concentration steps are 
#'  stopped, whenever two consecutive steps lead to the same data partition.
#' @param nkeep The number of iterated initializations (after niter1 concentration 
#'  steps) with the best values in the target function that are kept for further iterations
#' @param iter.max (deprecated, use the combination \code{nkeep, niter1 and niter2}) 
#'  The maximum number of concentration steps to be performed.
#'  The concentration steps are stopped, whenever two consecutive steps lead
#'  to the same data partition.
#' @param zero_tol The zero tolerance used. By default set to 1e-16.
#' @param points Optional initial mean vectors, \code{NULL} or a matrix with \code{k} 
#'  vectors used as means to initialize the algorithm. If initial mean vectors are 
#'  specified, \code{nstart} should be 1 (otherwise the same initial means are 
#'  used for all runs).
#' @param center Optional centering of the data: a function or a vector of length p 
#'  which can optionally be specified for centering x before calculation
#' @param scale Optional scaling of the data: a function or a vector of length p 
#'  which can optionally be specified for scaling x before calculation
#' @param store_x A logical value, specifying whether the data matrix \code{x} shall be 
#'  included in the result object. By default this value is set to \code{TRUE}, because 
#'  some of the plotting functions depend on this information. However, when big data 
#'  matrices are handled, the result object's size can be decreased noticeably 
#'  when setting this parameter to \code{FALSE}.
#' @param parallel A logical value, specifying whether the nstart initializations should be done in parallel.
#' @param n.cores The number of cores to use when paralellizing, only taken into account if parallel=TRUE.
#' @param drop.empty.clust Logical value specifying, whether empty clusters shall be 
#'  omitted in the resulting object. (The result structure does not contain center 
#'  estimates of empty clusters anymore. Cluster names are reassigned 
#'  such that the first l clusters (l <= k) always have at least one observation.
#' @param trace Defines the tracing level, which is set to 0 by default. Tracing level 1 
#'  gives additional information on the stage of the iterative process.
#'
#' @return The function returns the following values:
#' \itemize{
#'     \item cluster - A numerical vector of size \code{n} containing the cluster assignment 
#'          for each observation. Cluster names are integer numbers from 1 to k, 0 indicates 
#'          trimmed observations. Note that it could be empty clusters with no observations 
#'          when \code{equal.weights=FALSE}.
#'     \item obj - The value of the objective function of the best (returned) solution.
#'     \item size - An integer vector of size k, returning the number of observations contained by each cluster.
#'     \item centers - A matrix of size p x k containing the centers (column-wise) of each cluster. 
#'     \item code - A numerical value indicating if the concentration steps have 
#'          converged for the returned solution (2).
#'     \item cluster.ini - A matrix with nstart rows and number of columns equal to 
#'          the number of observations and where each row shows the final clustering 
#'          assignments (0 for trimmed observations) obtained after the \code{niter1} 
#'          iteration of the \code{nstart} random initializations.
#'     \item obj.ini - A numerical vector of length \code{nstart} containing the values 
#'          of the target function obtained after the \code{niter1} iteration of the 
#'          \code{nstart} random initializations.
#'     \item x - The input data set.
#'     \item k - The input number of clusters.
#'     \item alpha - The input trimming level.
#' }
#' 
#' @author Valentin Todorov, Luis Angel Garcia Escudero, Agustin Mayo Iscar.
#' 
#' @references 
#' 
#' Cuesta-Albertos, J. A.; Gordaliza, A. and Matrán, C. (1997), "Trimmed k-means: 
#'  an attempt to robustify quantizers". Annals of Statistics, Vol. 25 (2), 553-576. 
#' 
#' @export
#'
#' @examples
#' 
#'  \dontshow{
#'      set.seed(0)
#'  }
#'  ##--- EXAMPLE 1 ------------------------------------------
#'  sig <- diag(2)
#'  cen <- rep(1,2)
#'  x <- rbind(MASS::mvrnorm(360, cen * 0,   sig),
#'             MASS::mvrnorm(540, cen * 5,   sig),
#'             MASS::mvrnorm(100, cen * 2.5, sig))
#'  
#'  ## Two groups and 10\% trimming level
#'  (clus <- tkmeans(x, k = 2, alpha = 0.1))
#'
#'  plot(clus)
#'  plot(clus, labels = "observation")
#'  plot(clus, labels = "cluster")
#'
#'  #--- EXAMPLE 2 ------------------------------------------
#'  data(geyser2)
#'  (clus <- tkmeans(geyser2, k = 3, alpha = 0.03))
#'  plot(clus)
#'  
tkmeans <- function(x, k, alpha=0.05, nstart=500, niter1=3, niter2=20, nkeep=5, iter.max, points=NULL, 
                   center=FALSE, scale=FALSE, store_x=TRUE, 
                   parallel=FALSE, n.cores=-1, 
                   zero_tol=1e-16, drop.empty.clust=TRUE, trace=0)  {
    
    if(!missing(iter.max)) {
        warning("The parameter 'iter.max' is deprecated, please read the help and use the combination of 'niter1', 'niter2', 'nkeep'.")
        niter1 <- iter.max
    }
        

    if(!is.null(points)) {
        if(nstart > 1)
            warning("If initial mean vectors are specified, 'nstart' should be 1 (otherwise the same initial means are used for all runs)")
        nstart <- nkeep <- 1
        parallel <- FALSE
    }
        
	parlist <- list(k=k, alpha=alpha, nstart=nstart, niter1=niter1, niter2=niter2, nkeep=nkeep, 
        center=center, scale=scale,
        zero_tol=zero_tol, drop.empty.clust=drop.empty.clust, trace=trace, store_x=store_x,
        equal.weights=TRUE)
              
    # Initial checks
    
    if(is.data.frame(x))
        x <- (data.matrix(x))
    if(!is.matrix (x))
        x <- as.matrix(x)
    if(any(is.na(x)))
        stop ("x cannot contain NA")
    if(!is.numeric (x)) 
        stop ("Parameter x: numeric matrix/vector expected")
        
    scaled <- myscale(x, center=center, scale=scale)
    x <- scaled$x
    if(store_x)
        parlist$x <- x     
  
      if (k < 1) 
        stop ("Parameter k: must be >= 1")
      if (alpha < 0 || alpha > 1) 
        stop ("Parameter alpha: must be in [0,1]")
      if (niter1 < 0 || as.integer(niter1) != niter1) 
        stop ("Parameter niter1: must be an integer >= 0")
      if (niter2 < 0 || as.integer(niter2) != niter2) 
        stop ("Parameter niter2: must be an integer >= 0")
      if (nkeep < 0 || as.integer(nkeep) != nkeep) 
        stop ("Parameter nkeep: must be an integer >= 0")
      if (nkeep > nstart) 
        stop ("Parameter nkeep: must be <= nstart")
      if(!is.logical(parallel))
        stop ("Parameter parallel: must be a logical TRUE or FALSE")
      if (n.cores < -2 || as.integer(n.cores) != n.cores) 
        stop ("Parameter n.cores: must be an integer >= 0, -1 or -2")
      if(zero_tol < 0)
        stop ("Parameter zero_tol: must be >= 0")
      if(trace != 0 & trace !=1)
        stop ("Parameter trace: must be 0 or 1")

    ###
    # FIRST STEP: get nstart solutions with niter1 concentration steps
    ###
    if(trace){
        cat(paste("\nPhase 1: obtaining ", nstart, " solutions.\n", sep = ""))
        if(parallel) cat(paste("\n Parallelizing initializations using ", n.cores, 
            " cores. Progress bar will not display accurate information.\n", sep = ""))
        pb <- txtProgressBar(min = 0, max = nstart, style = 3)
    }
  
    if(!parallel) {
        cluster.ini <- vector("list", nstart)   ##	for containing the best values for the parameters after several random starts
        obj.ini <- rep(0, nstart)               ## for containing best objective values
    
        for(j in 1:nstart) {
            assig_obj <- tkmeans_c1(x, k, alpha, niter1, zero_tol, points)                                                                   # niter1 steps!
            cluster.ini[[j]] <- assig_obj$cluster
            obj.ini[j] <- assig_obj$obj
      
            if(trace){
                setTxtProgressBar(pb, j)
            }
        }
    } else {
    
        # Setup parallel cluster
        if(n.cores == -1){
            n.cores <- detectCores()
        } else if (n.cores == -2){
            n.cores <- detectCores() - 1
        }
        parclus <- makeCluster(n.cores)
        registerDoParallel(parclus)

        count <- 0
        comb <- function(...) { # Custom combination function for the foreach loop
            count <<- count + length(list(...)) - 1
            setTxtProgressBar(pb, count)
            flush.console()
            c(...) # this can feed into .combine option of foreach
        }
        
        init.results <- foreach(j = 1:nstart,
                                .packages = "tclust",
                                .combine = ifelse(trace, "comb", "c"),
                                .multicombine = TRUE,
                                .inorder = FALSE) %dopar% {
            assig_obj <- tkmeans_c1(x, k, alpha, niter1, zero_tol)
            list(assig_obj$cluster, assig_obj$obj)
        }

        stopCluster(parclus)

        cluster.ini <- init.results[0:(nstart-1) * 2 + 1]       # Impair positions of cluster.ini
        obj.ini <- unlist(init.results[1:nstart * 2])           # Pair positions of cluster.ini
    }
  
    ###
    # SECOND STEP: get nkeep best solutions so far
    ###
    if(trace){
        cat(paste("\n\nPhase 2: obtaining ", nkeep, " best solutions out of the intial ", nstart," solutions.\n", sep = ""))
    }
    best_index  <- order(obj.ini)[1:nkeep]
    best_assig_list <- cluster.ini[best_index]

  
    ###
    # THIRD STEP: apply niter2 concentration steps to nkeep best solutions, return the best solution
    ###
    if(trace){
        cat(paste("\nPhase 3: applying ", niter2, " concentration steps to each of the ", nkeep," best solutions.\n", sep = ""))
        pb2 <- txtProgressBar(min = 0, max = nkeep, style = 3)
    }
    
    best_iter <- NULL
    best_iter_obj <- Inf
  
    for(j in 1:nkeep){
        iter <- tkmeans_c2(x, k, best_assig_list[[j]], alpha, niter2, zero_tol=1e-16)
    
        if(iter$obj <= best_iter_obj){
            best_iter <- iter
            best_iter_obj <- iter$obj
        }
        
        if(trace){
            setTxtProgressBar(pb2, j)
        }
    }
    if(trace){
        cat("\n\n")
    }
  
    ## Adjust the returned object to be similar to 'tkmeans': 
    
    ## Handle empty clusters
    if(drop.empty.clust)
        idxuse <- which(best_iter$size > 0)
    else
        idxuse <- 1:k
    idxuse <- idxuse[order(best_iter$size[idxuse], decreasing = TRUE)]
    
    ClusterIDs <- rep(0, k + 1)
    ClusterIDs[idxuse + 1] <- 1:length(idxuse)

    k.real <- length(idxuse)
    
    ##  - A list of values internally used by functions related to tclust objects.
    int <- list(
    	iter.successful=0, iter.converged=0,
    	dim=dim(x))
    
    ## - Transpose the 'centers' matrix
    best_iter$centers <- t(best_iter$centers)
    best_iter$centers <- best_iter$centers[, idxuse, drop = FALSE]
    
    best_iter$cluster <- as.vector(best_iter$cluster)
    best_iter$cluster <- ClusterIDs[best_iter$cluster + 1]
    best_iter$size <- as.vector(best_iter$size)
    best_iter$size <- best_iter$size[idxuse]

    best_iter$cov <- NULL
    best_iter$posterior <- NULL

    best_iter$NlogL <- NULL     
    best_iter$disttom <- NULL
    
    ret <- c(best_iter, list(cluster.ini=matrix(unlist(cluster.ini), byrow=TRUE, nrow=length(cluster.ini)),
        obj.ini=obj.ini, int=int, par=parlist, k=sum(best_iter$size > 0)))

    ##  - Dimnames of center and cov
    dn.x <- dimnames(x)
    rownames(ret$centers) <- if(is.null(colnames(x))) paste ("X", 1:ncol(x)) else colnames (x)
    colnames(ret$centers) <- paste("C", 1:k.real)

    ## Back transform centers, cov and x
    cmm <- scaled$scale %*% t (scaled$scale)
    for(i in 1:ret$k) {
        ret$centers[,i] <- ret$centers[,i] * scaled$scale + scaled$center
    }
    x <- sweep(x, 2L, scaled$scale, `*`, check.margin = FALSE)
    x <- sweep(x, 2L, scaled$center, `+`, check.margin = FALSE)
    if(store_x)
        ret$par$x <- x
    
    class(ret) <- "tkmeans"    
    return(ret)
}

