##  roxygen2::roxygenise("C:/users/valen/onedrive/myrepo/r/tclust", load_code=roxygen2:::load_installed)

#'
#' TCLUST method for robust clustering
#' 
#' @name tclust
#' @aliases print.tclust
#' @description This function searches for \code{k} (or less) clusters with 
#'  different covariance structures in a data matrix \code{x}. Relative cluster 
#'  scatter can be restricted when \code{restr="eigen"} by constraining the ratio 
#'  between the largest and the smallest of the scatter matrices eigenvalues 
#'  by a constant value \code{restr.fact}. Relative cluster scatters can be also 
#'  restricted with \code{restr="deter"} by constraining the ratio between the 
#'  largest and the smallest of the scatter matrices' determinants. 
#'
#'  For robustifying the estimation, a proportion \code{alpha} of observations is trimmed. 
#'  In particular, the trimmed k-means method is represented by the \code{tclust()} method,
#'  by setting parameters \code{restr.fact=1}, \code{opt="HARD"} and \code{equal.weights=TRUE}. 
#'
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
#' @param equal.weights A logical value, specifying whether equal cluster weights 
#'  shall be considered in the concentration and assignment steps.
#' @param restr Restriction type to control relative cluster scatters. 
#'  The default value is \code{restr="eigen"}, so that the maximal ratio between 
#'  the largest and the smallest of the scatter matrices eigenvalues is constrained 
#'  to be smaller then or equal to \code{restr.fact} 
#'  (Garcia-Escudero, Gordaliza, Matran, and Mayo-Iscar, 2008). 
#'  Alternatively, \code{restr="deter"} imposes that the maximal ratio between 
#'  the largest and the smallest of the scatter matrices determinants is smaller 
#'  or equal than \code{restr.fact} (see Garcia-Escudero, Mayo-Iscar and Riani, 2020).
#'  As a third option the 14 Gaussian Parsimonious Clustering Models (GPCM) can be applied
#'      by setting \code{restr="GPCM"}. All necesary parameters are taken from 
#'  the optional parameter \code{pars}.
#'
#' @param restr.fact The constant \code{restr.fact >= 1} constrains the allowed 
#'  differences among group scatters in terms of eigenvalues ratio
#'  (if \code{restr="eigen"}) or determinant ratios (if \code{restr="deter"}). Larger values 
#'  imply larger differences of group scatters, a value of 1 specifies the 
#'  strongest restriction.
#' @param cshape constraint to apply to the shape matrices, \code{cshape >= 1}, 
#'  (see Garcia-Escudero, Mayo-Iscar and Riani, 2020)). 
#'  This options only works if \code{restr=='deter'}. In this case the default 
#'  value is \code{cshape=1e10} to ensure the procedure is (virtually) affine equivariant. 
#'  On the other hand, \code{cshape} values close to 1 would force the clusters to 
#'  be almost spherical (without necessarily the same scatters if \code{restr.fact} 
#'  is strictly greater than 1).
#' @param pars a character string or a list containg the parameters for the GPCM model. The only 
#'   mandatory element of the list is \code{pars$pars} which specifies the type of Gaussian 
#'  Parsimonious Clustering Model and can be one of:
#'  \code{"VVE","EVE","VVV","EVV","VEE","EEE","VEV","EEV","VVI", "EVI","VEI","EEI","VII","EII"}.
#'  If \code{par} is a list it can have the following elements:
#'  \itemize{
#'     \item pars - type of the Gaussian Parsimonious Clustering Model
#'      \item cdet - restriction to be applied to the determinants, default is \code{cdet=100}
#'      \item shw - restriction to be applied to the elements of the shape matrices 
#'          inside each group, default is \code{shw=100}
#'      \item shb - restriction to be applied to the elements of the shape matrices 
#'          across each group, default is \code{shb=100}
#'      \item maxiterDSR - maximum number of iterations to obtain the requested restricted 
#'        determinants, shape matrices and rotation, default is \code{maxiterDRS=20}
#'      \item tolDSR - tolerance to use to exit the loop for obtaining the requested 
#'        restricted determinants, shape matrices and rotation, default is \code{maxiterDRS=20}
#'      \item maxiterS - maximum number of iterations to obtain the restricted 
#'          shape matrix, default is \code{maxiterS=20}
#'      \item tolS - tolerance to use to exit the iterative procedure for estimating 
#'          the shape, default is \code{tolS=1e-5}
#'      \item maxiterR- maximum number of iterations to obtain the common rotation 
#'          matrix in presence of varying shape, default is \code{maxiterR=20}
#'      \item tolR - tolerance to use to exit the iterations to obtain the common 
#'          rotation matrix in presence of varying shape, default is \code{tolR=1e-5}
#'      \item zerotol - tolerance value to declare all input values equal to 0 in 
#'          the eigenvalues restriction routine, default is \code{zerotol=1e-10}
#'  }
#' @param zero_tol The zero tolerance used. By default set to 1e-16.
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
#' @param n.cores The number of cores to use when paralellizing, only taken into account if \code{parallel=TRUE}.
#' @param opt Define the target function to be optimized. A classification likelihood 
#'  target function is considered if \code{opt="HARD"} and a mixture classification 
#'  likelihood if \code{opt="MIXT"}.
#' @param drop.empty.clust Logical value specifying, whether empty clusters shall be 
#'  omitted in the resulting object. (The result structure does not contain center 
#'  and covariance estimates of empty clusters anymore. Cluster names are reassigned 
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
#'     \item NlogL - A value related to the classification log-likelihood of the best (returned) solution. 
#'          If \code{opt=="HARD"}, \code{NlogL = -2*obj}.
#'     \item size - An integer vector of size k, returning the number of observations contained by each cluster.
#'     \item weights - Vector of Cluster weights
#'     \item centers - A matrix of size p x k containing the centers (column-wise) of each cluster. 
#'     \item cov - 	An array of size p x p x k containing the covariance matrices of each cluster. 
#'     \item code - A numerical value indicating if the concentration steps have 
#'          converged for the returned solution (2).
#'     \item posterior - A matrix with k columns that contains the posterior 
#'          probabilities of membership of each observation (row-wise) to the \code{k} 
#'          clusters. This posterior probabilities are 0-1 values in the 
#'          \code{opt="HARD"} case. Trimmed observations have 0 membership probabilities 
#'          to all clusters.
#'      \item{MIXMIX} - BIC which based on the parameters estimated through the mixture log-likelihood and 
#'          the maximized mixture likelihood as goodness of fit measure. This output 
#'          is present only  if \code{opt="MIXT"}.
#'      \item{MIXCLA} - BIC which uses the classification likelihood based on
#'          parameters estimated through the mixture likelihood (In some books 
#'          this quantity is called ICL). This output 
#'          is present only  if \code{opt="MIXT"}.
#'     \item{CLACLA} - BIC which uses the classification likelihood based on
#'          parameters estimated using the classification likelihood. This output 
#'          is present only  if \code{opt="HARD"}.
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
#' @details The procedure allows to deal with robust clustering with an \code{alpha}
#'  proportion of trimming level and searching for \code{k} clusters. We are considering 
#'  classification trimmed likelihood when using \code{opt=”HARD”} so that “hard” or “crisp” 
#'  clustering assignments are done. On the other hand, mixture trimmed likelihood 
#'  are applied when using \code{opt=”MIXT”} so providing a kind of clusters “posterior” 
#'  probabilities for the observations. 

#'  Relative cluster scatter can be restricted when \code{restr="eigen"} by constraining 
#'  the ratio between the largest and the smallest of the scatter matrices eigenvalues 
#'  by a constant value \code{restr.fact}. Setting \code{restr.fact=1}, yields the 
#'  strongest restriction, forcing all clusters to be spherical and equally scattered. 
#'  Relative cluster scatters can be also restricted with \code{restr="deter"} by 
#'  constraining the ratio between the largest and the smallest of the scatter 
#'  matrices' determinants. 
#'
#'  This iterative algorithm performs "concentration steps" to improve the current 
#'  cluster assignments. For approximately obtaining the global optimum, the procedure 
#'  is randomly initialized \code{nstart} times and \code{niter1} concentration steps are performed for 
#'  them. The \code{nkeep} most “promising” iterations, i.e. the \code{nkeep} iterated solutions with 
#'  the initial best values for the target function, are then iterated until convergence 
#'  or until \code{niter2} concentration steps are done. 
#'
#' The parameter \code{restr.fact} defines the cluster scatter matrices restrictions, 
#'  which are applied on all clusters during each concentration step. It restricts 
#'  the ratio between the maximum and minimum eigenvalue of 
#'  all clusters' covariance structures to that parameter. Setting \code{restr.fact=1}, 
#'  yields the strongest restriction, forcing all clusters to be spherical and equally scattered. 
#'
#' Cluster components with similar sizes are favoured when considering \code{equal.weights=TRUE} 
#'  while \code{equal.weights=FALSE} admits possible different prior probabilities for 
#'  the components and it can easily return empty clusters when the number of 
#'  clusters is greater than apparently needed.
#' 
#' @author Javier Crespo Guerrero, Luis Angel Garcia Escudero, Agustin Mayo Iscar.
#' 
#' @references 
#' 
#' Fritz, H.; Garcia-Escudero, L.A.; Mayo-Iscar, A. (2012), "tclust: An R Package 
#'  for a Trimming Approach to Cluster Analysis". Journal of Statistical Software, 
#'  47(12), 1-26. URL http://www.jstatsoft.org/v47/i12/
#' 
#' Garcia-Escudero, L.A.; Gordaliza, A.; Matran, C. and Mayo-Iscar, A. (2008), 
#'  "A General Trimming Approach to Robust Cluster Analysis". Annals of Statistics, 
#'  Vol.36, 1324--1345.  
#'
#' García-Escudero, L. A., Gordaliza, A. and Mayo-Íscar, A. (2014). A constrained 
#'  robust proposal for mixture modeling avoiding spurious solutions. 
#'  Advances in Data Analysis and Classification, 27--43. 
#' 
#' García-Escudero, L. A., and Mayo-Íscar, A. and Riani, M. (2020). Model-based 
#'  clustering with determinant-and-shape constraint. Statistics and Computing, 
#'  30, 1363--1380.] 
#' @export
#'
#' @examples
#' 
#'  \dontshow{
#'      set.seed (0)
#'  }
#'  ##--- EXAMPLE 1 ------------------------------------------
#'  sig <- diag(2)
#'  cen <- rep(1,2)
#'  x <- rbind(MASS::mvrnorm(360, cen * 0,   sig),
#'             MASS::mvrnorm(540, cen * 5,   sig * 6 - 2),
#'             MASS::mvrnorm(100, cen * 2.5, sig * 50))
#'  
#'  ## Two groups and 10\% trimming level
#'  clus <- tclust(x, k = 2, alpha = 0.1, restr.fact = 8)
#'  
#'  plot(clus)
#'  plot(clus, labels = "observation")
#'  plot(clus, labels = "cluster")
#'  
#'  ## Three groups (one of them very scattered) and 0\% trimming level
#'  clus <- tclust(x, k = 3, alpha=0.0, restr.fact = 100)
#'  
#'  plot(clus)
#'  
#'  ##--- EXAMPLE 2 ------------------------------------------
#'  data(geyser2)
#'  (clus <- tclust(geyser2, k = 3, alpha = 0.03))
#'  
#'  plot(clus)
#'  
#' \donttest{
#'
#'  ##--- EXAMPLE 3 ------------------------------------------
#'  data(M5data)
#'  x <- M5data[, 1:2]
#'  
#'  clus.a <- tclust(x, k = 3, alpha = 0.1, restr.fact =  1,
#'                    restr = "eigen", equal.weights = TRUE)
#'  clus.b <- tclust(x, k = 3, alpha = 0.1, restr.fact =  50,
#'                     restr = "eigen", equal.weights = FALSE)
#'  clus.c <- tclust(x, k = 3, alpha = 0.1, restr.fact =  1,
#'                    restr = "deter", equal.weights = TRUE)
#'  clus.d <- tclust(x, k = 3, alpha = 0.1, restr.fact = 50,
#'                    restr = "deter", equal.weights = FALSE)
#'  
#'  pa <- par(mfrow = c (2, 2))
#'  plot(clus.a, main = "(a)")
#'  plot(clus.b, main = "(b)")
#'  plot(clus.c, main = "(c)")
#'  plot(clus.d, main = "(d)")
#'  par(pa)
#'  
#'  ##--- EXAMPLE 4 ------------------------------------------
#'
#'  data (swissbank)
#'  ## Two clusters and 8\% trimming level
#'  (clus <- tclust(swissbank, k = 2, alpha = 0.08, restr.fact = 50))
#'  
#'  ## Pairs plot of the clustering solution
#'  pairs(swissbank, col = clus$cluster + 1)
#'  ## Two coordinates
#'  plot(swissbank[, 4], swissbank[, 6], col = clus$cluster + 1,
#'       xlab = "Distance of the inner frame to lower border",
#'       ylab = "Length of the diagonal")
#'  plot(clus)
#'  
#'  ## Three clusters and 0\% trimming level
#'  clus<- tclust(swissbank, k = 3, alpha = 0.0, restr.fact = 110)
#'  
#'  ## Pairs plot of the clustering solution
#'  pairs(swissbank, col = clus$cluster + 1)
#'  
#'  ## Two coordinates
#'  plot(swissbank[, 4], swissbank[, 6], col = clus$cluster + 1, 
#'        xlab = "Distance of the inner frame to lower border", 
#'        ylab = "Length of the diagonal")
#'  
#'  plot(clus)
#'  
#'  ##--- EXAMPLE 5 ------------------------------------------
#'   data(M5data)
#'   x <- M5data[, 1:2]
#'   
#'   ## Classification trimmed likelihood approach
#'   clus.a <- tclust(x, k = 3, alpha = 0.1, restr.fact =  50,
#'                      opt="HARD", restr = "eigen", equal.weights = FALSE)
#'  ## Mixture trimmed likelihood approach
#'   clus.b <- tclust(x, k = 3, alpha = 0.1, restr.fact =  50,
#'                      opt="MIXT", restr = "eigen", equal.weights = FALSE)
#'  
#'  ## Hard 0-1 cluster assignment (all 0 if trimmed unit)
#'  head(clus.a$posterior)
#'  
#'  ## Posterior probabilities cluster assignment for the
#'  ##  mixture approach (all 0 if trimmed unit)
#'  head(clus.b$posterior)
#'  
#'  ##--- EXAMPLE 6 ------------------------------------------
#'  sig <- diag(2)
#'  cen <- rep(1,2)
#'  x <- rbind(MASS::mvrnorm(360, cen * 0,   sig),
#'             MASS::mvrnorm(540, cen * 5,   sig * 6 - 2),
#'             MASS::mvrnorm(100, cen * 2.5, sig * 50))
#'  
#'  ## Two groups and 10\% trimming level
#'  (clus <- tclust(x, k = 2, alpha = 0.1, restr = "GPCM", pars = "EEE"))
#'  plot(clus, main="/p")
#'
#' }
#'

tclust <- function(x, k, alpha=0.05, nstart=500, niter1=3, niter2=20, nkeep=5, iter.max,
                   equal.weights=FALSE, restr=c("eigen", "deter", "GPCM"), restr.fact=12, cshape=1e10,
                   pars, 
                   opt=c("HARD", "MIXT"),
                   center=FALSE, scale=FALSE, store_x=TRUE, 
                   parallel=FALSE, n.cores=-1, 
                   zero_tol=1e-16, drop.empty.clust=TRUE, trace=0)  {
    
    restr <- match.arg(restr)
    opt <- match.arg(opt)
    restrC <- ifelse(restr=="GPCM", 1, 0)
    deterC <- restr == "deter"

    if(!missing(iter.max)) {
        warning("The parameter 'iter.max' is deprecated, please read the help and use the combination of 'niter1', 'niter2', 'nkeep'.")
        niter1 <- iter.max
    }
        
    if(equal.weights && opt == "MIXT") {
        warning("The parameter 'equalweights' cannot be TRUE if mixture model approach is assumed (opt='MIXT')")
        equal.weights <- FALSE
    }

    ## process 'pars': CPCM parameters, if any
    defaults <- list(
        p          = NA,
        k          = NA,
        pars       = "",        # type of Gaussian Parsimonious Clustering Model
        cdet       = 100,       # restriction to be applied to the determinants
        shw        = 100,       # restriction to be applied to the elements of the shape matrices inside each group
        shb        = 100,       # restriction to be applied to the elements of the shape matrices across each group
        maxiterDSR = 20,        # maximum number of iterations to obtain the requested restricted determinants, shape matrices and rotation
        tolDSR     = 1e-5,      # ?
        maxiterS   = 20,        # maximum number of iterations to obtain the restricted shape matrix
        tolS       = 1e-5,      # tolerance to use to exit the iterative procedure for estimating the shape
        maxiterR   = 20,        # maximum number of iterations to obtain the common rotation matrix in presence of varying shape
        tolR       = 1e-5,      # ?
        zerotol    = 1e-10,     # tolerance value to declare all input values equal to 0 in the eigenvalues restriction routine
        sortsh     = 0
    )
    
    if(missing(pars)) {
        if(restr == "GPCM")
            stop("For the GPCM model the parameter 'pars' cannot be missing!")
        else
            pars <- defaults
    }
    if(is.character(pars)) {
        if(length(pars) != 1)
            stop("If 'pars' is a character type, its length must be 1!")
        pars <- list(pars=pars)
    }
    pars <- utils::modifyList(defaults, pars)
    if(restr == "GPCM") {
        if(pars$pars == "")  
            stop("For the GPCM model the 'pars' entry in the pars list cannot be empty!")
        pars$pars <- match.arg(pars$pars, choices=c("VVE","EVE","VVV","EVV","VEE","EEE","VEV","EEV","VVI", "EVI","VEI","EEI","VII","EII"))
    }

	parlist <- list(k=k, alpha=alpha, nstart=nstart, niter1=niter1, niter2=niter2, nkeep=nkeep, 
        restr=restr, restr.C=restrC, deter.C=deterC, restr.fact=restr.fact, cshape=cshape, pars=pars, opt=opt,
        equal.weights=equal.weights, center=center, scale=scale,
#           fuzzy=fuzzy, m=m, 
        zero_tol=zero_tol, drop.empty.clust=drop.empty.clust, trace=trace, store_x=store_x)
              
    ## Initial checks
    
    if(is.data.frame(x))
        x <- (data.matrix(x))
    if(!is.matrix (x))
        x <- as.matrix(x)
    if(any(is.na(x)))
        stop ("x cannot contain NA")
    if(!is.numeric (x)) 
        stop ("Parameter x: numeric matrix/vector expected")

    parlist$n <- n <- nrow(x)
    parlist$p <- p <- ncol(x)
            
    scaled <- myscale(x, center=center, scale=scale)
    x <- scaled$x
    if(store_x)
        parlist$x <- x     
  
      if (!k >= 1) 
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
      if(!is.logical(equal.weights))
        stop ("Parameter equal.weights: must be a logical TRUE or FALSE")
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
    if(parallel) cat("\n Running parallel initializations.", 
        "\nProgress bar will not display accurate information.\n")
    pb <- txtProgressBar(min = 0, max = nstart, style = 3)
  }
  
  if(!parallel){
    cluster.ini <- vector("list", nstart)   ##	for containing the best values for 
                                            ##      the parameters after several random starts
    obj.ini <- rep(0, nstart)               ## for containing best objective values
    
    for(j in 1:nstart) {
      assig_obj <- tclust_c1(x, k, alpha, restrC=restrC, deterC=deterC, restr.fact, cshape=cshape, pars,
        niter1, opt, equal.weights, zero_tol)                                                                   # niter1 steps!
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
                            .inorder = F) %dopar% {
      assig_obj <- tclust_c1(x, k, alpha, restrC=restrC, deterC=deterC, restr.fact, cshape=cshape, pars,
        niter1, opt, equal.weights, zero_tol)
      
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
  best_index  <- order(obj.ini, decreasing = TRUE)[1:nkeep]
  best_assig_list <- cluster.ini[best_index]
  
  ###
  # THIRD STEP: apply niter2 concentration steps to nkeep best solutions, return the best solution
  ###
  if(trace){
    cat(paste("\nPhase 3: applying ", niter2, " concentration steps to each of the ", nkeep," best solutions.\n", sep = ""))
    pb2 <- txtProgressBar(min = 0, max = nkeep, style = 3)
  }
  
    best_iter <- NULL
    best_iter_obj <- -Inf
    
    for(j in 1:nkeep){
        iter <- tclust_c2(x, k, best_assig_list[[j]], alpha, restrC=restrC, deterC=deterC, restr.fact, cshape=cshape, pars,
            niter2, opt, equal.weights, zero_tol=1e-16)
        
        if(iter$obj > best_iter_obj){
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

    ## Adjust the returned object to be similar to 'tclust': 
    
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
 
    best_iter$cov <- best_iter$cov[, , idxuse, drop = FALSE]
   
    best_iter$cluster <- as.vector(best_iter$cluster)
    best_iter$cluster <- ClusterIDs[best_iter$cluster + 1]
    
    best_iter$size <- as.vector(best_iter$size)
    best_iter$size <- best_iter$size[idxuse]
    best_iter$weights <- as.vector(best_iter$weights)
    best_iter$weights <- best_iter$weights[idxuse]
     
    best_iter$disttom <- NULL       # remove disttom which is used only in tkmeans objects
    
    ret <- c(best_iter, list(cluster.ini=matrix(unlist(cluster.ini), byrow=TRUE, nrow=length(cluster.ini)),
                obj.ini=obj.ini, int=int, par=parlist, k=sum(best_iter$size > 0)))

    ##  - Dimnames of center and cov
    dn.x <- dimnames(x)
    rownames(ret$centers) <- if(is.null(colnames(x))) paste ("X", 1:ncol(x)) else colnames (x)
    colnames(ret$centers) <- paste("C", 1:k.real)
    dimnames(ret$cov) <- dimnames(ret$centers)[c (1, 1, 2)]

	## calculate the "unrestr.fact"
    get.Mm <- if(parlist$deter.C) .get.Mm.det else .get.Mm.eigen
    EV <- sapply(1:ret$k, get.Mm, ret=ret, x=x)
    if(all(is.na(EV)))  ret$unrestr.fact <- 1
    else	            ret$unrestr.fact <- ceiling(max(EV, na.rm=TRUE)/min(EV, na.rm=TRUE))

    ## Back transform centers, cov and x
    cmm <- scaled$scale %*% t (scaled$scale)
    for(i in 1:ret$k) {
        ret$centers[,i] <- ret$centers[,i] * scaled$scale + scaled$center
        ret$cov[,,i] <- ret$cov[,,i] * cmm
    }
    x <- sweep(x, 2L, scaled$scale, `*`, check.margin = FALSE)
    x <- sweep(x, 2L, scaled$center, `+`, check.margin = FALSE)
    if(store_x)
        ret$par$x <- x
    
    ## Calculate mahalanobis distances
    if(trace)
        cat("\nCalculate mahalanobis distances...")
    
    ret$mah <- array (dim = nrow(x))
    ret$mah[!ret$cluster] <- NA

    for(i in 1:ret$k) {
        idx <- ret$cluster == i
        if(trace){
            cat("\ni=", i, "idx=", idx, "\n")
            cat("\ni=", i, "size=", ret$size, "\n")
            print(ret$cov[,, i])
        }
        if(is_singular(ret$cov[,, i])) {
            warning(paste0("Class ", i, " has a singular covariance matrix and size=", ret$size[i], "."))
            ret$mah[idx] <- NA
        } else
            ret$mah[idx] <- mahalanobis(x[idx, , drop=FALSE], center=ret$centers[, i], cov=ret$cov[,, i])
    }

    ## VT::25.09.2024
    ## Compute information criteria
    ic <- .getIC(ret)
    if(opt == "HARD") 
        ret$CLACLA <- ic$CLACLA
    else {
        ret$MIXMIX <- ic$MIXMIX
        ret$MIXCLA <- ic$MIXCLA
    }
    
    class(ret) <- "tclust"    
    return(ret)
}

.get.Mm.eigen <- function(ii, ret, x) {
	idx <- which(ii==ret$cluster)
	if(length(idx) <= 1)
		return(rep(NA, ncol(x)))
	eigen(cov(x[idx, , drop=FALSE]))$value
}

.get.Mm.det <- function (ii, ret, x) {
	idx <- which(ii==ret$cluster)
	if(length(idx) <= 1)
		return(NA)
	det(cov(x[idx, , drop=FALSE]))
}

.getIC <- function(obj) {

    n <- obj$par$n
    p <- obj$par$p
    MIXMIX <- MIXCLA <- CLACLA <- NULL
    h <- floor((1 - obj$par$alpha) * n)
    
    NlogL <- obj$NlogL
    NlogLmixt <- if(obj$par$opt == "HARD") NULL else -2*obj$obj
    
    npar <- p * obj$k            # p * k
    if(!obj$par$equal.weights)   # if equalweights = false the k-1 mixture proportions parameters must be added
        npar = npar + (obj$k-1)
    
    if(obj$par$restr == "eigen")
        nParam <- npar + 0.5 * p * (p-1) * obj$k + (p * obj$k - 1) * (1 - 1/obj$par$restr.fact) +1
    else if(obj$par$restr == "deter")
        nParam <-   npar + 0.5 * p * (p-1) * obj$k  + (obj$k - 1) * (1 - 1/(obj$par$restr.fact^(1/p))) + 1 + obj$k * (p-1) * (1 - 1/obj$par$cshape)
    else if(obj$par$restr == "GPCM")
        nParam <- .getIC_GPCM(obj)
         
    logh <- log(h)
    if(obj$par$opt == "HARD")
        CLACLA <- NlogL + nParam * logh
    else {
        MIXCLA <- NlogL + nParam * logh
        MIXMIX <- NlogLmixt + nParam * logh
    }

    ##  cat("\nNlogL=", NlogL, " nParam=", nParam, " h=", h, " logh=", logh, "CLALA=", CLACLA, "\n")
    
    list(nParam=nParam, CLACLA=CLACLA, MIXMIX=MIXMIX, MIXCLA=MIXCLA)
}

.getIC_GPCM <- function(obj) {

    # 1. Extract the model type string
    modeltype <- obj$par$pars$pars
    p <- obj$par$p
    k <- obj$par$k
    
    ## Number of estimated parameters
    ## k centroids of size v
    ## 0.5*p*(p+1) estimates for each of the k covariance matrices
    npar <- p * k                   # p * k
    if(!obj$par$equal.weights)      # if equalweights = false the k-1 mixture proportions parameters must be added
        npar = npar + (k-1)

    # 2. Use switch to assign the structural parameters based on modeltype
    params <- switch(modeltype,
      "EII" = list(detpar = 1, shapepar = 0, rotpar = 0),
     
      "VII" = list(detpar = (k - 1) * (1 - 1 / (obj$par$pars$cdet^(1 / p))) + 1, shapepar = 0, rotpar = 0),
      
      "EEI" = list(detpar = 1, shapepar = p - 1, rotpar = 0),
      
      "VEI" = list(detpar = (k - 1) * (1 - 1 / (obj$par$pars$cdet^(1 / p))) + 1, shapepar = p - 1, rotpar = 0),
      
      "EVI" = list(detpar = 1, rotpar = 0, shapepar = (p - 1) * (1 - 1 / obj$par$pars$shw) * ((k - 1) * (1 - 1 / obj$par$pars$shb) + 1)),
      
      "VVI" = list(detpar = (k - 1) * (1 - 1 / (obj$par$pars$cdet^(1 / p))) + 1, shapepar = (p - 1) * (1 - 1 / obj$par$pars$shw) * ((k - 1) * (1 - 1 / obj$par$pars$shb) + 1), rotpar = 0),
      
      "EEE" = list(detpar = 1, shapepar = p - 1, rotpar = 0.5 * p * (p - 1)),
      
      "VEE" = list(detpar = (k - 1) * (1 - 1 / (obj$par$pars$cdet^(1 / p))) + 1, shapepar = p - 1, rotpar = 0.5 * p * (p - 1)),
      
      "EVE" = list(detpar = 1, shapepar = (p - 1) * (1 - 1 / obj$par$pars$shw) * ((k - 1) * (1 - 1 / obj$par$pars$shb) + 1), rotpar = 0.5 * p * (p - 1)),
      
      "EEV" = list(detpar = 1, shapepar = p - 1, rotpar = 0.5 * k * p * (p - 1)),
      
      "VVE" = list(detpar = (k - 1) * (1 - 1 / (obj$par$pars$cdet^(1 / p))) + 1, shapepar = (p - 1) * (1 - 1 / obj$par$pars$shw) * ((k - 1) * (1 - 1 / obj$par$pars$shb) + 1), rotpar = 0.5 * p * (p - 1)),
      
      "VEV" = list(detpar = (k - 1) * (1 - 1 / (obj$par$pars$cdet^(1 / p))) + 1, shapepar = p - 1, rotpar = 0.5 * k * p * (p - 1)),
      
      "EVV" = list(detpar = 1, shapepar = (p - 1) * (1 - 1 / obj$par$pars$shw) * ((k - 1) * (1 - 1 / obj$par$pars$shb) + 1), rotpar = 0.5 * k * p * (p - 1)),
      
      "VVV" = list(detpar = (k - 1) * (1 - 1 / (obj$par$pars$cdet^(1 / p))) + 1, shapepar = (p - 1) * (1 - 1 / obj$par$pars$shw) * ((k - 1) * (1 - 1 / obj$par$pars$shb) + 1), rotpar = 0.5 * k * p * (p - 1)),
      
      # Default fallback case if no text matches match
      stop("FSDA:tclust:WrongModel: Wrong model for cov matrices, must be one of the 14 GPCM")
    )
    
    # 3. Unpack variables into the global environment
    detpar   <- params$detpar
    shapepar <- params$shapepar
    rotpar   <- params$rotpar
    
    # 4. Calculate total parameters
    nParam   <- npar + detpar + shapepar + rotpar    
    
    nParam
} 
