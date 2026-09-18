##  roxygen2::roxygenise("C:/users/valen/onedrive/myrepo/rrdev/robClus", load_code=roxygen2:::load_installed)

#' Robust Linear Grouping
#' 
#' @description The function \code{rlg()} searches for clusters around affine subspaces of dimensions given by 
#'  vector \code{d} (the length of that vector is the number of clusters). For instance \code{d=c(1,2)} 
#'  means that we are clustering around a line and a plane. For robustifying the estimation, 
#'  a proportion \code{alpha} of observations is trimmed. In particular, the trimmed k-means 
#'  method is represented by the rlg method, if \code{d=c(0,0,..0)} (a vector of length 
#'  \code{k} with zeroes).
#'
#' @param x A matrix or data.frame of dimension n x p, containing the observations (rowwise).
#' @param d A numeric vector of length equal to the number of clusters to be detected. 
#'  Each component of vector \code{d} indicates the intrinsic dimension of the affine subspace 
#'  where observations on that cluster are going to be clustered. All the elements 
#'  of vector \code{d} should be smaller than the problem dimension minus 1.
#' @param alpha The proportion of observations to be trimmed. 
#' @param nstart The number of random initializations to be performed.
#' @param niter1 The number of concentration steps to be performed for the nstart initializations.
#' @param nkeep The number of iterated initializations (after niter1 concentration 
#'  steps) with the best values in the target function that are kept for further iterations
#' @param niter2 The maximum number of concentration steps to be performed for 
#'  the nkeep solutions kept for further iteration. The concentration steps 
#'  are stopped, whenever two consecutive steps lead to the same data partition.
#' @param scale A robust centering and scaling (using the median and MAD) is done if TRUE.
#' @param parallel A logical value, specifying whether the nstart initializations should be done in parallel.
#' @param n.cores The number of cores to use when paralellizing, only taken into account if parallel=T.
#' @param trace Defines the tracing level, which is set to 0 by default. Tracing level 1 gives additional information on the stage of the iterative process.
#'
#' @return Returns an object of class \code{rlg} which is basically a list with the following elements:
#' \itemize{
#'     \item centers - A matrix of size p x k containing the location vectors (column-wise) of each cluster. 
#'     \item U - A list with k elements where each element is p x d_j matrix whose d_j columns are unitary and orthogonal vectors generating the affine subspace (after subtracting the corresponding cluster’s location parameter in centers). d_j is the intrinsic dimension of the affine subspace approximation in the j-th cluster, i.e., the elements of vector d.
#'     \item cluster - A numerical vector of size n containing the cluster assignment 
#'          for each observation. Cluster names are integer numbers from 1 to k, 
#'          0 indicates trimmed observations.
#'     \item obj - The value of the objective function of the best (returned) solution.
#'     \item cluster.ini - A matrix with nstart rows and number of columns equal to the number of observations and where each row shows the final clustering assignments (0 for trimmed observations) obtained after the niter1 iteration of the nstart random initializations.
#'     \item obj.ini -A numerical vector of length nstart containing the values of the target function obtained after the niter1 iteration of the nstart random initializations.
#'     \item x - The input data set. 
#'     \item dimensions - The input d vector with the intrinsic dimensions. The number of clusters is the length of that vector. 
#'     \item alpha - The input trimming level.
#' }
#' 
#' @details The procedure allows to deal with robust clustering around affine subspaces 
#'    with an alpha proportion of trimming level by minimizing the trimmed sums of squared 
#'    orthogonal residuals. Each component of vector \code{d} indicates the intrinsic dimension of 
#'    the affine subspace where observations on that cluster are going to be clustered. 
#'    Therefore a component equal to 0 on that vector implies clustering around centres, 
#'    equal to 1 around lines, equal to 2 around planes and so on. The procedure so 
#'    allows simultaneous clustering and dimensionality reduction. 
#'
#' This iterative algorithm performs "concentration steps" to improve the current 
#'  cluster assignments. For approximately obtaining the global optimum, the procedure 
#'  is randomly initialized \code{nstart} times and \code{niter1} concentration steps are performed 
#'  for them. The \code{nkeep} most “promising” iterations, i.e. the \code{nkeep} iterated solutions 
#'  with the initial best values for the target function, are then iterated until 
#'  convergence or until \code{niter2} concentration steps are done.
#' 
#' @author Javier Crespo Guerrero, Jesús Fernández Iglesias, Luis Angel Garcia Escudero, Agustin Mayo Iscar.
#' 
#' @references
#' 
#' García‐Escudero, L. A., Gordaliza, A., San Martin, R., Van Aelst, S., & Zamar, R. (2009). 
#'  Robust linear clustering. Journal of the Royal Statistical Society: 
#'  Series B (Statistical Methodology), 71, 301-318. 
#' 
#' @export
#'
#' @examples
#' ##--- EXAMPLE 1 ------------------------------------------
#' data (LG5data)
#' x <- LG5data[, 1:10]
#' clus <- rlg(x, d = c(2,2,2), alpha=0.1)
#' plot(x, col=clus$cluster+1)
#' plot(clus, which="eigenvalues") 
#' plot(clus, which="scores") 
#'
#' ##--- EXAMPLE 2 ------------------------------------------
#'  data (pine) 
#'  clus <- rlg(pine, d = c(1,1,1), alpha=0.035)
#'  plot(pine, col=clus$cluster+1)
#'  


rlg <- function(x, d, alpha=0.05, nstart=500, niter1=3, niter2=20, nkeep=5, scale=FALSE, 
    parallel=FALSE, n.cores=-1, trace=FALSE){
  x <- as.matrix(x)
  if(is.null(nstart)){nstart <- sum(d)*40}
  
  if(scale){
    standard <- function(x){(x-median(x,na.rm=TRUE))/mad(x,na.rm=TRUE)}
    x <- apply(x,2,standard)
  }
  
  clusters<-matrix(nrow = nstart, ncol = dim(x)[1])
  ecms<-rep(0,nstart)
  
  if(trace){
    cat(paste("\nPhase 1: obtaining ", nstart, " solutions.\n", sep = ""))
    if(parallel) cat(paste("\n Parallelizing initializations using ", n.cores, 
        " cores. Progress bar will not display accurate information.\n", sep = ""))
    pb <- txtProgressBar(min = 0, max = nstart, style = 3)
  }
  
  if(!parallel){
    for(inicializacion in 1:nstart){
      fit<-rlg_c1(x,d,alpha,niter1=niter1)
      clusters[inicializacion,]<-fit$cluster
      ecms[inicializacion]<-fit$obj
      
      if(trace){
        setTxtProgressBar(pb, inicializacion)
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
                            .multicombine = T,
                            .inorder = F) %dopar% {
                              fit <- rlg_c1(x,d,alpha,niter1=niter1)
                              list(fit$cluster, fit$obj)
                            }
    stopCluster(parclus)
    
    clusters <- init.results[0:(nstart-1) * 2 + 1] # Impair positions of cluster.ini
    clusters <- do.call(rbind, clusters) # Turn list into matrix
    ecms <- unlist(init.results[1:nstart * 2]) # Pair positions of cluster.ini
  }
  
  if(trace){
    cat(paste("\n\nPhase 2: obtaining ", nkeep, " best solutions out of the intial ", 
        nstart," solutions.\n", sep = ""))
  }
  num<-min(nkeep,length(na.omit(ecms)))
  ecms_ord<-sort(ecms,decreasing = FALSE)
  candidatos<-ecms_ord[1:num]
  
  grupos_optimos<-matrix(nrow = 0, ncol = dim(x)[1])
  minecms<-c()
  
  for(i in 1:num){
    for(j in 1:length(ecms)){
      if(!is.na(ecms[j])){
        if(candidatos[i] == ecms[j]){
          if(sum(unique(clusters[j,]) != 0) == length(d)){
            grupos_optimos<-rbind(grupos_optimos,clusters[j,])
            minecms<-c(minecms,ecms[j])
          }
          break
        }
      }
    }
  }
  
  n_optimos<-length(minecms)
  resultado_ecm<-c()
  resultado_grupo<-matrix(nrow = 0, ncol = dim(x)[1])
  
  if(trace){
    cat(paste("\nPhase 3: applying ", niter2, " concentration steps to each of the ", 
        n_optimos," best solutions.\n", sep = ""))
    pb2 <- txtProgressBar(min = 0, max = nkeep, style = 3)
  }
  
  if(n_optimos == 1){
    fit<-rlg_c2(x,d,grupos_optimos,alpha,niter2=niter2)
    
    if(trace){
      setTxtProgressBar(pb2, n_optimos)
    }
    
    grupos <- fit$cluster
    optimo <- fit$obj
    centers <- matrix(NA,nrow=dim(x)[2],ncol=length(d))
    u <- vector("list", length(d))
    
    for (k in 1:length(d)){
      nk <-  dim(x[grupos==k,])[1]
      centers[,k] <- apply(x[grupos==k,],2,mean)
      u[[k]] <- svd(t( x[grupos==k,]-rep(1,nk)%*%t(centers[,k])),nu=d[k])$u
    }
    
    ret <- list(obj = optimo, cluster = grupos , cluster.ini = clusters, obj.ini = ecms , dimensions = d , centers=centers , U = u, x=x, alpha=alpha)
    
    
    return(ret)
  }else{
    for(omega in 1:n_optimos){
      fit<-rlg_c2(x,d,grupos_optimos[omega,],alpha,niter2=niter2)
      resultado_grupo<-rbind(resultado_grupo,fit$cluster)
      resultado_ecm<-c(resultado_ecm,fit$obj)
      
      if(trace){
        setTxtProgressBar(pb2, omega)
      }
    }
  }
  if(trace){
    cat("\n\n")
  }
  optimo<-min(resultado_ecm)
  flag<-which.min(resultado_ecm)
  grupos<-resultado_grupo[flag,]
  
  centers <- matrix(NA,nrow=dim(x)[2],ncol=length(d))
  u <- vector("list", length(d))
  for (k in 1:length(d)){
    nk <-  dim(x[grupos==k,])[1]
    centers[,k] <- apply(x[grupos==k,],2,mean)
    #u[[k]] <- princomp(x[grupos==k,])$loadings[,1:d[k],drop=F]
    u[[k]] <- svd(t( x[grupos==k,]-rep(1,nk)%*%t(centers[,k])),nu=d[k])$u
  }
  ret <- list(obj = optimo, cluster = grupos , cluster.ini = clusters, obj.ini = ecms , dimensions = d , centers=centers , U = u, x=x, alpha=alpha)
  class(ret) <- "rlg"
  
  return(ret)
}
