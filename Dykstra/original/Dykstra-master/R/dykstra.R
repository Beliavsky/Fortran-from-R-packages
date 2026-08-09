dykstra <-
  function(Dmat, dvec, Amat, bvec, meq = 0, factorized = FALSE,
           maxit = NULL, eps = NULL){
    # Solve a Quadratic Programming Problem via Dykstra's Algorithm
    #   Dykstra, Richard L. (1983). An algorithm for restricted 
    #   least squares regression. JASA, 78(384), 837-842.
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: Jan 12, 2018
    
    # finds the vector beta that minimizes
    # - t(dvec) %*% beta + 0.5 * t(beta) %*% Dmat %*% beta
    # subject to  t(Amat) %*% beta >= bvec
    
    # Note 1: the first meq constraints are equality constraints
    # Note 2: if factorized = TRUE, first input is solve(Rmat) 
    #         where Dmat = t(Rmat) %*% Rmat
    
    ### check Dmat and dvec
    Dmat <- as.matrix(Dmat)
    dvec <- as.numeric(dvec)
    pts <- length(dvec)
    if(nrow(Dmat) != pts | ncol(Dmat) != pts) stop("Inputs 'Dmat' and 'dvec' are incompatible.")
    
    ### check Amat
    Amat <- as.matrix(Amat)
    if(nrow(Amat) != pts) stop("Input 'Amat' must satisfy:  nrow(Amat) == length(dvec)")
    ncon <- ncol(Amat)
    
    ### check bvec
    if(missing(bvec)){
      bvec <- rep(0, ncon)
    } else {
      bvec <- as.numeric(bvec)
      if(length(bvec) != ncon) stop("Input 'bvec' must satisfy:  length(bvec) == ncol(Amat)")
    }
    
    ### check meq
    meq <- as.integer(meq[1])
    if(meq < 0 | meq > ncon) stop("Input 'meq' must be between 0 and length(bvec).")
    
    ### check maxit
    if(is.null(maxit)){
      maxit <- 30 * pts
    } else {
      maxit <- as.integer(maxit[1])
      if(maxit < 1) stop("Input 'maxit' must be a positive integer.")
    }
    
    ### check eps
    if(is.null(eps)) {
      eps <- .Machine$double.eps * pts
    } else {
      eps <- as.numeric(eps[1])
      if(eps < 0) stop("Input 'eps' must be a non-negative scalar.")
    }
    tol <- .Machine$double.eps * pts
    
    ### check factorized
    factorized <- as.logical(factorized[1])
    if(!factorized && !isSymmetric(Dmat)){
      stop("Input 'Dmat' must be a symmetric matrix when 'factorized = FALSE'.")
    }
    
    ### check if Dmat is diagonal (and make Rinv)
    uti <- upper.tri(Dmat)
    if(max(abs(Dmat[uti])) <= tol){
      diag.flag <- TRUE
      if(factorized){
        Rinv <- diag(Dmat)
      } else {
        Ddiag <- diag(Dmat)
        if(any(Ddiag < -tol)) stop("Input 'Dmat' must be positive definite (or semidefinite).")
        nvals <- sum(Ddiag > (tol * max(Ddiag)))
        if(nvals < pts) Ddiag <- Ddiag + (tol * max(Ddiag) - min(Ddiag))
        Rinv <- 1 / sqrt(Ddiag)
      }
    } else {
      diag.flag <- FALSE
      if(factorized){
        Rinv <- Dmat
      } else {
        Deig <- eigen(Dmat, symmetric = TRUE)
        if(any(Deig$values < -tol)) stop("Input 'Dmat' must be positive definite (or semidefinite).")
        nvals <- sum(Deig$values > (tol * Deig$values[1]))
        if(nvals < pts) Deig$values <- Deig$values + (tol * Deig$values[1] - Deig$values[pts])
        Rinv <- matrix(0, pts, pts)
        for(i in 1:pts) Rinv[,i] <- Deig$vectors[,i] / sqrt(Deig$values[i])
      }
    } # end if(max(abs(Dmat[lti])) <= tol)
    
    ### transform dvec and Amat
    if(diag.flag){
      gvec <- rep(0, pts)
      for(i in 1:pts){
        gvec[i] <- Rinv[i] * dvec[i]
        Amat[i,] <- Rinv[i] * Amat[i,]
      }
    } else {
      gvec <- as.numeric(crossprod(Rinv, dvec))
      Amat <- crossprod(Rinv, Amat)
    }
    Acss <- colSums(Amat^2)
    
    ### unconstrained solution
    if(diag.flag){
      beta.unconstrained <- as.numeric(dvec / Rinv^2)
    } else {
      beta.unconstrained <- as.numeric(Rinv %*% crossprod(Rinv, dvec))
    }
    
    ### rescale eps
    maxb0 <- max(abs(beta.unconstrained))
    if(maxb0 > 1) eps <- eps * maxb0
    
    ### initializations
    zeros <- rep(0, pts)
    beta.change <- matrix(0, pts, ncon)
    beta.solution <- beta.old <- gvec
    iter <- 0L
    ctol <- eps + 1
    
    ### cyclic projection
    while(ctol > eps && iter < maxit){
      
      # loop through constraints
      for(i in 1:ncon){
        beta.work <- beta.solution - beta.change[,i]
        Ai <- sum(Amat[,i] * beta.work)
        passive <- ifelse(i <= meq, Ai == bvec[i], Ai >= bvec[i])
        if(passive){
          beta.change[,i] <- zeros
          beta.solution <- beta.work
        } else {
          beta.change[,i] <- (bvec[i] - Ai) * Amat[,i] / Acss[i]
          beta.solution <- beta.work + beta.change[,i]
        }
      }
      
      # check for convergence
      ctol <- max(abs(beta.solution - beta.old))
      beta.old <- beta.solution
      iter <- iter + 1L
      
    } # end while(ctol > eps && iter < maxit)
    
    # retransform solution
    if(diag.flag){
      beta.solution <- Rinv * beta.solution
    } else {
      beta.solution <- as.numeric(Rinv %*% beta.solution)
    }
    
    # get results
    converged <- ifelse(ctol <= eps, TRUE, FALSE)
    if(factorized){
      value <- NA
    } else {
      value.Q <- 0.5 * crossprod(beta.solution, Dmat %*% beta.solution)
      value.P <- sum(beta.solution * dvec)
      value <- as.numeric(value.Q - value.P)
    }
    results <- list(solution = beta.solution, value = value,
                    unconstrained = beta.unconstrained,
                    iterations = iter, converged = converged)
    class(results) <- "dykstra"
    return(results)
    
  } # end dykstra