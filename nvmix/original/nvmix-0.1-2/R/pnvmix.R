### pnvmix() ###################################################################

##' @title Evaluate Integrand of pnvmix()
##' @param U (n, d) matrix of uniforms (evaluation points)
##' @param qmix see ?pnvmix
##' @param rmix see ?pnvmix 
##' @param lower see ?pnvmix
##' @param upper see ?pnvmix
##' @param scale see ?pnvmix. Has to be full rank here though!
##' @param precond see ?get_set_param()
##' @param mean.sqrt.mix see ?get_set_param()
##' @param return.all logical if all function evaluations should be returned.
##' @param ... additional parameters passed to qmix()
##' @return if return.all = TRUE a n-vector with g(U) values, otherwise a 2-vector
##' with estimated mean(g(U)) and estimated var(g(U))
##' @note if return.all = TRUE, all computations are done in R, so can be slow.
##'       if return.all = FALSE, this function calls underlying C code.
##' @author Erik Hintz
##' @note - This function is *only* needed for numerical experiments.
##'       - It corresponds to 'g' in the paper.
pnvmix_g <- function(U, qmix, rmix = NULL, upper, lower = rep(-Inf, d), 
                     groupings = rep(1, d), scale, precond, mean.sqrt.mix = NULL, 
                     return.all = FALSE, verbose = TRUE, do.ant = TRUE, ...)
{
   d <- dim(scale)[1]
   ## Define the quantile function of the mixing variable
   mix_list <- 
      get_mix_(qmix = qmix, groupings = groupings, callingfun = "pnvmix", ...) # function(u)
   mix_     <- mix_list[[1]]
   use.q    <- mix_list$use.q 
   ## Check if 'do.ant' and 'use.q' are compatible
   if(do.ant & !use.q){
      warning("Antithetic variates only available when 'qmix' provided")
      do.ant <- FALSE
   }
   ## Handle 'mean.sqrt.mix' 
   if(precond & d > 2){
      if(is.null(mean.sqrt.mix)){
         ## 'mean.sqrt.mix' not provided => check if it's in mix_list
         if(!is.null(mix_list$mean.sqrt.mix)){
            mean.sqrt.mix <- mix_list$mean.sqrt.mix
         } else {
            ## Approximate it 
            mean.sqrt.mix <- if(use.q) 
               mean(sqrt(mix_(qrng::sobol(n = 2^9, d = 1, randomize = TRUE)))) else
                  mean(sqrt(mix_(2^9)))
         }
      } else { # 'mean.sqrt.mix' was provided 
         stopifnot(length(mean.sqrt.mix) == length(unique(groupings)))
      } 
      mean.sqrt.mix <- mean.sqrt.mix[groupings]   
      ## Check if provided/approximated 'mean.sqrt.mix' is strictly positive
      if(any(mean.sqrt.mix <= 0))
         stop("'mean.sqrt.mix' has to be positive (possibly after being generated in pnvmix())")
   }
   if(!is.matrix(U)) U <- as.matrix(U)
   ## Generate realizations of sqrt(W) 
   rtW <- if(use.q) sqrt(mix_(U[, 1])) else sqrt(mix_(dim(U)[1]))
   rtWant <- if(do.ant) sqrt(mix_(1 - U[, 1])) else -1 
   ## Call pnvmix_g_() to do the calculation
   pnvmix_g_(U = U[, 2:d, drop = FALSE], rtW = rtW, rtWant = rtWant, 
             groupings = groupings, upper = upper, 
             lower = lower, scale = scale, precond = precond, 
             mean.sqrt.mix = mean.sqrt.mix, return.all = return.all, 
             verbose = verbose, do.ant = do.ant)
}

##' @title Evaluate Integrand of pnvmix()
##' @param U (n, d-1) matrix of uniforms (used to generate normals)
##' @param rtW (n, numgroups) matrix of realizations of sqrt(W)
##' @param rtWant antithetic realizations of rtW (or -1 if do.ant = F)
##' @param groupings see ?pnvmix 
##' @param lower see ?pnvmix
##' @param upper see ?pnvmix
##' @param scale see ?pnvmix. *must* be full rank here. 
##' @param precond see ?get_set_param()
##' @param mean.sqrt.mix see ?get_set_param()
##' @param return.all logical if all function evaluations should be returned.
##' @param verbose logical if wanrnings shall be returned
##' @param do.ant logical of antithetic variates are to be used 
##' @return if return.all = TRUE a n-vector with g(U) values, otherwise a 2-vector
##' with estimated mean(g(U)) and estimated var(g(U))
##' @note if return.all = TRUE, all computations are done in R, so can be slow.
##'       if return.all = FALSE, this function calls underlying C code.
##' @author Erik Hintz
##' @note - This function is *only* needed for numerical experiments and is not
##'         called by any other function in the package 'nvmix' 
pnvmix_g_ <- function(U, rtW, rtWant, groupings = rep(1, d), upper, 
                      lower = rep(-Inf, d), scale, precond, mean.sqrt.mix, 
                      return.all = FALSE, verbose = TRUE, do.ant = TRUE)
{
   if(!is.matrix(U)) U <- as.matrix(U)
   if(!is.matrix(rtW)) rtW <- as.matrix(rtW)
   if(!is.matrix(rtWant)) rtWant <- as.matrix(rtWant) 
   ## Dimension of the problem and number of evaluations
   d <- dim(scale)[1]
   n <- dim(U)[1]
   stopifnot(all.equal(dim(rtW), c(n, length(unique(groupings)))))
   if(do.ant)  stopifnot(all.equal(dim(rtWant), c(n, length(unique(groupings)))))
   ## Factor (lower triangular)
   factor   <- t(chol(scale))
   rank     <- d # only consider full rank case here
   k.factor <- rep(1, d)
   ## Precondition?
   if(precond & d > 2) {
      ## 'mean.sqrt.mix' must have length d and contain E(sqrt(W_j)), j=1,..,d
      temp <- precondition(lower, upper = upper, scale = scale, factor = factor,
                           mean.sqrt.mix = mean.sqrt.mix)
      if(is.null(temp)) {
         ## Preconditioning did not work, continue with original inputs
         if(verbose) warning("Preconditioning led to (numerically) singular 'scale',
                             continuing with original input.")
      } else {
         lower  <- temp$lower
         upper  <- temp$upper
         factor <- temp$factor # *vector* containing lower triangular part 
         groupings <- groupings[temp$perm] 
      }
   } else { # lower triangular part of 'factor' as vector 
      factor <- factor[lower.tri(factor, diag = TRUE)] 
   } 
   ## For evaluating qnorm() close to 0 and 1
   ONE <- 1-.Machine$double.neg.eps
   ZERO <- .Machine$double.eps
   
   if(return.all) {
      ## Here (and only here) do we need 'factor' as a matrix (readability)
      factor_ <- matrix(0, ncol = d, nrow = d)
      factor_[lower.tri(factor_, diag = TRUE)] <- factor
      factor <- factor_
      ## Matrix to store results (y_i from paper)
      Yorg <- matrix(NA, ncol = d - 1, nrow = dim(U)[1])
      ## First 'iteration' (d1, e1 in the paper)
      dorg <- pnorm(lower[1] / (rtW[, groupings[1]] * factor[1, 1]))
      eorg <- pnorm(upper[1] / (rtW[, groupings[1]] * factor[1, 1]))
      forg <- eorg - dorg
      if(do.ant){
         ## Antithetic values
         Yant <- matrix(NA, ncol = d - 1, nrow = dim(U)[1])
         dant <- pnorm(lower[1] / (rtWant[, groupings[1]] * factor[1, 1]))
         eant <- pnorm(upper[1] / (rtWant[, groupings[1]] * factor[1, 1]))
         fant <- eant - dant
      }
      ## Recursively calculate (e_i - d_i)
      for(i in 2:d) {
         ## Store realization:
         Yorg[,(i-1)] <- qnorm(pmax(pmin(dorg + U[,i-1]*(eorg-dorg), ONE), ZERO))
         ## Update d__, e__, f___:
         dorg <- pnorm((lower[i]/rtW[, groupings[i]] - Yorg[,1: (i-1)] %*% 
                           as.matrix(factor[i, 1:(i-1)]))/factor[i,i])
         eorg <- pnorm((upper[i]/rtW[, groupings[i]] - Yorg[,1: (i-1)] %*% 
                           as.matrix(factor[i, 1:(i-1)]))/factor[i,i])
         forg <- forg*(eorg-dorg)
         if(do.ant){
            ## Antithetic values
            Yant[, (i-1)] <- 
               qnorm( pmax( pmin(dant + (1-U[, i-1])*(eant-dant), ONE), ZERO))
            dant <- pnorm((lower[i]/rtWant[, groupings[i]] - Yant[,1: (i-1)] %*% 
                              as.matrix(factor[i, 1:(i-1)]))/factor[i,i])
            eant <- pnorm((upper[i]/rtWant[, groupings[i]] - Yant[,1: (i-1)] %*% 
                              as.matrix(factor[i, 1:(i-1)]))/factor[i,i])
            fant <- fant*(eant-dant)
         }
      }
      ## Return:
      if(do.ant) (forg+fant)/2 else forg 
   } else { # return,all = FALSE => call C function (faster) 
      .Call("eval_nvmix_integral",
            lower    = as.double(lower),
            upper    = as.double(upper),
            gropings = as.integer(groupings),
            numgroups= as.integer(length(unique(groupings))),
            U        = as.double(U),
            rtW      = as.double(rtW),
            rtWant   = as.double(rtWant), 
            n        = as.integer(n),
            d        = as.integer(d),
            r        = as.integer(rank),
            kfactor  = as.integer(k.factor),
            factor   = as.double(factor), # lower triangular part as vector 
            ZERO     = as.double(ZERO),
            ONE      = as.double(ONE),
            doant    = as.integer(do.ant))
   } 
}


##' @title Cholesky Factor for Positive-Semidefinite Matrices
##' @param mat (n,n) symmetric, positive-semidefinite matrix.
##' @return list of length 2:
##'         - C: (n,n) lower triangular matrix such that C % * % t(C) = mat
##'         - D: n-vector with diagonal elements of 'C'
##' @author Erik Hintz and Marius Hofert
##' @note Internal function being called by pnvmix().
cholesky_ <- function(mat, tol = 1e-12)
{
   n <- dim(mat)[1] # dimension
   stopifnot(dim(mat)[2] == n)
   ## First try 'chol()' (=>fast)
   C <- tryCatch(t(chol(mat)), error = function(e) e)
   if(is.matrix(C) && all.equal(dim(C), rep(n, 2))) {
      ## C is the desired Cholesky factor
      ## Grab diagonal
      diag.elements <- diag(C)
   } else {
      ## In this case, 't(chol(scale))' returned an error so that we manually
      ## compute the Cholesky factor of the *singular* matrix 'mat'
      C <- matrix(0, ncol = n, nrow = n) # initialize Cholesky factor
      diag.elements <- rep(NA, n)
      for(col in 1:n) {
         dsq <- mat[col, col] - sum(C[col, 1:(col-1)]^2) # C(col,col)^2
         if(dsq < 0) stop("Matrix not positive semi-definite")
         d <- if(dsq < abs(mat[col, col] * tol)) 0 else sqrt(dsq) # set 0 if too small
         C[col, col] <- d
         diag.elements[col] <- d # store diagnonal element
         if(col < n && d > 0) { # calculate the remaining elements in column 'col'
            for(row in (col+1):n) {
               C[row, col] <- (mat[row, col] -
                                  sum(C[row, 1:(row-1)]*C[col, 1:(row-1)]))/d
            }
         }
      }
   }
   list(C = C, D = diag.elements)
}


##' @title Swap Variables i and j in a, b and R
##' @param i variable to be switched with j
##' @param j variable to be switched with i
##' @param lower d-vector of lower evaluation limits
##' @param upper d-vector of upper evaluation limits
##' @param scale (d, d)-covariance matrix (scale matrix)
##' @return list with lower, upper and scale after components/rows/columns
##'         i and j have been switched
##' @author Erik Hintz and Marius Hofert
swap <- function(i, j, lower, upper, scale)
{
   ## Build vector
   ij <- c(i,j)
   ji <- c(j,i)
   ## Reorder lower, upper
   lower[ij] <- lower[ji]
   upper[ij] <- upper[ji]
   ## Reorder scale
   wo.ij <- setdiff(seq_len(nrow(scale)), ij)
   temp_i <- scale[wo.ij,i, drop = FALSE]
   temp_j <- scale[wo.ij,j, drop = FALSE]
   temp_ii <- scale[i,i]
   scale[wo.ij,i] <- temp_j
   scale[wo.ij,j] <- temp_i
   scale[i,wo.ij] <- temp_j
   scale[j,wo.ij] <- temp_i
   scale[i,i] <- scale[j,j]
   scale[j,j] <- temp_ii
   ## Return
   list(lower = lower, upper = upper, scale = scale)
}


##' @title Reorder Limits and Scale Matrix According to a Permutation
##' @param perm d-vector giving the desired permutation
##' @param upper see ?pnvmix
##' @param lower see ?pnvmix
##' @param scale see ?pnvmix
##' @return list with lower, upper and scale after components/rows/columns
##'         have been switched according to perm
##' @author Erik Hintz
##' @note No input checking is done. This function is mostly for experimenting.
reorder_limits_scale <- function(perm, upper, lower = rep(-Inf, d), scale = diag(d))
{
   d <- length(upper)
   ## Vector to save current positions of original variables
   org.perm <- 1:d
   for(i in 1:d) {
      ## Variable at position i before swapping
      curr.i <- org.perm[i]
      ## Variable at position i after swapping
      next.i <- perm[i]
      ## Do we need to swap?
      if(curr.i != next.i) {
         ## Index of the variable to be put to position i
         ind.next.i <- which(org.perm == next.i)
         ## Swap positions i and ind.next.i
         tmp <- swap(i, ind.next.i, lower = lower, upper = upper,
                     scale = scale)
         lower <- tmp$lower
         upper <- tmp$upper
         scale <- tmp$scale
         ## Update position vector
         org.perm[i] <- next.i
         org.perm[ind.next.i] <- curr.i
      }
   }
   list(lower = lower, upper = upper, scale = scale)
}


##' @title Preconditioning (Reordering Variables According to their Expected
##'        Integration Limits)
##' @param lower see ?pnvmix
##' @param upper see ?pnvmix
##' @param scale (d,d) positive definite 'scale' matrix.
##' @param factor Cholesky factor (lower triangular matrix) of 'scale'
##' @param mean.sqrt.mix in NVM case numeric(1) (giving E(sqrt(W))), otherwise
##'        numeric(d) with elements E(sqrt(W_i))
##' @param tol if a calculated diagonal element of factor is < sqrt(tol),
##'        factor is deemed singular and preconditioning fails.
##' @param use.C logical if preconditioning is to be performed in C
##' @param verbose logical if warnings shall be thrown        
##' @return list of length 5 with reordered integration limits, scale matrix and
##'         Cholesky factor as well as a d-vector 'perm' giving the ordering obtained..
##'         Only the lower triangular part of 'scale' and 'factor' are returned. 
##'         If preconditioning was unsuccessful, 'NULL' is returned
##' @author Erik Hintz and Marius Hofert
##' @note See Genz and Bretz (2002, p. 957).
##' It may happen that the original 'scale' admits a (numerically stable full rank)
##' 'factor' and that the reordered 'scale' is numerically singular so that a
##' full rank 'factor' cannot be found. This is detected by 'precondition()'
##' and if this happens, NULL is returned. pnvmix1() then throws a warning
##' and estimation is carried out with un-preconditioned inputs.
precondition <- function(lower, upper, scale, factor, mean.sqrt.mix, 
                         tol = 1e-16, use.C = TRUE, verbose = FALSE){
   
   d <- length(lower)
   ## If scalar 'mean.sqrt.mix' provided => repeat to vector of length d for
   ## compatibility with the gNVM case 
   if(length(mean.sqrt.mix) == 1) mean.sqrt.mix <- rep(mean.sqrt.mix, d) 
   stopifnot(all(mean.sqrt.mix > 0), length(mean.sqrt.mix) == d) # check
   if(!use.C){
      ## Use pure R to perform preconditioning 
      y <- rep(0, d - 1) # to store conditional expected values 
      perm <- 1:d # to record the final permuation used 
      ## Main
      for(j in 1:(d-1)) {
         if(j == 1) {
            denom <- sqrt(diag(scale))
            c <- 0
         } else {
            denom <- diag(scale)[j:d] - .rowSums(factor[j:d, 1:(j-1), drop = FALSE]^2,
                                                 m = d-j+1, n = j-1)
            if(all(denom > 0)) denom <- sqrt(denom) else {
               if(verbose) warning("Computation failed due to singularity; returning NULL.")
               return(NULL)
            }
            c <- factor[j:d, 1:(j-1), drop = FALSE] %*% y[1:(j-1)]
         }
         ## Transformed limits with indices greater than j
         next.uppers <- (upper[j:d] / mean.sqrt.mix[j:d] - c) / denom
         next.lowers <- (lower[j:d] / mean.sqrt.mix[j:d] - c) / denom
         ## Find i = argmin { <expected length of interval j> }
         i <- tail(which.min(pnorm(next.uppers) - pnorm(next.lowers))) + j - 1
         ## Swap variables i and j if they are different
         if(i != j) {
            tmp   <- swap(i = i, j = j, lower = lower, upper = upper, scale = scale)
            lower <- tmp$lower
            upper <- tmp$upper
            scale <- tmp$scale
            perm[c(i, j)] <- perm[c(j, i)]
            mean.sqrt.mix[c(i, j)] <- mean.sqrt.mix[c(j, i)]
            ## If j>1 and an actual swap has occured, need to reorder Cholesky factor:
            if(j > 1) {
               factor[c(i,j),]   <- factor[c(j,i),, drop = FALSE]
               factor[j,(j+1):i] <- matrix(0, ncol = i - j, nrow = 1)
            }
         }
         ## Update Cholesky factor
         if(j == 1) {
            factor[1, 1] <- sqrt(scale[1, 1])
            factor[2:d, 1] <- scale[2:d, 1, drop = FALSE] / factor[1, 1]
            ## Store y1
            low.j.up.j <- c(lower[1], upper[1]) / (mean.sqrt.mix[1]*factor[1, 1])
            y[1] <- (dnorm(low.j.up.j[1]) - dnorm(low.j.up.j[2])) /
               (max(pnorm(low.j.up.j[2]) - pnorm(low.j.up.j[1]), tol)) 
         } else {
            factorjjsq <- scale[j,j] - sum(factor[j,1:(j-1)]^2)
            if(all(factorjjsq > tol)) factor[j,j] <- sqrt(factorjjsq) else {
               warning("Factorjjsq < 0, NULL returned.")
               return(NULL)
            }
            factor[(j+1):d, j] <-
               if(j < d-1) {
                  (scale[(j+1):d, j] - factor[(j+1):d, 1:(j-1), drop = FALSE] %*%
                      factor[j, 1:(j-1)]) / factor[j, j]
               } else {
                  (scale[(j+1):d, j] - factor[(j+1):d, 1:(j-1)] %*%
                      factor[j, 1:(j-1)]) / factor[j, j]
               }
            ## Get yj
            scprod     <- sum(factor[j, 1:(j-1)] * y[1:(j-1)]) # needed twice
            low.j.up.j <- c(lower[j] / mean.sqrt.mix[j] - scprod,
                            upper[j] / mean.sqrt.mix[j] - scprod) / factor[j, j]
            y[j] <- (dnorm(low.j.up.j[1]) - dnorm(low.j.up.j[2])) /
               (max(pnorm(low.j.up.j[2]) - pnorm(low.j.up.j[1]), tol))
         }
      } # for()
      factorddsq <- scale[d, d] - sum(factor[d, 1:(d-1)]^2)
      if(factorddsq > tol) {
         factor[d,d] <- sqrt(factorddsq)
      } else {
         warning("Factorddsq <= tol => Trying chol()")
         ## Try 'chol()' for reordered 'scale' (more accurate)
         factor <- tryCatch(t(chol(scale)), error = function(e) e)
         if(!is.matrix(factor)){
            if(verbose) warning("chol() failed; returning NULL.")
            return(NULL)
         } # else 'factor' is correct => return
      }
      ## Return
      return(list(lower = lower, upper = upper, 
                  scale = scale[lower.tri(scale, diag = TRUE)], 
                  factor = factor[lower.tri(factor, diag = TRUE)], 
                  perm = perm))
   } else {
      ## Use C function "precond" to perform preconditioning
      ## Grab lower triangular part of 'scale' (length d*(d+1)/2)
      scale.tri <- scale[lower.tri(scale, diag = TRUE)]  
      precond_C <- .C("precond", as.double(lower), as.double(upper), 
                      as.double(scale.tri), as.double(rep(0, d*(d+1)/2)), 
                      as.double(mean.sqrt.mix), as.double(1e-16),
                      as.integer(d), as.integer(1:d), as.integer(1),
                      NAOK = TRUE) # if 'lower'/'upper' contain +/-Inf
      ## Arguments = return of 'void precond()': 
      ## [[1]]: lower; [[2]]: upper; [[3]]: scale.tri (vector);
      ## [[4]]: factor.tri (vector); [[5]]: mean.sqrt.mix; [[6]]: tol; [[7]]: d;
      ## [[8]]: perm; [[9]]: status (1: ok; 2: recompute chol in R, >= 10: error)
      if(precond_C[[9]] == 1){ # preconditioning terminated; return
         return(list(lower = precond_C[[1]], upper = precond_C[[2]], 
                     scale = precond_C[[3]], factor = precond_C[[4]], 
                     perm = precond_C[[8]]))
      } else if(precond_C[[9]] > 2){
         if(verbose) warning("Computation failed due to singularity; returning NULL.")
         return(NULL) # error 
      } else {
         ## Try recomputing cholesky factor 
         ## Compute 'scale' from triangular part
         scale <- matrix(NA, ncol = d, nrow = d)
         scale[lower.tri(scale, diag = TRUE)] <- precond_C[[3]]
         scale[upper.tri(scale, diag = TRUE)] <- 
            t(scale)[upper.tri(scale, diag = TRUE)]
         factor <- tryCatch(t(chol(scale)), error = function(e) e)
         if(!is.matrix(factor)){
            if(verbose) warning("chol() failed; returning NULL.")
            return(NULL) # error 
         } # else 'factor' is correct now => return
         return(list(lower = precond_C[[1]], upper = precond_C[[2]], 
                     scale = precond_C[[3]], 
                     factor = factor[lower.tri(factor, diag = TRUE)], 
                     perm = precond_C[[8]]))
      }
   }
}


##' @title Distribution Function of a Multivariate Normal Variance Mixture
##'        for a Single Observation
##' @param upper d vector
##' @param lower d vector (<= upper)
##' @param mix_ quantile function or RNG of the mixture distribution; built inside pnvmix();
##' @param use.q logical if 'mix_' is a quantile function 
##' @param is.const.mix logical, TRUE if qmix is constant (=> normal distutions)
##' @param mean.sqrt.mix see details in ?pnvmix
##' @param loc see details in ?pnvmix
##' @param scale see details in ?pnvmix
##' @param factor lower triangular part of 'factor' as a vector 
##' @param k.factor vector of length rank(scale) with heights of each step in
##'        'factor'
##' @param method see details in ?pnvmix
##' @param precond see details in ?pnvmix
##' @param tol absolute/relative error tolerance depending on 'do.reltol'
##' @param do.reltol logical if relative precision is required
##' @param CI.factor see details in ?pnvmix
##' @param fun.eval see details in ?pnvmix
##' @param increment see details in ?pnvmix
##' @param B see details in ?pnvmix
##' @param verbose see ?pnvmix
##' @param seeds  B - vector with integer seeds for 'sobol()'    
##' @param ... see details in ?pnvmix
##' @return list of length 5:
##'         - value: computed probability
##'         - error: error estimate (depending on do.reltol)
##'         - abserror: absolute error estimate 
##'         - relerror: relative error estimate
##'         - numiter: number of iterations needed
##' @author Erik Hintz and Marius Hofert
##' @note Internal function being called by pnvmix; all error estimates
##'       returned for backwards compatibility 
pnvmix1 <- function(upper, lower = rep(-Inf, d), groupings = rep(1, d), 
                    mix_ = NULL, use.q, do.ant,
                    is.const.mix = FALSE, 
                    mean.sqrt.mix, loc = rep(0, d), scale = diag(d), 
                    factor = NULL, k.factor = rep(1, d),
                    method = c("sobol", "ghalton", "PRNG"), precond = TRUE,
                    tol = 1e-3, do.reltol = FALSE, CI.factor,
                    fun.eval, max.iter.rqmc, 
                    increment = c("doubling", "num.init"),
                    B = 15,  verbose = TRUE, seeds = NULL, mix.stored = NULL,
                    mix.doStore = FALSE, maxiter.stored = 4, ...)
{
   numgroups <- length(unique(groupings))
   rank <- length(k.factor)
   d    <- sum(k.factor)
   stopifnot(length(lower) == d, lower <= upper, is.function(mix_))
   if(any(lower == upper))
      return(list(value = 0, error = 0, numiter = 0))

   ## Preconditioning (resorting the limits; only for d > 2 and non-singular case)
   if(precond & d > 2 & rank == d) {
      ## Note that 'mean.sqrt.mix' has already been calculated in pnvmix()
      temp <- precondition(lower = lower, upper = upper, scale = scale,
                           factor = factor, mean.sqrt.mix = mean.sqrt.mix)
      if(is.null(temp)) {
         ## Preconditioning did not work, continue with original inputs
         if(verbose)
            warning("Preconditioning led to (numerically) singular 'scale', continuing with original input.")
      } else {
         lower  <- temp$lower
         upper  <- temp$upper
         factor <- temp$factor # lower triangular part only
         groupings <- groupings[temp$perm] 
      }
   } else if(is.null(factor)){ 
      ## 'factor' was not provided (happens when pnvmix1() is *not* called from pgnmvix())
      factor <- t(chol(scale))
      factor <- factor[lower.tri(factor, diag = TRUE)] # lower triangular part only
   }

   ## 1 Basics for while loop below ###########################################
   
   ## Error is calculated as CI.factor * sd( estimates) / sqrt(B)
   ## For performance:
   CI.factor.sqrt.B <- CI.factor / sqrt(B)
   ## Grab the number of sample points for the first iteration
   current.n <- fun.eval[1] 
   ## Vector to store 'B' RQMC estimates
   rqmc.estimates <- rep(0, B)
   ## Initialize error to > 'tol' so that we can enter the while loop below
   error <- tol + 42
   ## Initialize the total number of function evals and number of iterations
   total.fun.evals <- 0
   numiter <- 0
   ## It may happen that qnorm(u) for u too close to 1 (or 0) is evaluated;
   ## in those cases, u will be replaced by ONE and ZERO which is the largest
   ## (smallest) number different from 1 (0) such that qnorm(u) is not +/- Inf
   ZERO <- .Machine$double.xmin
   ONE  <- 1 - .Machine$double.neg.eps
   ## Sample 'B' seeds for 'sobol(..., seed = seeds[b])' to get the same shifts 
   if(method == "sobol"){
      if(is.null(seeds)) seeds <- sample(1:(1e3*B), B) else stopifnot(length(seeds) == B)
   } else {
      mix.doStore <- FALSE # only works for sobol, because of the same 'seeds'
   }
   ## Additional variables needed if the increment chosen is "doubling"
   do.doubling <- (increment == "doubling")
   if(do.doubling){
      if(method == "sobol") useskip <- 0
      denom <- 1
   }
   ## Deal with 'mix.stored' and 'mix.doStore'
   if(mix.doStore){
      ## No stored mixings provided but storing required => create 'mix.stored'
      if(is.null(mix.stored)){
         max.nrow <- if(do.doubling) current.n* B *2^(maxiter.stored-1) else 
            maxiter.stored * B *current.n # maximum number of mix_ evaluations possible
         num.col <- if(do.ant) 2 * numgroups else numgroups 
         ## 'mix.stored' = list of length B, b'th element = matrix containing the
         ## mixings for the b'th seed (Col's 1-numgroups: original; 
         ## col's numgroups + 1 - 2*numgroups: antithetic, if applicable)
         mix.stored <- lapply(1:B, function(i) matrix(, ncol = num.col, nrow = max.nrow))
         numiter.stored <- -1 # number of iterations worth of mixings stored - 1
      } else {
         ## Stored mixings provided => grab how many 
         numiter.stored <- attr(mix.stored, "numiter.stored") 
      }
      ## Function returning the row indices in iteration 'numiter' 
      whichRows <- function(numiter){
         ind <- if(do.doubling){
            if(numiter == 0) 1:fun.eval[1] else 
               2^(numiter-1) * fun.eval[1] + (1 : (2^(numiter-1)*fun.eval[1]))
         } else numiter * fun.eval[1] + (1 : fun.eval[1])
         ind # return
      }
   } else {
      numiter.stored <- -1
   }
   
   ## 2 Major while() loop ####################################################
   
   ## while() runs until precision 'tol' is reached or the number of function
   ## evaluations exceed fun.eval[2] or until 'max.iter.rqmc' is exhausted.
   ## In each iteration, B RQMC estimates of the desired probability are calculated.
   while(error > tol & total.fun.evals < fun.eval[2] & numiter < max.iter.rqmc)
   {
      ## Get B RQCM estimates
      for(b in 1:B)
      {
         ## 2.1 Get the point set ###########################################
         ## Construct 'rtW' (realizations of sqrt(W)), 'rtWant' (antithetic realizations;
         ## only if(do.ant)), 'U' (uniforms being inverted to normals later)
         rtWant <- -1 # overwritten if(do.ant); o.w. any numeric value will do
         ## If is.const.mix = TRUE, we only need (rank - 1) (quasi) random numbers
         ## ('is.const.mix = TRUE' and 'rank = 1' has already been dealt with)
         U <- if(is.const.mix) {
            U <- switch(method, 
                        "sobol" = {
                           if(do.doubling) {
                              qrng::sobol(n = current.n, d = rank - 1,
                                          randomize = "digital.shift",
                                          seed = seeds[b],
                                          skip = (useskip * current.n))
                           } else {
                              qrng::sobol(n = current.n, d = rank - 1,
                                          randomize = "digital.shift",
                                          seed = seeds[b],
                                          skip = (numiter * current.n))
                           }
                        },
                        "ghalton" = {
                           qrng::ghalton(n = current.n, d = rank - 1,
                                         method = "generalized")
                        },
                        "PRNG" = {
                           matrix(runif( current.n * (rank - 1)), ncol = rank - 1)
                        })
            ## Realizations of sqrt(W) are all 1 
            rtW    <- rep(1, current.n)
            rtWant <- rep(1, current.n)
            U
         } else {
            U <- switch(method,
                        "sobol" = {
                           if(increment == "doubling") {
                              qrng::sobol(n = current.n, d = rank,
                                          randomize = "digital.shift",
                                          seed = seeds[b],
                                          skip = (useskip * current.n))
                           } else {
                              qrng::sobol(n = current.n, d = rank,
                                          randomize = "digital.shift",
                                          seed = seeds[b],
                                          skip = (numiter * current.n))
                           }
                        },
                        "ghalton" = {
                           qrng::ghalton(n = current.n, d = rank,
                                         method = "generalized")
                        },
                        "PRNG" = {
                           matrix(runif( current.n * rank), ncol = rank)
                        })
            if(!is.matrix(U)) U <- cbind(U)
            ## Compute/grab mixing realizations 
            ## If d == 1, U <- -1 as not needed anymore (no normals, only W)
            if(use.q){ # use quantile function (as opposed to RNG)
               if(!mix.doStore | 
                  (mix.doStore & (numiter > min(numiter.stored, 
                                                maxiter.stored)))){
                  ## No storing, or storing but not stored yet, or storing
                  ## but beyond numiters to store: Call 'mix_' 
                  rtW <- sqrt(mix_(U[, 1]))
                  if(do.ant) rtWant <- sqrt(mix_(1 - U[, 1]))
                  ## Store if needed 
                  if(mix.doStore & numiter < maxiter.stored){
                     if(!is.matrix(rtW)) rtW <- cbind(rtW)
                     if(!is.matrix(rtWant)) rtWant <- cbind(rtWant)
                     mix.stored[[b]][whichRows(numiter), 1:numgroups]<- rtW
                     if(do.ant) mix.stored[[b]][
                        whichRows(numiter), numgroups+(1:numgroups)] <- rtWant
                     if(b == B) numiter.stored <- numiter.stored + 1 
                  }
               } else {
                  ## Grab realizations from 'mix.stored'
                  rtW <- 
                     mix.stored[[b]][whichRows(numiter), 1:numgroups, drop = FALSE]
                  if(do.ant) rtWant <- mix.stored[[b]][
                     whichRows(numiter), numgroups+(1:numgroups), drop = FALSE]
               }
               ## Return (for 'U') correctly
               if(d == 1) -1 else U[, 2:rank] 
            } else {
               ## If use.q = FALSE, do.ant = FALSE automatically 
               rtW <- sqrt(mix_(current.n))
               if(d == 1) -1 else U[, 2:rank] 
            }
         }
         if(!is.matrix(rtW)) rtW <- cbind(rtW)
         if(!is.matrix(rtWant)) rtWant <- cbind(rtWant)
         
         ## 2.2 Evaluate the integrand at the (next) point set ##############
         
         next.estimate <-
            if(d == 1) {
               ## If d=1, use pnorm(); normal + t in d = 1 already addressed 
               if(do.ant){
                  mean( (pnorm(upper/rtW[, groupings]) - pnorm(lower/rtW[, groupings]) + 
                            pnorm(upper/rtWant[, groupings]) - pnorm(lower/rtWant[, groupings])) / 2)
               } else {
                  mean(pnorm(upper/rtW[, groupings]) - pnorm(lower/rtW[, groupings]))
               }
            } else {
               .Call("eval_nvmix_integral",
                     lower     = as.double(lower),
                     upper     = as.double(upper),
                     groupings = as.integer(groupings),
                     numgroups = as.integer(numgroups), 
                     U         = as.double(U),
                     rtW       = as.double(rtW),
                     rtWant    = as.double(rtWant), 
                     n         = as.integer(current.n),
                     d         = as.integer(d),
                     r         = as.integer(rank),
                     kfactor   = as.integer(k.factor),
                     factor    = as.double(factor),
                     ZERO      = as.double(ZERO),
                     ONE       = as.double(ONE),
                     doant     = as.integer(do.ant))[1]
            } 
         
         ## 2.3 Update RQMC estimates #######################################
         
         rqmc.estimates[b] <-
            if(do.doubling) {
               ## In this case both, rqmc.estimates[b] and
               ## next.estimate depend on n.current points
               (rqmc.estimates[b] + next.estimate) / denom
            } else {
               ## In this case, rqmc.estimates[b] depends on
               ## numiter * n.current points whereas next.estimate
               ## depends on n.current points
               (numiter * rqmc.estimates[b] + next.estimate) / (numiter + 1)
            }
         
      } # end for(b in 1:B)
      ## Update of various variables
      ## Number of function evaluations (* 2 depending on 'do.ant')
      total.fun.evals <- if(do.ant) total.fun.evals + 2 * B * current.n else
         total.fun.evals + B * current.n
      ## Double sample size and adjust denominator in averaging as well as useskip
      if(do.doubling) {
         ## Change denom and useksip (exactly once, in the first iteration)
         if(numiter == 0) {
            denom <- 2
            useskip <- 1
         } else {
            ## Increase sample size n. This is done in all iterations
            ## except for the first two
            current.n <- 2 * current.n
         }
      }
      ## Update error depending on 'do.reltol'
      error <- if(!do.reltol) {
         CI.factor.sqrt.B * sd(rqmc.estimates)
      } else {
         CI.factor.sqrt.B * sd(rqmc.estimates)/mean(rqmc.estimates)
      }
      numiter <- numiter + 1 # update counter.
   } # while()
   ## Finalize
   value <- mean(rqmc.estimates) # calculate the RQMC estimator
   ## Compute absolute and relative error for return 
   abserror <- if(do.reltol){
      relerror <- error
      error * value 
   } else {
      relerror <- error / value # 'error' is absolute error
      error 
   }
   ## Return with correctly set 'mix.stored'
   if(mix.doStore) attr(mix.stored, "numiter.stored") <- numiter.stored # update
   ## Return (if(!mix.doStore), 'mix.stored' is NULL by default) 
   list(value = value, error = error, abserror = abserror, 
        relerror = relerror, numiter = numiter, mix.stored = mix.stored)
}


##' @title Distribution Function of a Multivariate Normal Variance Mixture
##' @param upper (n,d) matrix of upper evaluation limits
##' @param lower (n,d) matrix of lower evaluation limits (<= upper)
##' @param qmix specification of the (mixture) distribution of W. This can be:
##'        1) a character string specifying a supported distribution (additional
##'           arguments of this distribution are passed via '...').
##'        2) a list of length at least one; the first argument specifies
##'           the base name of an existing distribution which can be sampled
##'           with prefix "q", the other elements denote additional parameters
##'           passed to this "rmix" random number generator.
##'        3) a function being interpreted as the quantile function F_W^-.
##' @param loc d-vector (location vector)
##' @param scale (d, d)-covariance matrix (scale matrix)
##' @param standardized logical indicating whether 'scale' is assumed to be a
##'        correlation matrix; if FALSE (default), 'upper', 'lower' and 'scale'
##'        will be normalized.
##' @param control list() of algorithm parameters; see details in ?pnvmix
##' @param verbose logical (or integer: 0 = FALSE, 1 = TRUE, 2 = more output)
##'        indicating whether a warning is given if the required precision
##'        'control$pnvmix.abstol'/'control$pnvmix.reltol' has not been reached.
##' @param ... additional arguments passed to the underlying mixing distribution
##' @return numeric vector with the computed probabilities and attributes "abs. error"
##'         and "rel. error" (error estimates of the RQMC estimator) and "numiter"
##'         (number of iterations)
##' @author Erik Hintz and Marius Hofert
pnvmix <- function(upper, lower = matrix(-Inf, nrow = n, ncol = d), 
                   qmix, rmix, 
                   loc = rep(0, d), scale = diag(d), standardized = FALSE,
                   control = list(), verbose = TRUE,  ...)
{
   ## Checks
   if(!is.matrix(upper)) upper <- rbind(upper) # 1-row matrix if upper is a vector
   n <- nrow(upper) # number of evaluation points
   d <- ncol(upper) # dimension
   if(!is.matrix(lower)) lower <- rbind(lower) # 1-row matrix if lower is a vector
   
   ## Call the more general 'pgnvmix()' with 'groupings = rep(1, d)' 
   pgnvmix(upper = upper, lower = lower, qmix = qmix, rmix = rmix, 
           groupings = rep(1, d), loc = loc, scale = scale, 
           standardized = standardized, control = control, verbose = verbose, ...)
}


##' @title Distribution Function of a Generalized Normal Variance Mixture
##' @param upper (n,d) matrix of upper evaluation limits
##' @param lower (n,d) matrix of lower evaluation limits (<= upper)
##' @param groupings d-vector giving group index of variable i 
##'        (rep(1, d) => NVM dist'n; 1:d => generalized NVM (each component has their own W)
##' @param qmix specification of the (mixture) distributions of W. This can be:
##'        1) a character string specifying a supported distribution (currently
##'        supported are "inverse.gamma" and "pareto" in which case the additional
##'        argument 'df' or 'alpha' needs to be provided via (...); the parameter
##'        must have length = the number of different groups). 
##'        2) a list of length = number of different groups. Each element of this
##'        list can be 
##'        2.1) a list of length at least one; the first argument specifies
##'           the base name of an existing distribution which can be sampled
##'           with prefix "q", the other elements denote additional parameters
##'           passed to this "qfun". 
##'       2,2) a function being interpreted as quantile function F_W_i^-. Additional
##'           arguments can be passed via (...) and will be matched. 
##' @param rmix smiliar to 'qmix' but RNG instead of quantile function; if used,
##'        groupings *must* be rep(1, d) ( => classical NVM dist'n) 
##' @param mean.sqrt.mix E(sqrt(W)) or NULL; the latter if not available
##'        in which case it is estimated by QMC
##' @param loc d-vector (location vector)
##' @param scale (d, d)-covariance matrix (scale matrix)
##' @param standardized logical indicating whether 'scale' is assumed to be a
##'        correlation matrix; if FALSE (default), 'upper', 'lower' and 'scale'
##'        will be normalized.
##' @param control list() of algorithm parameters; see details in ?pnvmix
##' @param verbose logical (or integer: 0 = FALSE, 1 = TRUE, 2 = more output)
##'        indicating whether a warning is given if the required precision
##'        'control$pnvmix.abstol'/'control$pnvmix.reltol' has not been reached.
##' @param ... additional arguments passed to the underlying mixing distribution
##' @return numeric vector with the computed probabilities and attributes "error"
##'         (error estimate of the RQMC estimator) and "numiter"
##'         (number of iterations)
##' @author Erik Hintz and Marius Hofert  
pgnvmix <- function(upper, lower = matrix(-Inf, nrow = n, ncol = d), 
                    groupings = 1:d, qmix, rmix,  
                    loc = rep(0, d), scale = diag(d), standardized = FALSE,
                    control = list(), verbose = TRUE, ...)
{   
   ## Checks
   if(!is.matrix(upper)) upper <- rbind(upper) # 1-row matrix if upper is a vector
   n <- nrow(upper) # number of evaluation points
   d <- ncol(upper) # dimension
   if(!is.matrix(lower)) lower <- rbind(lower) # 1-row matrix if lower is a vector
   if(is.data.frame(upper)) upper <- as.matrix(upper)
   if(is.data.frame(lower)) lower <- as.matrix(lower)
   if(!is.matrix(scale)) scale <- as.matrix(scale)
   stopifnot(dim(lower) == c(n, d), length(loc) == d, # 'mean.sqrt.mix' is tested in pnvmix1()
             dim(scale) == c(d, d), length(groupings) == d)
   ## Prepare mixing variable 
   ## Set [q/r]mix to NULL if not provided (needed for 'get_mix_()')
   if(!hasArg(qmix)) qmix <- NULL
   if(!hasArg(rmix)) rmix <- NULL
   mix_list      <- get_mix_(qmix = qmix, rmix = rmix, groupings = groupings, callingfun = "pnvmix", ... )
   mix_          <- mix_list$mix_ # function(u) or function(n) depeneding on 'use.q'
   special.mix   <- mix_list$special.mix
   use.q         <- mix_list$use.q
   mean.sqrt.mix <- mix_list$mean.sqrt.mix
   numgroups     <- length(unique(groupings)) # number of groups 
   ## Check 'method': The default depends on wether 'rmix' or 'qmix' was provided
   meth.prov <- if(!use.q & length(control)>0 & any(names(control) == "method")){
      control$method
   } else NA # provided method
   ## Get/overwrite defaults 
   control <- get_set_param(control)
   do.ant <- if(!use.q){
      if(!is.na(meth.prov) & meth.prov != "PRNG") # method other than 'PRNG' provided
         warning("When 'rmix' provided, currently only available 'method' is 'PRNG'")
      control$method <- "PRNG" # set method to "PRNG" 
      FALSE
   } else control$pnvmix.do.ant 
   ## Grab the following variables for readability
   method    <- control$method
   increment <- control$increment
   ## Generate 'B' seeds for sobol(.., seed = seeds[b])
   seeds <- if(method == "sobol") sample(1:1e3, control$B) else NULL 
   ## Absolute or relative precision required?
   tol <- if(is.null(control$pnvmix.abstol)) {
      ## Set tol to <0 so that algorithm runs until 'fun.eval[2]'
      ## function evaluations or 'max.iter.rqmc' iterations are exhausted
      do.reltol <- FALSE
      -42
   } else if(is.na(control$pnvmix.abstol)) { # if 'abstol = NA' use relative precision
      do.reltol <- TRUE
      control$pnvmix.reltol
   } else { # otherwise use absolute precision (default)
      do.reltol <- FALSE
      control$pnvmix.abstol
   }
   ## If d = 1, call 'pnvmix1d()' which is truly *vectorized* (=> no variable reordering)
   ## Case of normal and t dist'ns handled below with 'pnorm()' and 'pt()' 
   if(d == 1 & (is.na(special.mix) || special.mix == "pareto")) {
      return(pnvmix1d(upper = as.vector(upper), lower = as.vector(lower),
                      mix_ = mix_, use.q = use.q, do.ant = do.ant, 
                      loc = loc, scale = as.numeric(scale),
                      standardized = standardized, method = method,tol = tol, 
                      do.reltol = do.reltol, CI.factor = control$CI.factor,
                      fun.eval = control$fun.eval, max.iter.rqmc = control$max.iter.rqmc,
                      increment = increment, B = control$B, verbose = verbose,
                      seeds = seeds))
   }
   ## Grab / approximate mean.sqrt.mix, which will be needed for preconditioning
   ## in pnvmix1(). This only depends on 'qmix', hence it is done (once) here in pnvmix.
   if(control$precond & d > 2) {
      if(is.null(mean.sqrt.mix) & !is.null(control$mean.sqrt.mix))
         mean.sqrt.mix <- control$mean.sqrt.mix
      if(is.null(mean.sqrt.mix)){
         ## Not provided => estimate it 
         mean.sqrt.mix <- if(use.q){
            colMeans(as.matrix(sqrt(mix_(qrng::sobol(n = 2^12, d = 1, randomize = TRUE)))))
         } else colMeans(as.matrix(sqrt(mix_(2^12))))
         mean.sqrt.mix <- mean.sqrt.mix[groupings] 
      } else {
         ## Check if provided 'mean.sqrt.mix' has the correct dimension
         ## (length == 1 for compatibility with 'pnvmix()')
         stopifnot(length(mean.sqrt.mix) == d | length(mean.sqrt.mix) == 1)
      }
      ## Check if provided/approximated 'mean.sqrt.mix' is strictly positive
      if(any(mean.sqrt.mix <= 0))
         stop("'mean.sqrt.mix' has to be positive (possibly after being generated in pnvmix())")
   }
   ## Build temporary result list
   res1 <- vector("list", length = n) # results from calls of pnvmix1()
   ## Deal with NA
   NAs <- apply(is.na(lower) | is.na(upper), 1, any) # at least one NA => use NA => nothing left to do
   ## Remove 'loc':
   if(any(loc != 0)) {
      lower <- lower - matrix(rep(loc, n), ncol = d, byrow = TRUE)
      upper <- upper - matrix(rep(loc, n), ncol = d, byrow = TRUE)
   }
   ## Get (lower triangular) Cholesky factor of 'scale'
   factor.obj   <- cholesky_(scale, tol = control$cholesky.tol)
   factor       <- factor.obj$C
   factor.diag  <- factor.obj$D
   ## Obtain rank
   D.zero       <- (factor.diag==0)
   rank         <- d - sum(D.zero)
   ## In case of a singular matrix, need to reorder 'factor':
   if(rank < d) {
      if(verbose) warning("Provided 'scale' is singular.")
      ## In each row i, get minimal index j such that factor[i,k]=0 for all k>j
      length.rows <- apply(factor, 2, function(i) which.max( i != 0))
      order.length.rows <- order(length.rows) # needed later to sort 'low' and 'up'
      length.rows.sorted <- length.rows[order.length.rows]
      factor <- factor[order.length.rows, ] # sort factor
      lower <- lower[, order.length.rows, drop = FALSE]
      upper <- upper[, order.length.rows, drop = FALSE]
      ## Now factor has the form ( * non-zero element, ? any element)
      ##    * 0 0 ......... 0
      ##    ...
      ##    * 0 0 ......... 0
      ##    ? * 0 ......... 0
      ##    ...
      ##    ? * 0 ......... 0
      ##    ...
      ##    ? ... ? * 0 ... 0
      ## Need to standardize: Divide each row by its rightmost ' * ' (which is !=0)
      row.scales   <- factor[cbind(1:d, length.rows.sorted)] # vector of ' * ' elements
      factor       <- diag(1/row.scales) %*% factor
      ## Standardize 'lower' and 'upper' accordingly (also works for +/- Inf)
      lower <- matrix(rep(1/row.scales, each = n), ncol = d, nrow = n, byrow = FALSE) *
         lower
      upper <- matrix(rep(1/row.scales, each = n), ncol = d, nrow = n, byrow = FALSE) *
         upper
      ## Need to swap columns in 'lower'/'upper' multiplied by negative 'row.scales':
      if(any(row.scales < 0)) {
         which.row.scales.neg <- which(row.scales < 0)
         temp.lower <- lower[, which.row.scales.neg]
         lower[, which.row.scales.neg] <- upper[, which.row.scales.neg]
         upper[, which.row.scales.neg] <- temp.lower
      }
      ## Now we have the form as above, with '*' replaced by '1'. Remove 0 columns:
      factor <- factor[, -which(D.zero)]
      ## Get 'height' of each step:
      k.factor <- sapply(unique(length.rows.sorted), function(i) sum(length.rows == i))
      ## no update of'scale' as in singular case, no preconditioning is performed!
   } else {
      ## In case of non-singular 'scale', each of the d steps in factor has height 1:
      k.factor <- rep(1, d)
      ## Standardize 'scale', 'lower', 'upper', 'factor' ('loc' was already taken care of)
      if(!standardized) {
         row.scales    <- diag(scale) # diagonal of cholesky factor
         rt.scales.inv <- sqrt(1/row.scales) # d-vector
         ## The following is equivalent but faster than
         ## 'scale <- diag(row.scales.inv) %*% scale %*% diag(row.scales.inv)'
         ## [http://comisef.wikidot.com/tip:diagonalmatrixmultiplication]
         scale  <- outer(rt.scales.inv, rt.scales.inv) * scale # standardize scale
         factor <- matrix(rep(rt.scales.inv, each = d),
                          ncol = d, byrow = TRUE) * factor # standardized cholesky
         ## Now standardize scale as above: Also works for +/- Inf
         lower <- matrix(rep(rt.scales.inv, each = n), ncol = d, nrow = n, byrow = FALSE) *
            lower
         upper <- matrix(rep(rt.scales.inv, each = n), ncol = d, nrow = n, byrow = FALSE) *
            upper
      }
   }
   ## Logical indicating if we're dealing with a multivariate normal dist'n
   is.const.mix = (!is.na(special.mix) & special.mix == "constant")
   ## If n>1, 'mix.realixations' is a matrix with B columns storing evalauations of mix_,
   ## will be created by 'pnvmix1()' below, if needed 
   mix.stored <- NULL 
   ## Used only when more than one estimation required for a non-normal dist'n
   ## and when sobol(..., seed = ...) is used 
   mix.doStore <- !is.const.mix & (n > 1) & (method == "sobol")
   
   ## Loop over observations ##################################################
   reached <- rep(TRUE, n) # indicating whether 'tol' has been reached in the ith integration bounds (needs default TRUE)
   for(i in seq_len(n)) {
      if(NAs[i]) {
         res1[[i]] <- 
            list(value = NA, error = 0, abserror = 0, relerror = 0, numiter = 0)
         next # => 'reached' already has the right value
      }
      ## Pick out ith row of lower and upper
      low <- lower[i,]
      up  <- upper[i,]
      ## Deal with equal bounds (result is 0 as integral over null set)
      if(any(low == up)) {
         res1[[i]] <- 
            list(value = 0, error = 0, abserror = 0, relerror = 0, numiter = 0)
         next # => 'reached' already has the right value
      }
      ## Deal with case where components of both low and up are Inf
      lowFin <- is.finite(low)
      upFin  <- is.finite(up)
      lowupFin <- lowFin | upFin # at least one finite limit
      ## Only for ungrouped full rank case for efficiency as ow rank, k.factor etc change
      ## completely and would need to be recalculated
      if(any(!lowupFin) & rank == d) {
         ## Update low, up
         low <- low[lowupFin]
         up  <- up [lowupFin]
         ## Grab (new) dimension. If 0, then all upper are +Inf, all lower are -Inf
         ## => Return 0
         dFin <- length(low) # Update dimension and 'k.factor'
         k.factorFin <- rep(1, dFin )
         if(dFin  == 0) {
            res1[[i]] <- list(value = 1, error = 0, abserror = 0, relerror = 0, numiter = 0)
            next
         }
         ## Update groupings
         groupingsFin <- groupings[lowupFin] 
         numgroupsFin <- length(unique(groupingsFin))
         ## Update scale and factor (for preconditioning)
         scaleFin <- scale[lowupFin, lowupFin, drop = FALSE]
         factorFin <- t(chol(scaleFin)) # Cholesky factor changes
      } else {
         ## If no such component exists, variables correctly 
         dFin <- d
         scaleFin <- scale 
         factorFin <- factor
         k.factorFin <- k.factor 
         groupingsFin <- groupings 
         numgroupsFin <- numgroups
      }
      ## If d = 1, deal with multivariate normal, and t via pnorm() and pt()
      ## Note that everything has been standardized.
      if(dFin == 1 & !is.na(special.mix) & numgroupsFin == 1) {
         if(is.const.mix) {
            value <- pnorm(up) - pnorm(low)
            res1[[i]] <- list(value = value, error = 0, abserror = 0, 
                              relerror = 0, numiter = 0)
            next
         } else if(special.mix == "inverse.gamma") {
            df_ <- mix_list$param[groupingsFin[1]] 
            value <- pt(up, df = df_ ) - pt(low, df = df_ )
            res1[[i]] <- list(value = value, error = 0, abserror = 0, 
                              relerror = 0, numiter = 0)
            next
         }
      }
      ## Compute result for ith row of lower and upper (in essence,
      ## because of the preconditioning, one has to treat each such
      ## row separately)
      ## Note: Only lower triangular part of 'factor' needed in pnvmix1()
      tmp <- 
         pnvmix1(up, lower = low, groupings = groupingsFin, mix_ = mix_,
                 use.q = use.q, do.ant = do.ant,
                 is.const.mix = is.const.mix,
                 mean.sqrt.mix = mean.sqrt.mix, loc = loc, 
                 scale = scaleFin,
                 factor = factorFin[lower.tri(factorFin, diag = TRUE)], 
                 k.factor = k.factorFin, method = method,
                 precond = control$precond, tol = tol, do.reltol = do.reltol,
                 CI.factor = control$CI.factor, fun.eval = control$fun.eval,
                 max.iter.rqmc = control$max.iter.rqmc, increment = increment,
                 B = control$B, verbose = as.logical(verbose), seeds = seeds, 
                 mix.stored = mix.stored, mix.doStore = mix.doStore, 
                 maxiter.stored = control$maxiter.stored, df = df)  
      res1[[i]] <- if(mix.doStore){
         mix.stored <- tmp[[6]]
         tmp[-6]
      } else tmp 
      ## Check if desired precision was reached
      ## Note: 'error' is either relative or absolute, depending on 'do.reltol' 
      reached[i] <- res1[[i]]$error <= tol 
      if(verbose >= 2 && !reached[i])
         warning(sprintf("Tolerance not reached for pair %d of integration bounds; consider increasing 'fun.eval[2]' and 'max.iter.rqmc'", i))
   } # for()
   if(verbose == 1 && any(!reached)) { # <=> verbose == 1
      ii <- which(!reached) # (beginning of) indices
      strng <- if(length(ii) > 6) {
         paste0(paste(head(ii), collapse = ", "), ",...")
      } else {
         paste(ii, collapse = ", ")
      }
      warning("Tolerance not reached for pair(s) ",strng," of integration bounds; consider increasing 'fun.eval[2]' and 'max.iter.rqmc'")
   }
   ## Return
   res <- vapply(res1, function(r) r$value, NA_real_)
   attr(res, "abs. error")   <- vapply(res1, function(r) r$abserror, NA_real_)
   attr(res, "rel. error")   <- vapply(res1, function(r) r$relerror, NA_real_)
   attr(res, "numiter") <- vapply(res1, function(r) r$numiter, NA_real_)
   res
}


##' @title Distribution Function of a 1d Normal Variance Mixture
##' @param upper n vector of upper evaluation limits
##' @param lower n vector of lower evaluation limits (<= upper)
##' @param mix_ function(u) or function(n) interpreted as quantile function 
##'        or RNG for the mixing variable
##' @param use.q logical if 'mix_' is a quantile function 
##' @param do.ant logical if antithetic variates are being used 
##' @param loc numeric (location)
##' @param scale numeric (scale), >0
##' @param standardized logical indicating whether 'scale' is assumed to be a
##'        correlation matrix; if FALSE (default), 'upper', 'lower' and 'scale'
##'        will be normalized.
##' @param method see details in ?pnvmix
##' @param tol absolute/relative error tolerance
##' @param do.reltol logical if relative precision is required
##' @param CI.factor see details in ?pnvmix
##' @param fun.eval see details in ?pnvmix
##' @param max.iter.rqmc see details in ?pnvmix
##' @param increment see details in ?pnvmix
##' @param B see details in ?pnvmix
##' @param verbose see details in ?pnvmix
##' @param ... see details in ?pnvmix
##' @return numeric vector with the computed probabilities and attributes
##'         "error" (error estimate of the RQMC estimator) and "numiter"
##'         (number of iterations)
##' @author Erik Hintz and Marius Hofert
pnvmix1d <- function(upper, lower = rep(-Inf, n), mix_, use.q = TRUE, do.ant = TRUE,
                     loc = 0, scale = 1,
                     standardized = FALSE, method = "sobol",
                     tol = 1e-3, do.reltol = FALSE,
                     CI.factor = 3.5, fun.eval = c(2^6, 1e8), max.iter.rqmc = 15,
                     increment = "doubling", B = 15, verbose = FALSE, 
                     seeds = NULL)
{
   ## 1 Setup #################################################################
   
   n <- length(upper)
   if(!is.vector(upper)) upper <- as.vector(upper)
   if(!is.vector(lower)) lower <- as.vector(lower)
   ## Standardize
   if(!standardized) {
      upper <- (upper - loc) / sqrt(scale)
      lower <- (lower - loc) / sqrt(scale)
   }
   ## Error is calculated as CI.factor * sd( estimates) / sqrt(B).
   CI.factor.sqrt.B <- CI.factor / sqrt(B)
   ## Grab the number of sample points for the first iteration
   current.n <- fun.eval[1] #
   ## Matrix to store the B * length(upper) RQMC estimates
   rqmc.estimates <- matrix(0, ncol = n, nrow = B)
   ## Initialize error to > 'tol' so that we can enter the while loop below
   error <- tol + 42
   ## Initialize the total number of function evaluations
   total.fun.evals <- 0
   ## Initialize counter of the number of iterations in the while loop below
   numiter <- 0
   ## It may happen that qnorm(u) for u too close to 1 (or 0) is evaluated; in those
   ## cases, u will be replaced by ONE and ZERO which is the largest (smallest) number
   ## different from 1 (0) such that qnorm(u) is not +/- Inf
   ZERO <- .Machine$double.eps # for symmetry reasons (-8/+8), use this instead of .Machine$double.xmin
   ONE  <- 1-.Machine$double.neg.eps
   ## If method == sobol, we want the same random shifts in each iteration below
   if(method == "sobol") {
      if(is.null(seeds)) seeds <- sample(1:1e3, B) else stopifnot(length(seeds) == B) 
   }
   ## Additional variables needed if the increment chosen is "doubling"
   if(increment == "doubling") {
      if(method == "sobol") useskip <- 0
      denom <- 1
   }
   
   ## 2 Main while() loop #####################################################
   
   ## while() runs until precision 'tol' is reached or the number of function
   ## evaluations exceed fun.eval[2]. In each iteration, B RQMC estimates of
   ## the desired probability are calculated.
   while(max(error) > tol & total.fun.evals < fun.eval[2] & numiter < max.iter.rqmc)
   {

            ## Get B RQCM estimates
      for(b in 1:B)
      {
         ## 2.1 Get the point set ###########################################
         ## U will contain realizations of 1 / sqrt(W):
         U <- if(use.q){
            U <- switch(method,
                        "sobol" = {
                           if(increment == "doubling") {
                              qrng::sobol(n = current.n, d = 1,
                                          randomize = "digital.shift",
                                          seed = seeds[b], 
                                          skip = (useskip * current.n))
                           } else {
                              qrng::sobol(n = current.n, d = 1,
                                          randomize = "digital.shift",
                                          seed = seeds[b], 
                                          skip = (numiter * current.n))
                           }
                        },
                        "ghalton" = {
                           qrng::ghalton(n = current.n, d = 1,
                                         method = "generalized")
                        },
                        "PRNG" = {
                           matrix(runif( current.n ), ncol = 1)
                        })
            if(do.ant) 1 / cbind( sqrt(mix_(U)), sqrt(mix_(1-U))) else sqrt(mix_(U))
         } else {
            ## 'mix_' is a RNG; do.ant = FALSE in this case 
            1 / sqrt(mix_(current.n))
         }
         
         ## 2.2 Evaluate the integrand at the (next) point set ##############
         next.estimate <- 
            if(do.ant){
               colMeans( 
                  (pnorm(outer( U[,1], upper)) - pnorm(outer(U[,1], lower)) +
                      pnorm(outer(U[,2], upper)) - pnorm(outer(U[,2], lower)))/2)
            } else {
               colMeans(pnorm(outer( U[,1], upper)) - pnorm(outer(U[,1], lower)))
            }
         
         ## 2.3 Update RQMC estimates #######################################
         rqmc.estimates[b, ] <-
            if(increment == "doubling") {
               ## In this case both, rqmc.estimates[b] and
               ## next.estimate depend on n.current points
               (rqmc.estimates[b, ] + next.estimate) / denom
            } else {
               ## In this case, rqmc.estimates[b] depends on
               ## numiter * n.current points whereas next.estimate
               ## depends on n.current points
               (numiter * rqmc.estimates[b, ] + next.estimate) / (numiter + 1)
            }
      } # end for(b in 1:B)
      
      ## Update of various variables
      ## Number of function evaluations
      ## (* 2 since antithetic variates are used in eval_nvmix_integral())
      total.fun.evals <- if(do.ant) total.fun.evals + 2 * B * current.n else
         total.fun.evals + B * current.n
      ## Double sample size and adjust denominator in averaging as well as useskip
      if(increment == "doubling") {
         ## Change denom and useksip (exactly once, in the first iteration)
         if(numiter == 0) {
            denom <- 2
            useskip <- 1
         } else {
            ## Increase sample size n. This is done in all iterations
            ## except for the first two
            current.n <- 2 * current.n
         }
      }
      ## Update error depending on 'do.reltol'
      error <- if(!do.reltol) { # absolute error
         CI.factor.sqrt.B * apply(rqmc.estimates, 2, sd)
      } else { # relative error
         CI.factor.sqrt.B * apply(rqmc.estimates, 2, sd)/colMeans(rqmc.estimates)
      }
      numiter <- numiter + 1 # update counter
   } # while()
   
   ## 3 Finalize ##############################################################
   
   ## Check if error tolerance reached and print warnings accordingly
   reached <- (error <= tol)
   if(any(!reached) && verbose > 0) {
      ii <- which(!reached)
      if(verbose == 1) {
         strng <- if(length(ii) > 6) {
            paste0(paste(head(ii), collapse = ", "), ",...")
         } else {
            paste(ii, collapse = ", ")
         }
         warning("Tolerance not reached for pair(s) ",strng," of integration bounds; consider increasing 'fun.eval[2]' and 'max.iter.rqmc' in the 'control' argument.")
      } else {
         for(i in 1:length(ii)) {
            warning(sprintf("Tolerance not reached for pair %d of integration bounds; consider increasing 'fun.eval[2]' and 'max.iter.rqmc' in the 'control' argument", ii[i]))
         }
      }
   }
   ## Obtain estimates
   res <- colMeans(rqmc.estimates)
   ## Compute absolute and relative error for return 
   abserror <- if(do.reltol){
      relerror <- error
      error * res 
   } else {
      relerror <- error / res # 'error' is absolute error
      error 
   }
   ## Assign attributes and return 
   attr(res, "abs. error") <- abserror
   attr(res, "rel. error") <- relerror
   attr(res, "numiter") <- numiter
   res
}