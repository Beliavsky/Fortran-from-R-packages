### dnvmix() ###################################################################

##' @title Merge two matrices that are sorted in one of their columns
##' @param A (n, k) matrix, column 'sort.by' sorted increasingly 
##' @param B (m, k) matrix, column 'sort.by' sorted increasingly 
##' @param sort.by column which is sorted and should be sorted by 
##' @return (n+m, k) matrix with rows from A and B so that 
##'          column 'sort.by' is sorted increasingly 
##' @author Erik Hintz 
##' @note k = 1 is allowed 
merge_m <- function(A, B, sort.by = 1){
   ## Checks
   if(!is.matrix(A)) A <- as.matrix(A)
   if(!is.matrix(B)) B <- as.matrix(B) 
   dimA <- dim(A)
   dimB <- dim(B)
   stopifnot(dimA[2] == dimB[2]) # A and B have the same number of columns 
   n <- dimA[1] # number of rows of A
   m <- dimB[1] # number of rows of B
   res <- matrix(NA, ncol = dimA[2], nrow = n + m) # result matrix
   p.A <- 1 # pointer to current element in 'A'
   p.B <- 1 # pointer to current element in 'B'
   for(i in 1:(n+m)){
      res[i, ] <- if(A[p.A, sort.by] < B[p.B, sort.by]){
         p.A <- p.A + 1 
         A[p.A-1, , drop = FALSE]
      } else {
         p.B <- p.B + 1
         B[p.B-1, , drop = FALSE]
      }
      ## Arrived at the end of 'A'
      if(p.A == n + 1){
         res[(i+1):(n+m), ] <- B[p.B:m, , drop = FALSE]
         break
      }
      ## Arrived at the end of 'B'
      if(p.B == m + 1){
         res[(i+1):(n+m), ] <- A[p.A:n, , drop = FALSE]
         break
      }
   }
   res
}

### dnvmix() ###################################################################

##' @title Merge two matrices that are sorted in one of their columns
##' @param A (n, k) matrix, column 'sort.by' sorted increasingly 
##' @param B (m, k) matrix, nothing sorted 
##' @param sort.by column which is sorted and should be sorted by 
##' @return (n+m, k) matrix with rows from A and B so that 
##'          column 'sort.by' is sorted increasingly 
##' @author Erik Hintz 
##' @note k = 1 is allowed 
merge_m2 <- function(A, B, sort.by = 1){
   ## Checks
   if(!is.matrix(A)) A <- as.matrix(A)
   if(!is.matrix(B)) B <- as.matrix(B) 
   dimA <- dim(A)
   dimB <- dim(B)
   stopifnot(dimA[2] == dimB[2]) # A and B have the same number of columns 
   n <- dimA[1] # number of rows of A
   m <- dimB[1] # number of rows of B
   res <- matrix(NA, ncol = dimA[2], nrow = n + m) # result matrix
   p.A <- 1 # pointer to current element in 'A'
   p.B <- 1 # pointer to current element in 'B'
   for(i in 1:(n+m)){
      res[i, ] <- if(A[p.A, sort.by] < B[p.B, sort.by]){
         p.A <- p.A + 1 
         A[p.A-1, , drop = FALSE]
      } else {
         p.B <- p.B + 1
         B[p.B-1, , drop = FALSE]
      }
      ## Arrived at the end of 'A'
      if(p.A == n + 1){
         res[(i+1):(n+m), ] <- B[p.B:m, , drop = FALSE]
         break
      }
      ## Arrived at the end of 'B'
      if(p.B == m + 1){
         res[(i+1):(n+m), ] <- A[p.A:n, , drop = FALSE]
         break
      }
   }
   res
}

##' @title Exp - log trick for log(sum_i exp(a_i))
##' @param M (n1, n2) matrix
##' @return n2-vector log(colSums(exp(M)))
##' @author Erik Hintz
##' @note NO checking is done for efficiency reasons
logsumexp <- function(M)
{
   n.rows <- nrow(M)
   n.cols <- ncol(M)
   if(!is.matrix(M)) M <- rbind(M) 
   cmax <- apply(M, 2, max)
   cmax + log(.colSums(exp( M - rep(cmax, each = n.rows)), m = n.rows, n = n.cols))
}

##' @title Adaptive RQMC Method to estimate the log-density of an NVMIX
##'        distribution
##' @param qW function of one variable specifying the quantile function of W.
##' @param maha2.2 squared maha-distances divided by 2. Assumed to be *sorted*
##' @param lconst see ?densmix_
##' @param d dimension of the underlying NVM_d distribution
##' @param k k = d for estimating log-density; k = d+2 for log-weights for fitnvmix
##' @param control see ?dnvmix
##' @param UsWs matrix; each row is of the form (u, qW(u)) where u in (0,1)
##' @return list of three ($ldensities, $error, $numiter); each vector of
##'         length(maha2.2)
##' @author Erik Hintz and Marius Hofert
densmix_adaptrqmc <- function(qW, maha2.2, lconst, d, k = d, control, UsWs)
{
   ## 1 Set-up #################################################################
   
   numiter   <- 0 # counter for the number of iterations
   ONE       <- 1-.Machine$double.neg.eps
   ZERO      <- .Machine$double.neg.eps
   tol.bisec <- control$dnvmix.tol.bisec # vector of length 3
   ## [1]/[2]/[3] => tolerance on 'u' / 'W' / 'log-integrand' for bisections
   ## Initialize output object
   n          <- length(maha2.2)
   ldensities <- rep(NA, n)
   error      <- rep(NA, n)
   numiters   <- rep(NA, n)
   ## Grab 'UsWs'; store them in vectors and sort
   stopifnot(is.matrix(UsWs), dim(UsWs)[2] == 2) # sanity check 
   ordering.U <- order(UsWs[, 1, drop = FALSE])
   numObs     <- dim(UsWs)[1] + 2 # will add 'ZERO' and 'ONE'
   ## Set up matrix of the form  U | qW(U) | l.integrand(qW((U))
   U.W.lint       <- matrix(NA, ncol = 3, nrow = numObs)
   U.W.lint[, 1]  <- c(ZERO, UsWs[ordering.U, 1], ONE)
   U.W.lint[, 2]  <- c(max(qW(ZERO), ZERO), UsWs[ordering.U, 2], qW(ONE))
   ## Check if W is bounded
   W.max <- qW(1) # <= inf
   W.min <- qW(0) # >= 0
   W.min_gen <- U.W.lint[1, 2] # qW(ZERO) >= 0; smallest 'W' we can generate
   W.max_gen <- U.W.lint[numObs, 2] # qW(ONE) <= Inf; largest 'W' we can generate
   isWfin <- is.finite(W.max) 
   
   ## 2 Main loop over Mahalanobis distances ##################################
   
   ## For *each* maha2.2 value separately, find optimal integration region 
   ## around the 'peek', integrate there via RQMC and outside via crude
   ## trapezoidal rules
   for(ind in 1:n) {
      curr.maha2.2 <- max(maha2.2[ind], ZERO) # avoid maha2.2 = 0
      curr.lconst  <- lconst[ind]
      ## Initialize various quantities
      curr.error    <- NA
      ldens.right   <- NA # log-density in (u_r, 1)
      ldens.left    <- NA # log-density in (0, u_l)
      ldens.stratum <- NA # log-density in (u_l, u_r)
      uLuR          <- rep(NA, 2) # c('u.left', 'u.right') for later
      ## Determine location of maximum
      int.argmax.u  <- NA # u* such that qmix(u*) = W* where integrand(u*) is max
      int.argmax.w  <- 2*curr.maha2.2/k # value of W = qmix(u*) at max assuming W is *unbounded*
      ## Check special cases
      peekRight <- FALSE # logical if peek beyond ONE (can't sample there)
      peekLeft  <- FALSE # logical if peek below ZERO (can't sample there)
      if(int.argmax.w <= W.min_gen) {
         ## Case 1: 'peek' too close to 0
         int.argmax.w <- if(int.argmax.w <= W.min){
            ## Can only happen when W.min > 0, ie when support(W) != (0, Inf) 
            W.min
         } else {
            peekLeft <- TRUE 
            W.min_gen
         }
         uLuR[1]      <- 0 # u.left = 0
         int.argmax.u <- 0
         ldens.left   <- -Inf # (0, u_l) = (0, 0) => empty region here
      } else if(int.argmax.w >= W.max_gen){
         ## Case 2: 'peek' too close to 1
         int.argmax.w <- if(int.argmax.w >= W.max){
            W.max # W has bounded support!
         } else {
            peekRight <- TRUE
            W.max_gen
         }
         uLuR[2]      <- 1 # u.right = 1
         int.argmax.u <- 1
         ldens.right  <- -Inf # (u_r, 1) = (1, 1) => empty region here 
      }
      outofreach.w  <- peekLeft ||  peekRight 
      ## Realizations of log-integrand
      U.W.lint[,3]  <- curr.lconst - log(U.W.lint[,2])*k/2 -
         curr.maha2.2/U.W.lint[,2] # sorted according to ordering(U)!
      ## Value of the theoretical maximum of the log-integrand
      l.int.theo.max <- curr.lconst - log(int.argmax.w)*k/2 -
         curr.maha2.2/int.argmax.w
      ## Threshold: Will only use RQMC where l.integrand > l.tol.int.lower
      l.tol.int.lower <- 
         min(max(log(control$dnvmix.tol.int.lower), l.int.theo.max - 
                    control$dnvmix.order.lower*log(10)), 0)
      
      ## 2.1 Find u* = argmax_u{integrand(u)} ################################
      
      ## If !NA => Already set above 
      if(is.na(int.argmax.u)) {
         ## Find u* such that g(u*) ~= g_max where g is the integrand
         ## <=> find u* such that qW(U*) = int.argmax.w via binary search
         ## => find starting values (u1, u2): qW(u1) < int.argmax.w & qW(u2) >= int.argmax.w
         ind.first.bigger <- which(U.W.lint[, 2] >= int.argmax.w)[1]
         u1u2 <- c(U.W.lint[ind.first.bigger - 1, 1], 
                   U.W.lint[ind.first.bigger, 1])
         ## Set up bisection 
         convd   <- FALSE
         numiter <- 0
         ## Matrix to store additional u's and W's generated:
         addvalues <- matrix(NA, ncol = 2, nrow = control$dnvmix.max.iter.bisec)
         index.first.u  <- which(U.W.lint[,1] >= u1u2[1])[1] # will insert after this index
         while(!convd && numiter < control$dnvmix.max.iter.bisec) {
            numiter <- numiter + 1
            ## Next point to check (midpoint of u1u2):
            u.next <- mean(u1u2)
            diff <- ((w.next <- max(ZERO, qW(u.next))) - int.argmax.w)
            ## Store u.next and quanitle:
            addvalues[numiter, ] <- c(u.next, w.next)
            ## Update u1u2 depending on sign of 'diff' and check convergence:
            if(diff > 0) u1u2[2] <- u.next else u1u2[1] <- u.next
            convd <- (abs(diff) < tol.bisec[2]) || (diff(u1u2) < tol.bisec[1])
         }
         int.argmax.u <- u.next
         ## Add additional values of 'U', 'W' and 'l.integrand' in 'U.W.lint'
         ## while *preserving the order*.
         ## Order of first column = order of all columns
         order.add.us <- order(addvalues[1:numiter, 1, drop = FALSE])
         WsToAdd <- addvalues[order.add.us, 2, drop = FALSE] # needed twice
         U.W.lint <- 
            rbind(U.W.lint[1:index.first.u,], 
                  cbind(addvalues[order.add.us, 1, drop = FALSE], WsToAdd,
                        curr.lconst - log(WsToAdd)*k/2 - curr.maha2.2/WsToAdd),
                  U.W.lint[(index.first.u+1):numObs,])
         ## Update length (=nrow) of 'U.W.lint'
         numObs <- numObs + numiter
      } 
      
      ## 2.2 Find stratum (u.left, u.right) over which RQMC is applied #######
      
      ## 2.2.1 Find "candidates" for bisection ###############################
      
      ## Need to distinguish multiple special cases:
      greater.thshold <- (U.W.lint[, 3] > l.tol.int.lower)
      if(any(greater.thshold)) {
         ## In this case there is at least one obs exceeding the threshold
         ind.greater <- which(greater.thshold)
         candid.left  <- c(max(0, U.W.lint[ind.greater[1] - 1, 1]), # < thshold
                           min(int.argmax.u, U.W.lint[ind.greater[1], 1]) ) # > thshold
         ## Note: If ind.greater[1] = 1 => U[ind.greater[1] - 1] = numeric(0)
         ## => max(0, U[ind.greater[1] - 1]) = 0
         ## Need several cases for 'candid.right'.
         last.ind.greater <- ind.greater[length(ind.greater)] # last index: g > thshold
         candid.right <- if(last.ind.greater < numObs) {
            ## Let u^ be largest u >= u* such that g(u^) > threshold (exists)
            ## Case 1:
            ## There is u' with u'>=u^ and g(u') < threshold
            ## =>  'u.right' in (u^, u')
            c(max(U.W.lint[last.ind.greater, 1], int.argmax.u), 
              U.W.lint[last.ind.greater + 1, 1] )
         } else {
            ## Case 2: the last point is such that g(u) > threshold => g(ONE) > threshold
            ## => u.right = 1 no matter where the peek is, but the handling of the region
            ## (ONE, 1) *does* depend on where the peek is.
            ldens.right <- if(peekRight) {
               ## Area with peek right of ONE => can't simulate there
               ## => stratify from u.left to u.right = 1 (=ONE)
               ## => use *crude* approximation for region (ONE, 1) where max is
               log(ZERO) + l.int.theo.max - log(2)
            } else {
               ## Area with the peek is left of ONE but theoretical 'u.right' is
               ## out of reach due to machine precision.
               ## => stratify from u.left to u.right = 1 (=ONE)
               ## => do nothing in  (ONE, 1)
               -Inf
            }
            ## See above: u.right = 1
            c(1, 1)
         }
      } else if(peekRight) {
         ## In this case, none of the observations is > threshold and peek close to 1:
         ## Will use *all* obs to estimate in (0, ONE) and crude approximation in (ONE, 1)
         ## where max is. NO stratification!
         candid.left    <- c(1, 1)
         candid.right   <- c(1, 1)
         ldens.right    <- log(.Machine$double.neg.eps) + l.int.theo.max - log(2)
         ldens.stratum  <- -Inf
      } else if(peekLeft) {
         candid.left    <- c(0, 0)
         candid.right   <- c(0, 0)
         ldens.left     <- log(.Machine$double.neg.eps) + l.int.theo.max - log(2)
         ldens.stratum <- -Inf
      } else {
         candid.left  <- c(0, int.argmax.u)
         candid.right <- c(int.argmax.u, 1)
      }
      
      ## 2.2.2 Bisection to find 'u.left' and 'u.right' ######################
      
      candids <- rbind(candid.left, candid.right)
      for(i in 1:2) {
         curr.candid <- candids[i, ]
         if(!is.na(uLuR[i])) next # we already set u.left or u.right
         uLuR[i] <- if(diff(curr.candid) <= tol.bisec[1]) {
            curr.candid[2]
         } else {
            ## Use bisection similar to the one used to find u*
            convd <- FALSE
            numiter <- 0
            ## Matrix to store additional u's, W's, l.integrand's generated:
            addvalues <- matrix(NA, ncol = 3,
                                nrow = control$dnvmix.max.iter.bisec)
            while(!convd && numiter < control$dnvmix.max.iter.bisec) {
               numiter <- numiter + 1
               ## Next point to check
               u.next <- mean(curr.candid)
               w.next <- max(qW(u.next), ZERO)
               diff   <- (l.int.next <- (curr.lconst -log(w.next)*k/2 -
                                            curr.maha2.2/w.next)) - l.tol.int.lower
               ## Store values generated
               addvalues[numiter, ] <- c(u.next, w.next, l.int.next)
               ## Update 'curr.candid' depending on sign of 'diff' and check convergence:
               if(diff > 0) curr.candid[2] <- u.next else curr.candid[1] <- u.next
               convd <- 
                  (abs(diff) < tol.bisec[3]) || (diff(curr.candid) < tol.bisec[1])
            }
            ## Update U.W.lint[]
            U.W.lint <- rbind(U.W.lint, addvalues[1:numiter,, drop = FALSE])
            ## Sort again
            U.W.lint <- U.W.lint[order(U.W.lint[, 1]), , drop = FALSE]
            ## Update length (=nrow) of 'U.W.lint'
            numObs <- numObs + numiter
            u.next
         }
      }
      ## For readability
      u.left  <- if(uLuR[1] <= ZERO) 0 else uLuR[1]
      u.right <- if(uLuR[2] >= ONE)  1 else uLuR[2]
      
      ## 2.3 Estimate integral in the three regions ##########################
      
      ## 2.3.1 Estimate log of integral over (0, u.left)
      ldens.left <- if(is.na(ldens.left)) {
         ## 'ldens.left' not set yet => peek is *not* in (0, u.left) => Riemann sums
         if(u.left == 0) {
            -Inf
         } else if(u.left == 1) {
            ## => Special case where no obs > threshold
            ## => Use all obs and Riemann for this region
            weights <- abs(c(U.W.lint[1, 1], U.W.lint[2:numObs, 1] -
                                U.W.lint[1:(numObs-1), 1]))
            logsumexp(as.matrix(log(weights) + 
                                   (c(-Inf, U.W.lint[1:(numObs-1), 3]) + 
                                       U.W.lint[1:numObs, 3])/2, ncol = 1))
         } else {
            ## 0 < u.left < 1 => Find obs in (0, u.left)
            u_sml <- (U.W.lint[, 1] <= u.left)
            sum_u_sml <- sum(u_sml)
            if(sum_u_sml > 1) {
               ## Case 1: We have >1 observations in (0, u.left)
               last_sml <- which(u_sml)[sum_u_sml]
               weights <- abs(c(U.W.lint[1, 1], 
                                U.W.lint[2:last_sml, 1] -
                                   U.W.lint[1:(last_sml-1), 1]))
               res <- logsumexp(as.matrix(log(weights) + 
                                      (c(-Inf, U.W.lint[1:(last_sml-1), 3]) + 
                                          U.W.lint[1:last_sml, 3])/2, ncol = 1))
               if(is.nan(res)) -Inf else res
            } else {
               ## Case 2: No observations in (0, u.left) 
               log(u.left) + l.tol.int.lower - log(2)
            }
         }
      } else {
         ldens.left
      }
      
      ## 2.3.2 Estimate log of integral over (u.right, 1)
      ldens.right <- if(is.na(ldens.right)) {
         ## 'ldens.right' not set yet => peek is *not* in (u.right, 1) => Riemann sums
         if(u.right == 1) {
            -Inf
         } else if(u.right == 0) {
            ## => Special case where no obs > threshold
            ## => Use all obs and Riemann for this region
            weights <- 
               abs(c(U.W.lint[1, 1], U.W.lint[2:numObs, 1]-U.W.lint[1:(numObs-1), 1]))
            logsumexp(as.matrix(log(weights) + 
                                   (c(-Inf, U.W.lint[1:(numObs-1), 3]) + 
                                       U.W.lint[1:numObs, 3])/2, ncol = 1))
         } else {
            ## 0 < u.right < 1 => Find obs in (u.right, 1)
            u_gtr <- (U.W.lint[,1] >= u.right)
            sum_u_gtr <- sum(u_gtr)
            if(sum_u_gtr > 1) { 
               ## Case 1: We have >1 observations in (u.right, 1)
               first_gtr <- which(u_gtr)[1]
               weights <- abs(c(U.W.lint[ (first_gtr+1):numObs, 1] -
                                   U.W.lint[first_gtr:(numObs-1), 1],
                                .Machine$double.neg.eps))
               res <- logsumexp(as.matrix(log(weights) + 
                                      (c(U.W.lint[(first_gtr+1):numObs, 3], -Inf) + 
                                          U.W.lint[first_gtr:numObs, 3])/2, ncol = 1))
               if(is.nan(res)) -Inf else res
            } else {
               ## Case 2: No or only one observations in (u.right, 1)
               log1p(-u.right) + l.tol.int.lower - log(2)
            }
         }
      } else {
         ldens.right
      }
      
      ## 2.3.3 Estimate log of integral over (u.left, u.right)
      ## We use RQMC in this region, unless 'ldens.strat' was already set.
      stratlength <- u.right - u.left
      rqmc.numiter <- 0
      ldens.stratum <- if(is.na(ldens.stratum)) {
         if(stratlength > control$dnvmix.tol.stratlength) {
            ldens.obj <- 
               densmix_rqmc(qW, maha2.2 = curr.maha2.2, lconst = curr.lconst,
                            d = d, k = k, control = control, u.left = u.left, 
                            u.right = u.right, return.all = FALSE,
                            max.iter.rqmc = control$max.iter.rqmc - 
                               control$dnvmix.max.iter.rqmc.pilot)
            curr.error        <- ldens.obj$error
            rqmc.numiter <- ldens.obj$numiter
            ldens.obj$ldensities + log(stratlength)
         } else {
            log(max(stratlength, ZERO)) + l.int.theo.max - log(2)
         }
      } else {
         ldens.stratum
      }
      
      ## 2.4 Combine ###########################################################
      ldensities[ind] <- 
         logsumexp(rbind(ldens.left, ldens.right, ldens.stratum, deparse.level = 0))
      error[ind]     <- curr.error
      numiters[ind]   <- rqmc.numiter
   }
   
   ## 3. Return ################################################################
   
   list(ldensities = ldensities, error = error, numiter = numiters)
}


##' @title Density of a Multivariate Normal Variance Mixture for restricted W
##' @param qW function of one variable specifying the quantile function of W.
##' @param maha2.2 squared maha-distances divided by 2. Assumed to be *sorted*
##' @param lrdet log(sqrt(det(scale))) where 'scale' is the scale matrix of
##'        the normal variance mixture distribution.
##' @param d dimension of the Normal Variance Mixture
##' @param k power of qW^{-1} in the root (k=d => density, k=d+2 =>unnormalized
##'        conditional expectation of 1/W)
##' @param control see ?dnvmix
##' @param u.left numeric in (0,1)
##' @param u.right numeric in (0,1), > u.left. Density will be estimated
##'         conditional on W being between its 'u.left' and 'u.right' quantile.
##' @param max.iter.rqmc maximum number of iterations
##' @param return.all logical; if true, matrix (U, qW(U)) also returned.
##' @return List of three:
##'         $ldensities n-vector with computed log-density values
##'         $numiter numeric, number of iterations needed
##'         $error n-vector of error estimates for log-densities; either
##'         relative error or absolte error depending on is.na(control$dnvmix.reltol)
##'         $UsWs (B, n) matrix (U, qW(U)) where U are uniforms
##'         (only if return.all = TRUE)
##' @author Erik Hintz and Marius Hofert
densmix_rqmc <- function(qW, maha2.2, lconst, d, k = d, control, u.left = 0, 
                         u.right = 1, max.iter.rqmc, return.all)
{
   ## 1. Setup #################################################################
   dblng           <- (control$increment == "doubling")
   B               <- control$B # number of randomizations
   n               <- length(maha2.2) # sample size
   current.n       <- control$fun.eval[1] # initial sample size
   numiter         <- 0 # counter for the number of iterations
   total.fun.evals <- 0
   ZERO            <- .Machine$double.neg.eps
   ## Absolte/relative precision?
   if(is.na(control$dnvmix.reltol)) {
      ## Use absolute error
      tol <- control$dnvmix.abstol
      do.reltol <- FALSE
   } else {
      ## Use relative error
      tol <- control$dnvmix.reltol
      do.reltol <- TRUE
   }
   ## Store seed if 'sobol' is used to get the same shifts later:
   if(control$method == "sobol") {
      seeds_ <- sample(1:(1e5*B), B) # B seeds for 'sobol()'
   }
   ## Additional variables needed if the increment chosen is "dblng"
   if(dblng) {
      if(control$method == "sobol") useskip <- 0
      denom <- 1
   }
   ## Matrix to store RQMC estimates
   rqmc.estimates <- matrix(-Inf, ncol = n, nrow = B)
   ## Will be needed a lot:
   CI.factor.sqrt.B <- control$CI.factor / sqrt(B)
   ## Define trafo-function that maps u to (q,1) or (1,q) depending on 'up'
   trafo <- function(u) u.left + (u.right - u.left)*u
   ## Initialize 'max.error' to > tol so that we can enter the while loop
   max.error <- tol + 42
   ## Matrix to store U, W values => nrows = maximal number of funevals
   if(return.all) {
      max.nrow <- if(dblng) current.n*B*2^(max.iter.rqmc-1) else 
         max.iter.rqmc*B*current.n
      UsWs <- matrix(NA, ncol = 2, nrow = max.nrow)
      curr.lastrow <- 0 # will count row-index additional points are being inserted after
   }
   
   ## 2. Main loop #############################################################
   
   ## while() runs until precision abstol is reached or the number of function
   ## evaluations exceed fun.eval[2]. In each iteration, B RQMC estimates of
   ## the desired log-densities are calculated.
   while(max.error > tol & numiter < max.iter.rqmc &
         total.fun.evals < control$fun.eval[2])
   {
      ## In each randomization ...
      for(b in 1:B) {
         ## Get the point set
         U <- sort(switch(
            control$method,
            "sobol" = {
               if(dblng) {
                  qrng::sobol(n = current.n, d = 1, randomize = "digital.shift", 
                              seed = seeds_[b], skip = (useskip * current.n))
               } else {
                  qrng::sobol(n = current.n, d = 1, randomize = "digital.shift", 
                              seed = seeds_[b], skip = (numiter * current.n))
               }
            },
            "ghalton" = {
               qrng::ghalton(n = current.n, d = 1, method = "generalized")
            },
            "PRNG" = {
               runif(current.n)
            })) # sorted for later!
         
         ## Evaluate the integrand at the (next) point set
         W <- qW(U <- trafo(U)) # realizations of the mixing variable; sorted!
         ## Need to replace values < ZERO by ZERO. W is *sorted*, so check using
         ## loop instd of 'pmax' (more efficient)
         for(ind in 1:current.n) if(W[ind] < ZERO) W[ind] <- ZERO else break
         ## Update 'UsWs'
         if(return.all) {
            UsWs[(curr.lastrow + 1) : (curr.lastrow + current.n), ] <- cbind(U, W)
            curr.lastrow <- curr.lastrow + current.n
         }
         next.estimate <- .Call("eval_densmix_integrand",
                                W          = as.double(W),
                                maha2_2    = as.double(maha2.2),
                                current_n  = as.integer(current.n),
                                n          = as.integer(n),
                                d          = as.integer(d),
                                k          = as.integer(k),
                                lconst     = as.double(lconst))
         ## Update RQMC estimates
         rqmc.estimates[b,] <-
            if(dblng) {
               ## In this case both, rqmc.estimates[b,] and
               ## next.estimate depend on n.current points
               .Call("logsumexp2",
                     a = as.double(rqmc.estimates[b, ]),
                     b = as.double(next.estimate),
                     n = as.integer(n)) - log(denom)
            } else {
               ## In this case, rqmc.estimates[b,] depends on
               ## numiter * current.n points whereas next.estimate
               ## depends on current.n points
               .Call("logsumexp2",
                     a = as.double(rqmc.estimates[b,] + log(numiter)),
                     b = as.double(next.estimate),
                     n = as.integer(n)) - log(numiter + 1)
            }
      } # end for(b in 1:B)
      ## Update of various variables
      ## Double sample size and adjust denominator in averaging as well as useskip
      if(dblng) {
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
      ## Total number of function evaluations
      total.fun.evals <- total.fun.evals + B * current.n
      numiter <- numiter + 1
      ## Update error. The following is slightly faster than 'apply(..., 2, var)'
      ldensities <- logsumexp(rqmc.estimates) - log(B) # performs better than .colMeans
      vars <- .colMeans((rqmc.estimates - rep(ldensities, each = B))^2, B, n, 0)
      error <- if(!do.reltol) { # absolute error
         sqrt(vars)*CI.factor.sqrt.B
      } else { # relative error
         sqrt(vars)/abs(ldensities)*CI.factor.sqrt.B
      }
      max.error <- max(error)
   } # while()
   
   ## 3. Return ################################################################
   
   if(return.all) {
      list(ldensities = ldensities, numiter = numiter, error = error,
           UsWs = UsWs[1:curr.lastrow,])
   } else {
      list(ldensities = ldensities, numiter = numiter, error = error)
   }
}


##' @title Integration Routine of log h(u) on the unit interval
##' @param qW function of one variable specifying the quantile function of W.
##' @param maha2.2 squared maha-distances divided by 2. Assumed to be *sorted*
##' @param lconst vector of same length as 'maha2.2', see below under 'Note' 
##' @param d dimension of the Normal Variance Mixture
##' @param control see ?get_set_param()
##' @param verbose see ?dnvmix()
##' @return List of three:
##'         $ldensities n-vector with computed log-density values
##'         $error n-vector of *absolute* error estimates for log-densities
##'         $numiter n-vector of number of iterations needed
##' @author Erik Hintz and Marius Hofert
##' @note Estimates int_0^1 h(u) du where log h(u) = lconst - d/2 log qW(u) - maha2.2/qW(u).
##' For density of 'nvmix': lconst = rep(-lrdet -d/2log(2pi), length(maha2.2));
##' For density of 'gammamix': lconst = -lgamma(d/2) - d/2 log(2) + (d/2-1)log(2*maha2.2)
densmix_ <- function(qW, maha2.2, lconst, d, control, verbose)
{
   ## Absolte/relative precision?
   tol <- if(is.na(control$dnvmix.reltol)) { # if 'reltol = NA' use absolute precision
      do.reltol <- FALSE
      control$dnvmix.abstol
   } else { # otherwise use relative precision (default)
      do.reltol <- TRUE
      control$dnvmix.reltol
   }
   ## Call RQMC procedure without any stratification
   rqmc.obj <- densmix_rqmc(qW, maha2.2 = maha2.2, lconst = lconst, d = d,
                            control = control, u.left = 0, u.right = 1,
                            max.iter.rqmc = control$dnvmix.max.iter.rqmc.pilot,
                            return.all = TRUE)
   ## Extract results
   ldens   <- rqmc.obj$ldensities
   numiter <- rep(rqmc.obj$numiter, length(maha2.2))
   error   <- rqmc.obj$error
   if(any(error > tol)) {
      ## Accuracy not reached for at least one 'maha2.2' value
      ## => Use adaptive approach 
      if(control$dnvmix.doAdapt) {
         notRchd <- which(error > tol)
         rqmc.obj <- densmix_adaptrqmc(
            qW, maha2.2 = maha2.2[notRchd], lconst = lconst[notRchd], d = d, 
            UsWs = rqmc.obj$UsWs,control = control)
         ldens[notRchd]    <- rqmc.obj$ldensities
         numiter[notRchd]  <- numiter[notRchd] + rqmc.obj$numiter
         error[notRchd]    <- rqmc.obj$error
      }
      ## Handle warnings:
      if(as.logical(verbose)) { # verbose can be passed by 'fitnvmix' => between 0 and 3
         if(any(is.na(error))) {
            ## At least one error is NA
            warning("Estimation unreliable, corresponding error estimate NA")
         }
         whichNA <- which(is.na(error))
         if(any(error[setdiff(1:length(error), whichNA)] > tol)) # 'setdiff' needed if 'whichNA' is empty
            warning("Tolerance not reached for all inputs; consider increasing 'max.iter.rqmc' in the 'control' argument.")
      }
   }
   ## Compute relative/absolute error for return 
   if(do.reltol){
      relerror <- error
      abserror <- error * abs(ldens)
   } else {
      abserror <- error
      relerror <- error/abs(ldens)
   } 
   ## Return
   list(ldensities = ldens, numiter = numiter, abserror = abserror, 
        relerror = relerror)
}


##' @title Density of a Multivariate Normal Variance Mixture
##' @param x (n, d)-matrix of evaluation points
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
##' @param factor Cholesky factor (lower triangular matrix) of 'scale';
##'        important here so that det(scale) is computed correctly!
##' @param control list; see ?get_set_param()
##' @param log logical indicating whether the logarithmic density is to be computed
##' @param verbose logical indicating whether warnings shall be thrown.
##' @param ... additional arguments passed to the underlying mixing distribution
##' @return n-vector with computed density values and attributes 'error'
##'         (error estimate) and 'numiter' (number of while-loop iterations)
##' @author Erik Hintz and Marius Hofert
dnvmix <- function(x, qmix, loc = rep(0, d), scale = diag(d),
                   factor = NULL, # needs to be lower triangular!
                   control = list(), log = FALSE, verbose = TRUE, ...)
{
   
   ## 1 Setup ##################################################################
   
   if(!is.matrix(x)) x <- rbind(x)
   d <- ncol(x) # dimension
   if(!is.matrix(scale)) scale <- as.matrix(scale)
   stopifnot(length(loc) == d, dim(scale) == c(d, d))
   verbose <- as.logical(verbose)
   ## Deal with algorithm parameters
   control <- get_set_param(control)
   ## If 'factor' is not provided, determine it here as a *lower* triangular matrix
   if(is.null(factor)) factor <- t(chol(scale)) # lower triangular
   ## Deal with 'qmix' 
   mix_list      <- get_mix_(qmix = qmix, callingfun = "dnvmix", ... ) 
   qW            <- mix_list[[1]] # function(u)
   special.mix   <- mix_list[[2]] # string or NA
   ## Build result object (log-density)
   lres <- rep(-Inf, (n <- nrow(x))) # n-vector of results
   abserror <- rep(NA, n)
   relerror <- rep(NA, n)
   notNA <- rowSums(is.na(x)) == 0
   lres[!notNA] <- NA
   x <- x[notNA,, drop = FALSE] # non-missing data (rows)
   
   ## 2 Actual computation ####################################################
   
   ## Recall that 'factor' is *lower triangular*. For short, let 'factor' = L
   ## Solve L * z = x_i - mu for z, so z = L^{-1} * (x_i - mu)   (d vector)
   ## => z^2 (=> componentwise) = z^T z = (x_i - mu)^T * (L^{-1})^T L^{-1} (x_i - mu)
   ##                           = z^T z = (x_i - mu)^T * (L L^T )^{-1} (x_i - mu)
   ##                           = (x_i - mu)^T * scale^{-1} * (x_i - mu)
   ##                           = quadratic form
   ## Now do this for *all* x_i simultaneously using that L is lower triangular:
   ## Forwardsolve: "right-hand sides" of equation must be in the *columns*, thus t(x)
   z <- forwardsolve(factor, t(x) - loc, transpose = FALSE)
   maha2 <- colSums(z^2) # = sum(z^T z); n-vector of squared Mahalanobis distances from x to mu w.r.t. scale
   ## Note: could probably be done with mahalanobis() but unclear how we would
   ##       get det(scale) then.
   ## log(sqrt(det(scale))) = log(det(scale))/2 = log(det(R^T R))/2 = log(det(R)^2)/2
   ##                       = log(prod(diag(R))) = sum(log(diag(R)))
   lrdet <- sum(log(diag(factor)))
   if(!is.finite(lrdet)) stop(paste("Density not defined for singular 'scale' "))
   ## Counter
   numiter <- 0 # initialize counter 
   ## Deal with the different distributions
   if(!is.na(special.mix)) {
      lres[notNA] <- switch(
         special.mix,
         "inverse.gamma" = {
            df <- mix_list$param
            lgamma((df + d) / 2) - lgamma(df/2) - (d/2) * log(df*pi) - lrdet - 
               (df+d)/2 * log1p(maha2/df)
         },
         "constant" = {
            -(d/2) * log(2 * pi) - lrdet - maha2/2
         },
         "pareto" = {
            alpha <- mix_list$param
            log(alpha) - d/2*log(2*pi) - lrdet - (alpha+d/2)*log(maha2/2) +
               pgamma(maha2/2, scale = 1, shape = alpha+d/2, log.p = TRUE) + 
               lgamma(alpha+d/2)
         })
      if(!log) lres <- exp(lres) # already exponentiate
      abserror <- rep(0, length(maha2))
      relerror <- rep(0, length(maha2)) # no error here 
   } else {
      ## General case of a multivariate normal variance mixture (=> RQMC)
      ## Prepare inputs for densmix_()
      ## Sort maha-distance and divide by 2; store ordering to recover original
      ## ordering later
      ordering.maha <- order(maha2)
      maha2.2 <- maha2[ordering.maha]/2
      ## Define log-constant for the integration
      lconst <- rep(-lrdet - d/2*log(2*pi) , length(maha2.2))
      ## Call internal densmix_ (which itself calls C-Code and handles warnings)
      ests <- densmix_(qW, maha2.2 = maha2.2, lconst = lconst,
                       d = d, control = control, verbose = verbose)
      ## Grab results
      lres[notNA] <- ests$ldensities[order(ordering.maha)]
      relerror[notNA] <- ests$relerror[order(ordering.maha)]
      abserror[notNA] <- ests$abserror[order(ordering.maha)]
      ## Correct results and error if 'log = FALSE'
      if(!log){
         lres[notNA] <- exp(lres[notNA])
         ## CI for mu: exp(logmu_hat +/- abserr(logmu_hat))) = (lower, upper)
         ## => compute max. error on mu_hat as max( (upper - mu), (mu - lower) ) 
         relerror[notNA] <- max( (exp(abserror[notNA]) - 1), (1 - exp(-abserror[notNA])) )
         abserror[notNA] <- lres[notNA] * relerror[notNA] 
      }
      numiter <- ests$numiter
   }
   
   ## 3. Return ################################################################
   
   ## Note that 'lres' was exponentiated already if necessary.
   attr(lres, "abs. error") <- abserror 
   attr(lres, "rel. error") <- relerror 
   attr(lres, "numiter") <- numiter
   lres
}
