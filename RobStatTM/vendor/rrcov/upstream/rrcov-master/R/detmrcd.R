## @title Compute the Minimum Regularized Covariance Determinant (MRCD) estimator
## @references Paper available at: http://dx.doi.org/10.2139/ssrn.2905259.
## @param x a numerical matrix. The columns represent variables, and rows represent observations.
## @param alpha the proportion of the contamination (between 0.5 and 1)
## @param h the size of the subset (between ceiling(n/2) and n)
## @param initHsets NULL or a K x h integer matrix of initial subsets of observations
##     of size h (specified by the indices in 1:n). If provided, then the initial
##     shape estimates are not calculated.
## @param save.hsets
## @param maxcsteps maximum number of generalized C-steps for each initial subset (default 200)
## @param maxcond maximum condition number allowed (see step 3.4 in algorithm 1) (default 50)
## @param minscale minimum scale allowed (default 0.001)
## @param target = c("identity", "equicorrelation"). Structure of the robust
##     positive definite target matrix: (default) "identity": target matrix is
##     diagonal matrix with robustly estimated univariate scales on the diagonal or
##     "equicorrelation": non-diagonal target matrix that incorporates an
##     equicorrelation structure (see (17) in paper)
## @param trace
## @return A list with the following elements:
## \describe{
## \item{icov}{inverse of the covariance matrix}
## \item{rho}{regularization parameter}
## \item{target}{the target matrix used}
## }
##
.detmrcd <- function(x, h=NULL, alpha=.75, rho=NULL,
            maxcond=50, minscale=0.001, target=0, maxcsteps = 200,
            hsets.init=NULL, save.hsets=missing(hsets.init), full.h = save.hsets,
            trace=FALSE)
{
## NOTES: VT
##
##  - use r6pack from robustbase
##  - X check subset #5 - something is wrong there
##  - X mX is not back transformed to compensate the SVD rescaling
##      ==> the distances returned by MWRCD are different from
##      the distances computed outside of MWRCD.
##
##  - X MRCD works with the transposed matrix X - it shoould be transposed inside the function
##  - Help of CovMcd - remove data.matrix from the example
##  - X What is doing the parameter 'bc'? Can be omitted - bc removed
##  - X parameter 'initrho' is never used: remove it
##  - X scfactor(alpha, p) is equivalent to robustbase:::.MCDcons(p, alpha) - replaced
##  - X condnumber () is never used - remove it - removed
##  - X robustbase:::r6pack() will not work for p > n, because there is a check svd$rank < p - comment out the check
##  - X using kendal's tau cor.fk from pcaPP - replace by spearman. Later we could move cor.fk() from pcaPP into robustbase
##  - X No SVD if the target matrix is the Identity; if equicorrelation, the eigenvectors are a helmert matrix
##----------------------------------------------------------------

## @title Robust Distance based observation orderings based on robust "Six pack"
## @param x  n x p data matrix
## @param h  integer
## @param full.h full (length n) ordering or only the first h?
## @param scaled is 'x' is already scaled?  otherwise, apply doScale(x, median, scalefn)
## @param scalefn function to compute a robust univariate scale.
## @return a h' x 6 matrix of indices from 1:n; if(full.h) h' = n else h' = h
r6pack <- function(x, h, full.h, adjust.eignevalues=TRUE, scaled=TRUE, scalefn=Qn)
{
    ## As the considered initial estimators Sk may have very
    ## inaccurate eigenvalues, we try to 'improve' them by applying
    ## a transformation similar to that used in the OGK algorithm.
    ##
    ## After that compute the corresponding distances, order them and
    ## return the indices
    initset <- function(data, scalefn, P, h)
    {
       stopifnot(length(d <- dim(data)) == 2, length(h) == 1, h >= 1)
       n <- d[1]
       stopifnot(h <= n)
       lambda <- doScale(data %*% P, center=median, scale=scalefn)$scale
       sqrtcov    <- P %*% (lambda * t(P)) ## == P %*% diag(lambda) %*% t(P)
       sqrtinvcov <- P %*% (t(P) / lambda) ## == P %*% diag(1/lambda) %*% t(P)
	   estloc <- colMedians(data %*% sqrtinvcov) %*% sqrtcov
       centeredx <- (data - rep(estloc, each=n)) %*% P
	   sort.list(mahalanobisD(centeredx, FALSE, lambda))[1:h]# , partial = 1:h
    }

    ##
    ##  Compute the raw OGK estimator. For m(.) and s(.) (robust
    ##  univariate estimators of location and scale) use the median
    ##  and Qn for reasons of simplicity (no choice of tuning parameters)
    ##  and to be consistent with the other components of DetMCD.
    ##
    ogkscatter <- function(Y, scalefn, only.P = TRUE)
    {
        stopifnot(length(p <- ncol(Y)) == 1, p >= 1)
        U <- diag(p)

        for(i in seq_len(p)[-1L]) {# i = 2:p
            sYi <- Y[,i]
            ii <- seq_len(i - 1L)
            for(j in ii) {
                sYj <- Y[,j]
                U[i,j] <- (scalefn(sYi + sYj)^2 - scalefn(sYi - sYj)^2) / 4
            }
            ## also set the upper triangle
            U[ii,i] <- U[i,ii]
        }

        ## now done above: U <- lower.tri(U) * U + t(U)    #    U <- tril(U, -1) + t(U)
        P <- eigen(U, symmetric=TRUE)$vectors
    	if(only.P)
    	    return(P)

        ## else :
        Z <- Y %*% t(P)
        sigz <- apply(Z, 2, scalefn)
        lambda <- diag(sigz^2)

        list(P=P, lambda=lambda)
    }

    stopifnot(length(dx <- dim(x)) == 2)
    n <- dx[1]
    p <- dx[2]

    ## If scalefn is missing or is NULL, use Qn for smaller data sets (n < 1000)
    ## and tau-scale of Yohai and Zamar (1988) otherwise.
    ## scalefn <- robustbase:::robScalefn(scalefn, n)

    ## If the data was not scaled already (scaled=FALSE), center and scale using
    ## the median and the provided function 'scalefn'.
    if(!scaled) { ## Center and scale the data to (0, 1) - robustly
        x <- doScale(x, center=median, scale=scalefn)$x
    }

    nsets <- 6
    hsets <- matrix(integer(), h, nsets)

    ## Determine 6 initial estimates (ordering of obs)
    ## 1. Hyperbolic tangent of standardized data
    y1 <- tanh(x)
    R1 <- cor(y1)
    P <- eigen(R1, symmetric=TRUE)$vectors
    hsets[,1] <- initset(x, scalefn=scalefn, P=P, h=h)

    ## 2. Spearmann correlation matrix
    R2 <- cor(x, method="spearman")
    P <- eigen(R2, symmetric=TRUE)$vectors
    hsets[,2] <- initset(x, scalefn=scalefn, P=P, h=h)

    ## 3. Tukey normal scores
    y3 <- qnorm((apply(x, 2L, rank) - 1/3)/(n + 1/3))
    R3 <- cor(y3, use = "complete.obs")
    P <- eigen(R3, symmetric=TRUE)$vectors
    hsets[,3] <- initset(x, scalefn=scalefn, P=P, h=h)

    ## 4. Spatial sign covariance matrix
    znorm <- sqrt(rowSums(x^2))
    ii <- znorm > .Machine$double.eps
    x.nrmd <- x
    x.nrmd[ii,] <- x[ii, ] / znorm[ii]
    SCM <- crossprod(x.nrmd)# / (n-1) not needed for e.vectors
    P <- eigen(SCM, symmetric=TRUE)$vectors
    hsets[,4] <- initset(x, scalefn=scalefn, P=P, h=h)

    ## 5. BACON
    ind5 <- order(znorm)
    half <- ceiling(n/2)
    Hinit <- ind5[1:half]
    covx <- cov(x[Hinit, , drop=FALSE])
    P <- eigen(covx, symmetric=TRUE)$vectors
    hsets[,5] <- initset(x, scalefn=scalefn, P=P, h=h)

    ## 6. Raw OGK estimate for scatter
    P <- ogkscatter(x, scalefn, only.P=TRUE)
    hsets[,6] <- initset(x, scalefn=scalefn, P=P, h=h)

    ## VT::15.11.2019
    ##  No need of the code below in MRCD.
    ##  Instead of computing md and sorting, just return hsets.
    ##
    if(!adjust.eignevalues)
        return (hsets)

    ## Now combine the six pack :
    if(full.h) hsetsN <- matrix(integer(), n, nsets)
    for(k in 1:nsets) ## sort each of the h-subsets in *increasing* Mah.distances
    {
        xk <- x[hsets[,k], , drop=FALSE]
	    svd <- classPC(xk, signflip=FALSE) # [P,T,L,r,centerX,meanvct] = classSVD(xk)

        ## VT::15.10.2018 - we do not need this check, because we are using r6pack
        ##  also for MRCD. Comment it out , maybe should move it to detmcd: FIXME
        ##  if(svd$rank < p) ## FIXME: " return("exactfit")  "
        ##      stop('More than half of the observations lie on a hyperplane.')

        score <- (x - rep(svd$center, each=n)) %*% svd$loadings
        ord <- order(mahalanobisD(score, FALSE, sqrt(abs(svd$eigenvalues))))
        if(full.h)
            hsetsN[,k] <- ord
        else hsets[,k] <- ord[1:h]
    }

    ## return
    if(full.h) hsetsN else hsets
}

# Return target correlation matrix with given structure
# input:
#   mX: the p by n matrix of the data
#   target: structure of the robust positive definite target matrix (default=1)
#           0: identity matrix
#           1: non-diagonal matrix with an equicorrelation structure (see (17) in paper)
# output:
#   the target correlation matrix
.TargetCorr <- function(mX, target=1, mindet=0)
{
    p <- dim(mX)[1]
    I <- diag(1, p)

    if(target == 0){
        R <- I
    } else if(target == 1) {
##        cortmp <- cor.fk(t(mX))                   # from pcaPP
##        cortmp <- cor(t(mX), method="kendal")     # very slow
        cortmp <- cor(t(mX), method="spearman")
        cortmp <- sin(1/2 * pi * cortmp)
        constcor <- mean(cortmp[upper.tri(cortmp, diag=FALSE)])

        # KB: add bound to ensure positive definiteness; see paper page 7, below (17)
        if(constcor <= min(c(0, (-1/(p-1) + 0.01)))){
            constcor <-  min(c(0,(-1/(p-1) + 0.01)))
        }
        J <- matrix(1, p, p)
        R <- constcor * J + (1-constcor) * I
    }
    return(R)
}

eigenEQ <- function(T)
{
    rho <- T[1,2]
    d <- ncol(T)

    helmert <- matrix(0, nrow=d, ncol=d)
    helmert[, 1] <- rep(1/sqrt(d), d)
    for(j in 2:ncol(helmert))
    {
        helmert[1:(j-1), j] <- 1/sqrt(j*(j-1))
        helmert[j, j] <- -(j-1)/sqrt(j*(j-1))
    }

    out <- NULL
    out$values <- c(1 + (d-1)*rho, rep(1-rho, d-1))
    out$vectors <- helmert
    out
}

# compute the regularized covariance
# input:
#   XX, the p by n (or h) matrix of the data (not necessarily demeaned)
#   vMu: the initial mean (p-vector)
#   rho: the regularization parameter
#   mT: the target matrix
#   scfac: the scaling factor
#   bcd: diagonal matrix used for rescaling (ratio of scales of target and mcd component)
#   target: structure of the robust positive definite target matrix (default=1)
#           0: identity matrix
#           1: non-diagonal matrix with an equicorrelation structure (see (17) in paper)
#   invert: if true, gives also inverted regularized covariance
#
# output (a list):
#   rho: the regularization parameter
#   mT: the target matrix
#   cov: the covariance matrix based on subset (without regularization step)
#   rcov: the regularized covariance matrix
#   inv_rcov: the inverse of the regularized covariance matrix (if invert=True)
.RCOV <- function(XX, vMu, rho=NULL, mT, scfac, target=1, invert=FALSE)
{
    mE <- XX-vMu
    n <- dim(mE)[2]
    p <- dim(mE)[1]
    mS <- mE %*% t(mE)/n
    rcov <- rho * mT + (1-rho) * scfac * mS

    if(invert) {
        if(p > n) {
          nu <- (1-rho) * scfac
          mU <- mE/sqrt(n)
          inv_rcov <- .InvSMW(rho=rho, mT=mT, nu=nu, mU=mU)
        }else {
          inv_rcov = chol2inv(chol(rcov))
        }
        return(list(rho=rho, mT=mT, cov=mS, rcov=rcov, inv_rcov=inv_rcov))
    }
    else
        return(list(rho=rho, mT=mT, cov=mS, rcov=rcov))

}

##  Compute inverse of covariance matrix using Sherman-Morrison-Woodbury
##  identity when dimension is larger than sample size
#
# input:
#   rho: the regularization parameter
#   mT: the target matrix
#	nu: the scaling factor multiplied with (1-rho)
#   mU: the scaled data
# output:
#	  the inverse of the covariance matrix
.InvSMW <- function(rho, mT, nu, mU)
{
    p = dim(mT)[1]
    pp = dim(mU)
    vD = sqrt(diag(mT))
    imD = diag(vD^(-1))
    R = imD %*% mT %*% imD
    constcor = R[2, 1]
    I = diag(1, p)
    J = matrix(1, p, p)
    imR = 1/(1-constcor) * (I - constcor/(1 + (p-1) * constcor) * J)
    imB = (rho)^(-1) * imD %*% imR %*% imD
    Temp <- base::chol2inv(base::chol(diag(pp[2]) + nu * (t(mU) %*% (imB %*% mU))))

    return(imB - (imB%*%mU) %*% (nu * Temp) %*% (t(mU)%*%imB))
}

# Apply generalized C-steps to obtain optimal subset
# input:
#   mX: the p by T matrix of the residuals or data, not necessarily demeaned
#   rho: the regularization parameter
#   mT: the target matrix
#   target: structure of the robust positive definite target matrix (default=1)
#           0: identity matrix
#           1: non-diagonal matrix with an equicorrelation structure (see (17) in paper)
#   vMu: the initial mean (as vector)
#   mIS: the p by p matrix of the initial inverted covariance
#   h: the size of the subset OR alpha: the proportion of the contamination
#   maxcsteps: the maximal number of iteration of the C-step algorithm
#   index: the initial subset H_0
# output (a list)
#   index: the optimal h-subset
#   numit: the number of iterations
#   mu: the vector with means
#   cov: the regularized covariance estimate
#   icov: the inverse of the regularized covariance matrix
#   rho: the regularization parameter
#   mT: the target matrix
#   dist: the Mahalanobis distances using the MRCD estimates
#   scfac: the scaling factor
.cstep_mrcd <- function(mX, rho=NULL, mT=NULL, target=1, vMu=NULL, mIS=NULL, h, scfac, index=NULL, maxcsteps=50)
{
    n <- dim(mX)[2]
    p <- dim(mX)[1]

    # random choice
    if(is.null(index))
        index <- sample(1:n, h) # if no index is given we sample one...
    XX <- mX[, index] # p x h

    if(is.null(vMu))
        vMu = rowMeans(XX)

    if(is.null(mIS)){
        ret = .RCOV(XX=XX, vMu=vMu, rho=rho, mT=mT, scfac=scfac, target=target, invert=T)
        mIS = ret$inv_rcov
    }

    vdst = diag(t(mX-vMu) %*% (mIS %*% (mX-vMu))) #vdst = apply(mX-vMu,2,ftmp)
    index = sort(sort.int(vdst, index.return=T)$ix[1:h])

    iter = 1

    while(iter < maxcsteps){
        XX <- mX[,index]
        vMu <- rowMeans(XX)
        ret <- .RCOV(XX=XX, vMu=vMu, rho=rho, mT=mT, target=target, scfac=scfac, invert=T)
        mIS <- ret$inv_rcov

        vdst <- diag(t(mX-vMu) %*% (mIS %*% (mX-vMu)))
        nndex <- sort(sort.int(vdst,index.return=T)$ix[1:h])

        if(all(nndex == index))
            break

        index <- nndex
        iter <- iter+1
    }

    return(list(index=index, numit=iter, mu=vMu, cov=ret$rcov,
              icov=ret$inv_rcov, rho=ret$rho, mT=ret$mT, dist=vdst, scfac=scfac))
}

    mX <- t(x)          # we want the transposed data matrix

    ## several parametrs which we do not want toexpose to the user.
    mindet <- 0         # minimum determinant allowed for target matrix
    objective <- "geom" # objective function to determine optimal subset, see (3) in paper
                        # 'det': typically one minimizes the determinant of the sample covariance based on the subset
                        # 'geom': p-th root of determinant or standardized generalized variance (for numerical reasons)

    n <- dim(mX)[2]
    p <- dim(mX)[1]

    if(!is.null(h))          alpha <- h/n
    else if(!is.null(alpha)) h <- ceiling(alpha*n)
    else
        stop("Either 'h' (number of observations in a subset) or 'alpha' (proportion of observations) has to be supplied!")
    if(alpha < 1/2 | alpha > 1)
        stop("'alpha' must be between 0.5 and 1.0!")

    ## VT::21.04.2021 - h has to be an integer (a non-integer h will break the sprintf at the end).
    h <- as.integer(h)

    # choose objective function to determine optimal subset
    if (objective == 'det'){
        obj <- function(x) det(x)
    }else if (objective == 'geom'){
        obj <- function(x)
        {
            det(x)^(1/p)
        } #geometric mean of eigenvalues
    }

    ## 1. Standardize the p variables: compute standardized observations u_i, see (6) in paper, using median and Qn estimator
    vmx <- apply(mX, 1, median)
    vsd <- apply(mX, 1, Qn)
    vsd[vsd < minscale] <- minscale
    Dx <- diag(vsd)
    mU <- scale(t(mX), center=vmx, scale=vsd)
    mX <- t(mU)
    mT <- .TargetCorr(mX, target=target, mindet=mindet)

    ## 2. Perform singular value decomposition of target matrix and compute observations w_i
    if(target == 1){
        mTeigen <- eigenEQ(mT)
        mQ <- mTeigen$vectors
        mL <- diag(mTeigen$values)
        msqL <- diag(sqrt(mTeigen$values))
        misqL <- diag(sqrt(mTeigen$values)^(-1))
        mW <- mU %*% mQ %*% misqL
        mX <- t(mW)
    }

    mT = diag(p)

    ## 3.1-3.2 Follow Hubert et al. (2012) to obtain 6 initial
    ##      scatter matrices (if scatter matrix is not invertible,
    ##      use its regularized version)
    ## 3.3 Determine subsets with lowest Mahalanobis distance

    ## Assume that 'hsets.init' already contains h-subsets: the first h observations each
    ## VT::15.11.2019 - added adjust.eignevalues=FALSE, this will set automatically full.h=FALSE
    if(is.null(hsets.init)) {
	   hsets.init <- r6pack(x=t(mX), h=h, full.h=FALSE, adjust.eignevalues=FALSE, scaled=FALSE, scalefn=Qn)
	   dh <- dim(hsets.init)
    } else { ## user specified, (even just *one* vector):
	   if(is.vector(hsets.init)) hsets.init <- as.matrix(hsets.init)
	   dh <- dim(hsets.init)
	   if(dh[1] < h || dh[2] < 1)
	       stop("'hsets.init' must be a  h' x L  matrix (h' >= h) of observation indices")
	   if(full.h && dh[1] != n)
	       warning("'full.h' is true, but 'hsets.init' has less than n rows")
	   if(min(hsets.init) < 1 || max(hsets.init) > n)
	       stop("'hsets.init' must be in {1,2,...,n}; n = ", n)
    }

    hsets.init <- hsets.init[1:h, ]
    scfac <- robustbase::.MCDcons(p, h/n)           # for consistency with MCD

    ## 3.4 Determine smallest value of rho_i for each subset
    rho6pack <- condnr <- c()
    nsets <- ncol(hsets.init)
    if(is.null(rho)) {
        for(k in 1:nsets){
            mXsubset <- mX[ , hsets.init[, k]]
            vMusubset <- rowMeans(mXsubset)
            mE <- mXsubset-vMusubset
            mS <- mE%*%t(mE)/(h-1)

            if(all(mT == diag(p))) {
                veigen <- eigen(scfac * mS)$values
                e1 <- min(veigen)
                ep <- max(veigen)

                ##  cat("\ncase T=I: ", k, e1, ep)

                fncond <- function(rho)
                {
                    condnr <-  (rho + (1-rho) * ep) / (rho + (1-rho) * e1)
                    ##  cat("\n ...... condnr: ", condnr, condnr-maxcond, "\n")
                    return(condnr - maxcond)
                }
            } else {
                fncond <- function(rho)
                {
                    rcov <- rho*mT + (1-rho) * scfac * mS
                    temp <- eigen(rcov)$values
                    condnr <- max(temp) / min(temp)
                    return(condnr - maxcond)
                }
            }

            out <- try(uniroot(f=fncond, lower=0.00001, upper=0.99), silent=TRUE)

            ## VT::11.08.2022: fix error "Found if() conditions comparing class() to string"
            ##  if(class(out) != "try-error") {
            if(!is(out, "try-error")) {
                rho6pack[k] <- out$root
                ##  cat("\nOK: ", k, out$root, "\n")
            }else {
                grid <- c(0.000001, seq(0.001, 0.99, by=0.001), 0.999999)
                if(all(mT == diag(p))) {
                    objgrid <- abs(fncond(grid))
                    irho <- min(grid[objgrid == min(objgrid)])
                }else {
                    objgrid <- abs(apply(as.matrix(grid), 1, "fncond"))
                    irho <- min(grid[objgrid == min(objgrid)])
                }
                rho6pack[k] <- irho
                ##  cat("\nNOT OK: ", k, irho, "\n")

            }
        }

        ## 3.5 Set rho as max of the rho_i's obtained for each subset in previous step
        cutoffrho <- max(c(0.1, median(rho6pack)))
        rho <- max(rho6pack[rho6pack <= cutoffrho])

        if(trace) {
            cat("\nSet rho as max of the rho_i obtained for each subset in previous step.")
            cat("\nrho, cutoffrho=", rho, cutoffrho, "\n")
            print(rho6pack)
        }

        Vselection <- seq(1, nsets)
        Vselection[rho6pack > cutoffrho] = NA
        if(sum(!is.na(Vselection)) == 0){
            stop("None of the initial subsets is well-conditioned")
        }

        initV <- min(Vselection, na.rm=TRUE)
        setsV <- Vselection[!is.na(Vselection)]
        setsV <- setsV[-1]
    }else{
        setsV <- 1:ncol(hsets.init)
        initV <- 1
    }

    ## 3.6 For each of the six initial subsets, repeat the generalized
    ##  C-steps (from Theorem 1) until convergence
    ##
    ## 3.7 Choose final subset that has lowest determinant among the
    ##  ones obtained from the six initial subsets
    hset.csteps <- integer(nsets)
    ret <- .cstep_mrcd(mX=mX, rho=rho, mT=mT, target=target, h=h, scfac=scfac, index=hsets.init[, initV], maxcsteps=maxcsteps)
    objret <- obj(ret$cov)
    hindex <- ret$index
    best6pack <- initV
    for(k in setsV){
        if(trace) {
            if(trace >= 2)
                cat(sprintf("H-subset %d = observations c(%s):\n-----------\n",
                    k, paste(hsets.init[1:h, k], collapse=", ")))
             else
                cat(sprintf("H-subset %d: ", k))
       }
       tmp <- .cstep_mrcd(mX=mX, rho=rho, mT=mT, target=target, h=h, scfac=scfac, index=hsets.init[,k], maxcsteps=maxcsteps)
        objtmp <- obj(tmp$cov)
        hset.csteps[k] <- tmp$numit
        if(trace)
            cat(sprintf("%3d csteps, obj=log(det|.|)=%g", k, objtmp))
        if(objtmp < objret){
            if(trace)
                cat(" = new optim.\n")
            ret <- tmp
            objret <- objtmp
            hindex <- tmp$index
            best6pack <- k
        } else if(objtmp == objret)     # store as well
            best6pack <- c(best6pack, k)
        else
            if(trace) cat("\n")
    }

    c_alpha <- ret$scfac #scaling factor
    mE <- mX[, hindex] - ret$mu
    weightedScov <-  mE %*% t(mE)/(h-1)
    D  <- c_alpha * diag(1, p)

    ## MRCD estimates of the standardized data W (inner part of (12) in paper)
    MRCDmu = rowMeans(mX[,hindex])
    MRCDcov = rho*mT + (1-rho) * c_alpha * weightedScov

    ## Computing inverse of scaled covariance matrix, using SMW identity
    ##      if data is fat (inner part of (14) and (15) in paper).
    if(p > n & target <= 1){    # !!!! formula InvSMW  is used when T is equicorrelation
        nu <- (1-rho) * c_alpha
        mU <- mE/sqrt(h-1)
        iMRCDcov <- .InvSMW(rho=rho, mT=mT, nu=nu, mU=mU)
    }else
        iMRCDcov <- chol2inv(chol(MRCDcov))

    ## Backtransforming the rescaling steps that we applied in
    ##  the beginning (outer part of (11) and (12) in paper)

    # transformations due to SVD on target matrix
    if(target == 1) {
        ##VT::12.12 - corrected the restoration of X after SVD transformation
        mX <- t(t(mX) %*% msqL %*% t(mQ))
##        mX <- t(t(mX) %*% mQ %*% msqL)
        MRCDmu <- mQ %*% msqL %*% MRCDmu
        MRCDcov <- mQ %*% msqL %*% MRCDcov %*% msqL %*% t(mQ)
        iMRCDcov <- mQ %*% misqL %*% iMRCDcov %*% misqL %*% t(mQ)
        mT <- mQ %*% msqL %*% mT %*% msqL %*% t(mQ)
    }

    # transformations due to rescaling median and Qn
    mX <- t(t(mX) %*% Dx) + vmx
    MRCDmu <- Dx %*% MRCDmu + vmx
    MRCDcov <- Dx %*% MRCDcov %*% Dx
    mT <- Dx %*% mT %*% Dx
    iDx <- diag(1/diag(Dx))
    iMRCDcov <- iDx %*% iMRCDcov %*% iDx

    ## Compute the Mahalanobis distances based on MRCD estimates
    dist <- mahalanobis(t(mX), center=MRCDmu, cov=iMRCDcov, inverted=TRUE)
    objret <- determinant(MRCDcov)$modulus[1]

    ret <- list(alpha=alpha, h=h,
              initmean=as.numeric(MRCDmu),
              initcovariance=MRCDcov,
              icov=iMRCDcov,
              rho=rho,
              best=hindex,
              mcdestimate=objret,
              mah=dist,
              target=mT,
              iBest = best6pack,
              n.csteps=hset.csteps,
              initHsets=if(save.hsets) hsets.init,
              calpha=c_alpha
            )

    return (ret)
}
