#### Started from "Version 8" (Ver8.R)

## Initialize the weight -- Find first free weight
initAlgo <- function(mu, lB, uB ) {
    ## New-ordered return, lB, uB with decreasing return
    w <- c()
    index.new <- order(mu,decreasing = TRUE) # new order with decreasing return
    lB.new <- lB[index.new]
    uB.new <- uB[index.new]
    ## free weight - starting solution
    i.new <- 0
    w.new <- lB.new # initially
    while(sum(w.new) < 1) {
        i.new <- i.new + 1
        w.new[i.new] <- uB.new[i.new]
    }
    w.new[i.new] <- 1 - sum(w.new[-i.new])
    w[index.new] <- w.new                # back to original order
    i <- index.new[i.new]
    ## return the index of first free asset and vector w :
    list(index = i, weights = w)
}

## getMatrices -----------------------------------------------------------------
getMatrices <- function(mu, covar, w, f) {
    ## Slice covarF,covarFB,covarB, muF,muB, wF,wB
    covarF <- covar[f,f , drop=FALSE]
    muF <- mu[f]
    b <- seq_along(mu)[-f]
    covarFB <- covar[f,b , drop=FALSE]
    wB <- w[b]
    list(covarF = covarF, covarFB = covarFB, muF = muF, wB = wB)
}

computeInv <- function(get) {
    solve(get$covarF, cbind(1, get$muF, get$covarFB %*% get$wB, deparse.level = 0L))
}

## computeW -----------------------------------------------------------------
computeW <- function(lam, inv, wB) {
    ## w2 <- inv[,1]; w3 <- inv[,2]; w1 <- inv[,3]
    inv.s <- colSums(inv) # g1 <- inv.s[2]; g2 <- inv.s[1]; g4 <- inv.s[3]
    ## 1) compute gamma
    g <- (-lam * inv.s[2] + (1- sum(wB) + inv.s[3]))/inv.s[1]
    ## 2) compute free weights
    list(wF = - inv[,3] + g * inv[,1] + lam * inv[,2],
         gamma = g)
}

## computeLambda --------------------------------------------------------------
computeLambda <- function(wB, inv, i, bi.input) {
    inv.s <- colSums(inv)
    ## c1 <- inv.s[1]; l2 <- inv.s[3]; c2i <- inv[i,2];
    ## c3 <- inv.s[2]; c4i <- inv[i, 1]; l1 <- sum(wB)
    c1 <- inv.s[1]
    if(length(bi.input) == 1L) {                    # 1.bound to free
        c4i <- inv[i, 1]
        Ci <- - c1 * inv[i, 2] + inv.s[2] * c4i
        ## return lambda :
        if(Ci == 0)
            0
        else
            ((1- sum(wB) + inv.s[3])* c4i- c1 * (bi.input + inv[i, 3]))/Ci
    }
    else {                                          # 2.free to bound
        c4i <- inv[, 1]
        Ci <- - c1 * inv[, 2] + inv.s[2] * c4i
        bi          <- bi.input[i, 1]         # bi.lB
        bi[Ci > 0]  <- bi.input[i[Ci > 0], 2] # bi.uB
        bi[Ci == 0] <- 0
        ## return lambda and boundary :
        list(lambda = ((1- sum(wB) + inv.s[3]) * c4i- c1 *(bi + inv[, 3]))/Ci,
             bi = bi)
    }
}

MS <- function(weights_set, mu, covar) {
    Sig2 <- colSums(weights_set *(covar %*% weights_set) )
    cbind(Sig = sqrt(Sig2), Mu = as.vector(t(weights_set) %*% mu))
}


CLA <- function(mu, covar, lB, uB, check.cov = TRUE, check.f = TRUE, tol.lambda = 1e-7,
                give.MS = TRUE, keep.names = TRUE, trace = 0) {
    ## minimal argument checks
    cl <- match.call()
    n <- length(mu)
    if(length(lB) == 1) lB <- rep.int(lB, n)
    if(length(uB) == 1) uB <- rep.int(uB, n)
    stopifnot(is.numeric(mu), is.matrix(covar), identical(dim(covar), c(n,n)),
              is.numeric(lB), length(lB) == n,
              is.numeric(uB), length(uB) == n, lB <= uB)# and in [0,1]
    if(check.cov) {
        ## if the covar matrix is *not* positive (semi)definite, CLA()
        ## may produce a "solution" which does not fulfill desired properties

        ## 1) Check  kappa() or rcond() ... if they are "large", check the eigen values

        " __ FIXME __ "

        ## 2)
        ev <- eigen(covar, only.values=TRUE)$values

        if(any(ev < 0))
            warning("covariance matrix 'covar' has negative eigenvalues")
    }
    ## Compute the turning points, free sets and weights
    ans <- initAlgo(mu, lB, uB)
    f <- ans$index
    w <- ans$weights
    ## initialize result parts
    lambdas <- gammas  <- numeric()
    weights_set <- array(dim = c(n,0L))
    eLambdas <- free_indices <- list()
    lam <- 1 # set non-zero lam
    while (lam > 0 && (nf <- length(f)) <= length(mu)) {
      if(trace) cat(sprintf("while(lam = %g > 0 && |f|=%d <= |mu|=%d)\n", lam, nf, length(mu)))
      ## 1) case a): Bound one free weight  F -> B
      l_in <- 0
      if(nf > 1L) {
        compl <- computeLambda(wB = w[-f], inv = inv, # inv from last step k (k > 1)
                               i = f, bi.input = cbind(lB, uB))
        if(trace >= 2)  { cat(" case a) : computeLambda(): "); str(compl) }
        lam_in <- compl$lambda
        bi     <- compl$bi
        k <- which.max(lam_in)
        i_in  <-      f[k]
        bi_in <-     bi[k]
        l_in  <- lam_in[k]
        if(trace >= 2)  cat(sprintf("  [a) cont.]: k = which.max(lam_in) = %d; l_in = %g\n",
                                    k, l_in))
      }

      ## 2) case b): Free one bounded weight  B -> F
      b <- seq_along(mu)[-f]
      inv_list <- lapply(b, function(bi) {
        get_i <- getMatrices(mu, covar, w, c(f,bi))
        computeInv(get_i)
      })
      if(trace >= 2)  { cat(sprintf(" case b) : \"B -> F\": b = (%s);    inv_list:\n",
                                    paste(b,collapse=", "))) ; str(inv_list) }
      if(nf < length(mu)) { # still have free weights
          fi <- nf + 1L
          lam_out <- sapply(seq_along(b), function(i)
              computeLambda(wB = w[b[-i]], inv = inv_list[[i]],
                            i = fi, bi.input = w[b[i]]))
          if(trace) {
              cat(sprintf("|f| < |mu|: computeLambda() => lam_out[1:%d]%s",
                          length(b), if(trace >= 2) "= " else ""))
              if(trace >= 2) print(lam_out)
          }
          if (length(lambdas) && !all(sml <- lam_out < lam*(1-tol.lambda))) {
              tol.l <- tol.lambda
              while((!any(sml)) && 2*tol.l >= .Machine$double.neg.eps) # empty
                  ## extreme: new lam_out are *not* smaller than lam(1-eps)
                  sml <- lam_out < lam*(1 - (tol.l <- tol.l/2))
              if(trace) cat("  new 'sml' case: which(sml) = ",
                            if(any(sml)) paste(which(sml), collapse=", ")
                            else "{} (i.e. empty)", "\n")
              lam_out <- lam_out[sml]
              b       <- b      [sml]
              inv_list <- inv_list[sml]
          }
          if((hasLam <- length(lam_out) > 0)) {
              k <- which.max(lam_out)
              if(trace) cat("   --> new k = which.max(lam_out): ", k, "\n")
              i_out <- b      [k] # one only !
              l_out <- lam_out[k]
              inv_out <- inv_list[[k]]
          } else { # 'empty' --- should not happen typically, but see 'mc3' ex.
              if(check.f)
                  warning("Had free weights but could not improve solution")
              l_out  <- -Inf
          }
      } else { ## length(f) == length(mu)  <==>  |b| = 0
          hasLam <- FALSE
          l_out  <- -Inf
      }
      ## 3) decide lambda
      lam <- max(l_in, l_out)
      if(trace) cat(sprintf("  l_{in,out}=(%g,%g) => new candidate lam=%g\n",
                            l_in, l_out, lam))
      if(lam > 0) { # remove i_in from f; or add i_out into f
        if(l_in > l_out) {
          w[i_in] <- bi_in  # set value at the correct boundary
          f <- f[f != i_in]
          getM <- getMatrices(mu, covar, w, f)
          inv <- computeInv(getM)
        }
        else {
          f <- c(f, i_out)
          inv <- inv_out
        }
      }
      else{ # 4) if max(l_in, l_out) <= 0, "stop" when at the min var solution!
        lam <- 0
        ## muF = 0 not necessary, get1 replaced by getM (ie getM from previous step)
      }
      compW <- computeW(lam, inv = inv, wB = w[-f])
      g    <- compW$gamma
      w[f] <- compW$wF[seq_along(f)]

      lambdas <- c(lambdas, lam)
      if(!hasLam)
          eLambdas <- c(eLambdas, lam)
      gammas <- c(gammas, g)
      weights_set <- cbind(weights_set, w, deparse.level = 0L) # store solution
      free_indices <- c(free_indices, list(sort(f)))
    }# end While

    if(keep.names) rownames(weights_set) <- names(mu)
    structure(class = "CLA",
              list(weights_set = weights_set,
                   free_indices = free_indices,
                   gammas = gammas, lambdas = lambdas,
                   emptyLambdas =  eLambdas,
                   MS_weights = if(give.MS)
                           MS(weights_set = weights_set, mu = mu, covar = covar),
                   call = cl))
}


## print method
print.CLA <- function(x, ...) {
    cat("Call: ", paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "\n", sep = "")
    wts <- x$weights_set
    n <- nrow(wts)
    nT <- ncol(wts)
    cat(gettextf("Critical Line Algorithm for n = %d assets, ", n),
        gettextf("resulting in  %d  turning points", nT),"\n", sep="")
    ## For now; we can do better later:
    cat("Overview of result parts:\n")
    utils::str(x[1:5], max.level = 1, give.attr=FALSE)
    ## TODO better, e.g., summarizing the number "active assets"
    ## and those with non-0 weights -- only if lower bounds were (mostly) 0
    invisible(x)
}

## TODO: see ../inst/issues/nonPD/many_NA_data__non-PD-covmat.R : cbind(...)
if(FALSE)
summary.CLA <- function(object, ...) {
}


## As basically from  .../YanhaoShi/R/Functions/Plot.R :
MS_plot <- function(ms, type = "o",
                    main = "Efficient Frontier",
                    xlab = expression(sigma(w)),
                    ylab = expression(mu(w)),
                    col = adjustcolor("blue", alpha.f = 0.5),
                    pch = 16, ...) {
    ## list of weights_set, legend...
    stopifnot(is.matrix(ms), ncol(ms) == 2)
    plot(ms[,"Sig"], ms[,"Mu"], type=type, pch=pch, col=col,
         xlab = xlab, ylab=ylab, main=main, ...)
}


## FIXME: --> see also in ../man/plot.CLA.Rd
## -----
## 0) Use findMu() and findSig() to draw the lines *between*
## 1) Learn from Tobias Setz to plot the lower part of the feasible region
## 2) Better title, using 'call'
## 3) mark some critical points particularly
## 4) give information about the *number* critical points / weights sets
## 5) consider using a  'add = FALSE' argument and then use 'lines()'
## 6) "label" turning points by #{assets}, or plot these additionally, or use axis 4
plot.CLA <- function(x, type = "o", main = "Efficient Frontier",
                    xlab = expression(sigma(w)),
                    ylab = expression(mu(w)),
                    col = adjustcolor("blue", alpha.f = 0.5),
                    pch = 16, ...) {
    stopifnot(is.matrix(ms <- x$MS_weights))
    plot(ms[,"Sig"], ms[,"Mu"], type=type, pch=pch, col=col,
         xlab=xlab, ylab=ylab, main=main, ...)
}

## TODO 2nd plot: --> ../inst/issues/nonPD/many_NA_data__non-PD-covmat.R
