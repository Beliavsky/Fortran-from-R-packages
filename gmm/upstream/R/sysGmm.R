#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/


five <- function(g, h, commonCoef = FALSE, data = NULL)
    {
        res  <- sysGmm(g, h, wmatrix = "optimal", vcov = "CondHom", commonCoef=commonCoef, data=data)
        attr(res$dat, "sysInfo")$typeDesc <- "Full-Information IV (FIVE)"
        res$call <- match.call()
        res
    }

threeSLS <- function(g, h, commonCoef = FALSE, data = NULL)
    {
        if (!is(h, "formula"))
            if (length(h) != 1)
                stop("In 3SLS, h is a single equation since the same instruments are used in each equation")
        res <- sysGmm(g, h, vcov = "CondHom", commonCoef=commonCoef, data=data)
        attr(res$dat, "sysInfo")$typeDesc <- "3-stage Least Squares"
        res$call <- match.call()
        res
    }
sur <- function(g, commonCoef = FALSE, data = NULL)
    {
        if (!is.list(g))
            stop("g must be list of formulas")
        if (length(g) == 1)
            stop("For single equation GMM, use the function gmm()")
        if (!all(sapply(1:length(g), function(i) is(g[[i]], "formula"))))
            stop("g must be a list of formulas")
        tm <- lapply(1:length(g), function(i) terms(g[[i]]))
        reg <- lapply(1:length(g), function(i) attr(tm[[i]], "term.labels"))
        reg <- unique(do.call(c,reg))
        h <- paste(reg, collapse = "+")
        if (all(sapply(1:length(g), function(i) attr(tm[[i]], "intercept") == 0)))
            h <- paste(h, "-1")
        h <- as.formula(paste("~", h))
        res <- sysGmm(g, h, vcov = "CondHom", commonCoef=commonCoef, data=data)
        attr(res$dat, "sysInfo")$typeDesc <- "Seemingly Unrelated Regression (SUR)"
        res$call <- match.call()        
        res
        }
            
randEffect <- function(g, data = NULL)
    {
        res <- sur(g, commonCoef = TRUE, data = data)
        attr(res$dat, "sysInfo")$typeDesc <- "Random Effect Estimator"
        res$call <- match.call()        
        res
    }


sysGmm <- function(g, h, wmatrix = c("optimal","ident"),
                   vcov=c("MDS", "HAC", "CondHom", "TrueFixed"), 
                   kernel=c("Quadratic Spectral","Truncated", "Bartlett", "Parzen", "Tukey-Hanning"),
                   crit=10e-7,bw = bwAndrews, prewhite = FALSE, ar.method = "ols", approx="AR(1)",tol = 1e-7,
                   model=TRUE, X=FALSE, Y=FALSE, centeredVcov = TRUE, weightsMatrix = NULL,
                   data, crossEquConst = NULL, commonCoef = FALSE)
    {
        kernel <- match.arg(kernel)
        vcov <- match.arg(vcov)
        wmatrix <- match.arg(wmatrix)
        TypeGmm = "sysGmm"

        if(vcov=="TrueFixed" & is.null(weightsMatrix))
            stop("TrueFixed vcov only for fixed weighting matrix")
        if(!is.null(weightsMatrix))
            wmatrix <- "optimal"
        if(missing(data))
            data<-NULL
        all_args<-list(data = data, g = g, h = h, wmatrix = wmatrix, vcov = vcov, kernel = kernel,
                       crit = crit, bw = bw, prewhite = prewhite, ar.method = ar.method, approx = approx, 
                       weightsMatrix = weightsMatrix, centeredVcov = centeredVcov, tol = tol, 
                       model = model, X = X, Y = Y, call = match.call(), commonCoef=commonCoef, crossEquConst = crossEquConst)
        class(all_args)<-TypeGmm
        Model_info<-getModel(all_args)
        z <- momentEstim(Model_info)
        z <- FinRes(z, Model_info)
        return(z)
    }

.momentFct_Sys <- function(tet, dat)
    {
        q <- length(dat)
        f <- function(i, dat)
            {
                d <- dat[[i]]
                attr(d, "eqConst") <- attr(dat, "eqConst")
                attr(d, "ModelType") <- attr(dat, "ModelType")
                attr(d,"momentfct")  <- attr(dat,"momentfct")
                attr(d, "smooth") <- attr(dat, "smooth")
                .momentFct(tet[[i]], d)
            }
        mom <- lapply(1:q, function(i) f(i, dat))
        do.call(cbind, mom)
    }


.DmomentFct_Sys <- function(tet, dat, pt=NULL)
    {
        q <- length(dat)
        f <- function(i, dat)
            {
                d <- dat[[i]]
                attr(d, "eqConst") <- attr(dat, "eqConst")
                attr(d, "ModelType") <- attr(dat, "ModelType")
                attr(d,"momentfct")  <- attr(dat,"momentfct")
                attr(d, "smooth") <- attr(dat, "smooth")
                .DmomentFct(tet[[i]], d, pt)
            }        
        dmom <- lapply(1:q, function(i) f(i, dat))
        if (attr(dat, "sysInfo")$commonCoef)
            do.call(rbind,dmom)
        else
            .diagMatrix(dmom)
    }

.weightFct_Sys<- function(tet, dat, type=c("MDS", "HAC", "CondHom", "ident", "fct", "fixed")) 
    {
        type <- match.arg(type)
        if (type == "fixed")
            {
                w <- attr(dat, "weight")$w
                attr(w, "inv") <- FALSE
            } else if (type == "ident") {
                w <- diag(attr(dat, "q"))
                attr(w, "inv") <- FALSE
            } else {
                if (type == "HAC")
                    {
                        gt <- .momentFct_Sys(tet,dat)
                        if(attr(dat, "weight")$centeredVcov)
                            gt <- residuals(lm(gt~1))
                        n <- NROW(gt)
                        obj <- attr(dat, "weight")
                        obj$centeredVcov <- FALSE
                        w <- .myKernHAC(gt, obj)
                        attr(w, "inv") <- TRUE
                    }
                if (type == "MDS")
                    {
                        gt <- .momentFct_Sys(tet,dat)
                        n <- NROW(gt)
                        if(attr(dat, "weight")$centeredVcov)
                            gt <- scale(gt, scale=FALSE)
                        w <- crossprod(gt)/n
                        attr(w, "inv") <- TRUE
                    }
                if (type == "CondHom")
                    {
                        e <- lapply(1:length(dat), function(i) .residuals(tet[[i]], dat[[i]])$residuals)
                        e <- do.call(cbind, e)
                        Sig <- crossprod(scale(e, scale=FALSE))/nrow(e)
                        Z <- lapply(1:length(dat), function(i) dat[[i]]$x[,(2+dat[[i]]$k):ncol(dat[[i]]$x)])
                        Z <- do.call(cbind, Z)
                        w <- crossprod(Z)/nrow(e)
                        for (i in 1:length(dat))
                            for (j in 1:length(dat))
                                {
                                    s1 <- 1+(i-1)*dat[[i]]$nh
                                    e1 <- i*dat[[i]]$nh
                                    s2 <- 1+(j-1)*dat[[j]]$nh
                                    e2 <- j*dat[[j]]$nh
                                    w[s1:e1, s2:e2] <- w[s1:e1, s2:e2]*Sig[i,j]
                                }
                        attr(w, "inv") <- TRUE
                    }
            }
        return(w)
    }

.diagMatrix <- function(xlist, which=NULL)
    {
        # Create block diagonal matrix from matrices with the same number of rows.
        m <- length(xlist)
        n <- NROW(xlist[[1]])
        l <- sapply(1:m, function(i) dim(as.matrix(xlist[[i]])))
        if (!is.null(which))
            {
                if (any(l[2,1]!=l[2,]))
                    stop("diagMatrix with which given is for X with the same number of columns")
                which <- sort(which)
                if (length(which) == l[2,1])
                    {
                        dimX <- rowSums(l)
                        which <- NULL
                    } else {
                        dimX <- c(sum(l[1,]), length(which)*m + l[2,1]-length(which))
                    }
            } else {
                dimX <- rowSums(l)
            }
        X <- matrix(0, dimX[1], dimX[2])
        if (is.null(which))
            {
                for (i in 1:m)
                    {
                        s1 <- 1 + (i-1)*n
                        e1 <- n*i
                        s2 <- 1 + sum(l[2,][-(i:m)])
                        e2 <- sum(l[2,][1:i])
                        X[s1:e1, s2:e2] <- xlist[[i]]
                    }
            } else {
                q <- match(1:l[2,1],which)
                q2 <- which(is.na(q))
                wai <- 1
                for (j in q2)
                    {
                        if (wai < j)
                            {
                                xlist2 <- lapply(1:length(xlist), function(i) xlist[[i]][,wai:(j-1), drop=FALSE])
                                k <- j-wai
                                Xtmp <- .diagMatrix(xlist2)
                                X[,wai:(wai+length(xlist)*k-1)] <- Xtmp
                                wai <- wai + length(xlist)*k
                            }
                        X[,wai] <- do.call(c, lapply(1:length(xlist), function(i) xlist[[i]][,j]))
                        wai <- wai+1
                    }
                if (max(q2) < l[2,1])
                    {
                        xlist2 <- lapply(1:length(xlist), function(i) xlist[[i]][,-(1:max(q2))])
                        Xtmp <- .diagMatrix(xlist2)
                        X[,wai:ncol(X)] <- Xtmp
                    }
            }
        X
    }


.getThetaList <- function(tet, dat)
    {
        neq <- length(dat)
        if (attr(dat, "sysInfo")$commonCoef)
            {
                k <- attr(dat, "k")[[1]]
                if (length(tet) != k)
                    stop("Wrong length of tet")
                tet2 <- rep(list(c(tet)), neq)
            } else if (!is.null(attr(dat, "sysInfo")$crossEquConst)) { 
                k <- attr(dat, "k")[[1]]
                cst <- attr(dat, "sysInfo")$crossEquConst
                nk <- (k-length(cst))*neq + length(cst)
                if (nk != length(tet))
                    stop("Wrong length of tet")
                tet2 <- matrix(NA, k, neq)
                wai <- 1
                for (j in cst)
                    {
                        if (wai < j)
                            {
                                ind <- wai:(j-1)
                                tet2[ind,] <- tet[1:(neq*length(ind))]
                                tet <- tet[-(1:(neq*(j-wai)))]
                                wai <- wai + length(ind) 
                            }
                        tet2[wai,] <- tet[1]
                        tet <- tet[-1]
                        wai <- wai+1    
                    }
                if (max(cst) < k)
                        tet2[(max(cst)+1):k,] <- tet
                attr(dat, "sysInfo")$crossEquConst <- NULL
                tet2 <- .getThetaList(c(tet2), dat)
            } else {
                k2 <- do.call(c, attr(dat, "k"))
                if (sum(k2) != length(tet))
                    stop("Wrong length of tet")
                tet2 <- list()
                for (j in 1:length(k2))
                    {
                        tet2[[j]] <- tet[1:k2[j]]
                        tet <- tet[-(1:k2[j])]
                    }
            }
        tet2
    }

.chkPerfectFit <- function(obj)
    {
        r <- as.matrix(obj$residuals)
        f <- as.matrix(obj$fitted.values)
        rdf <- obj$df.residual
        rss <- colSums(r^2)
        mf <- colMeans(f)
        mv <- apply(f, 2, var)
        resvar <- rss/rdf
        resvar < (mf^2 + mv) * 1e-30 
    }
