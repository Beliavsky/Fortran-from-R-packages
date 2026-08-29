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


FinRes <- function(z, object, ...)
    {
# object is computed by the getModel method #
        UseMethod("FinRes")
    }

FinRes.baseGmm.res <- function(z, object, ...)
    {
        P <- object
        x <- z$dat
        n <- ifelse(is.null(nrow(z$gt)),length(z$gt),nrow(z$gt))
        if (!is.null(attr(x,"eqConst")) & P$allArg$eqConstFullVcov)
            {
                eqConst <- attr(x,"eqConst")$eqConst
                coef <- rep(0,length(eqConst[,1])+length(z$coefficients))
                ncoef <- rep("",length(eqConst[,1])+length(z$coefficients))
                coef[-eqConst[,1]] <- z$coefficients
                ncoef[-eqConst[,1]] <- names(z$coefficients)
                coef[eqConst[,1]] <- eqConst[,2]
                ncoef[eqConst[,1]] <- rownames(eqConst)
                names(coef) <- ncoef
                z$coefficients <- coef
                if (!is.null(z$initTheta))
                    {
                        initTheta <- rep(0,length(z$coefficients))
                        initTheta[eqConst[,1]] <- eqConst[,2]
                        initTheta[-eqConst[,1]] <- z$initTheta
                        z$initTheta <- initTheta
                    }
                z$k <- z$k+nrow(eqConst)
                z$k2 <- z$k2+nrow(eqConst)
                attr(x, "eqConst") <- NULL
                z$specMod <- paste(z$specMod, "** Note: Covariance matrix computed for all coefficients based on restricted values \n   Tests non-valid**\n\n")
            }
        z$G <- z$gradv(z$coefficients, x)
        G <- z$G
        if (P$vcov == "TrueFixed")
            v <- .weightFct(z$coefficient, x, "fixed")
        else
            v <- .weightFct(z$coefficient, x, P$vcov)
        z$v <- v
        if (P$vcov == "TrueFixed") 
            {
                z$vcov=try(solve(crossprod(G, P$weightsMatrix) %*% G)/n, silent = TRUE)
                if(any(class(z$vcov) == "try-error"))
                    {
                        z$vcov <- matrix(Inf,length(z$coef),length(z$coef))
                        warning("The covariance matrix of the coefficients is singular")
                    }
            } else if ( (is.null(P$weightsMatrix)) & (P$wmatrix != "ident") ) {
                if (dim(G)[1] == dim(G)[2])
                    {
                        T1 <- try(solve(G), silent=TRUE)
                        z$vcov <- try(T1%*%v%*%t(T1)/n, silent=TRUE)
                    } else {                                                
                        z$vcov <- try(solve(crossprod(G, solve(v, G)))/n, silent = TRUE)
                    }
                if(any(class(z$vcov) == "try-error"))
                    {
                        z$vcov <- matrix(Inf,length(z$coef),length(z$coef))
                        warning("The covariance matrix of the coefficients is singular")
                    }
            } else {
                if (is.null(P$weightsMatrix))
                    w <- diag(ncol(z$gt))
                else
                    w <- P$weightsMatrix
                if (dim(G)[1] == dim(G)[2]){
                    T1 <- try(solve(G), silent=TRUE)
                } else {
                    T1 <- try(solve(t(G)%*%w%*%G,t(G)%*%w), silent = TRUE)
                }
                if(any(class(T1) == "try-error"))
                    {
                        z$vcov <- matrix(Inf,length(z$coef),length(z$coef))
                        warning("The covariance matrix of the coefficients is singular")
                    } else {
                        z$vcov <- T1%*%v%*%t(T1)/n
                    }
            }
        dimnames(z$vcov) <- list(names(z$coefficients), names(z$coefficients))
        z$call <- P$call
        
        if(is.null(P$weightsMatrix))
            {
                if(P$wmatrix == "ident")
                    {
                        z$w <- diag(ncol(z$gt))
                    } else {
                        z$w <- try(solve(v), silent = TRUE)
                        if(any(class(z$w) == "try-error"))
                            warning("The covariance matrix of the moment function is singular")
                    }
            } else {
                z$w <- P$weightsMatrix
            }
        z$weightsMatrix <- P$weightsMatrix
        z$infVcov <- P$vcov
        z$infWmatrix <- P$wmatrix
        z$allArg <- P$allArg
        if (P$wmatrix=="ident")
            z$met <- "One step GMM with W = identity"
        else
            z$met <- P$type
        z$kernel <- P$kernel
        z$coefficients <- c(z$coefficients)
        class(z) <- "gmm"
        return(z)
    }

FinRes.sysGmm.res <- function(z, object, ...)
    {
        P <- object
        x <- z$dat
        z$G <- z$gradv(z$coefficients, x)
        n <- z$n
        G <- z$G
        v <- .weightFct_Sys(z$coefficient, x, P$vcov)
        nk <- z$k
        z$v <- v
        if (P$vcov == "TrueFixed") 
            {
                z$vcov=try(solve(crossprod(G, P$weightsMatrix) %*% G)/n, silent = TRUE)
                if(any(class(z$vcov) == "try-error"))
                    {
                        z$vcov <- matrix(Inf,nk,nk)
                        warning("The covariance matrix of the coefficients is singular")
                    }
            } else if ( (is.null(P$weightsMatrix)) & (P$wmatrix != "ident") ) {
                z$vcov <- try(solve(crossprod(G, solve(v, G)))/n, silent = TRUE)
                if(any(class(z$vcov) == "try-error"))
                    {
                        z$vcov <- matrix(Inf,nk,nk)
                        warning("The covariance matrix of the coefficients is singular")
                    }
            } else {
                if (is.null(P$weightsMatrix))
                    w <- .weightFct_Sys(z$coefficients, x, "ident")
                else
                    w <- P$weightsMatrix
                if (dim(G)[1] == dim(G)[2]){
                    T1 <- try(solve(G), silent=TRUE)
                } else {
                    T1 <- try(solve(t(G)%*%w%*%G,t(G)%*%w), silent = TRUE)
                }
                if(any(class(T1) == "try-error"))
                    {
                        z$vcov <- matrix(Inf, nk, nk)
                        warning("The covariance matrix of the coefficients is singular")
                    } else {
                        z$vcov <- T1%*%v%*%t(T1)/n
                    }
            }        
        if (attr(x, "sysInfo")$commonCoef)
            dimnames(z$vcov) <- list(P$namesCoef[[1]], P$namesCoef[[1]])
        else
            dimnames(z$vcov) <- list(do.call(c, P$namesCoef), do.call(c, P$namesCoef))
        z$call <- P$call
        if(is.null(P$weightsMatrix))
            {
                if(P$wmatrix == "ident")
                    {
                        z$w <- .weightFct_Sys(z$coefficient, x, "ident")
                    } else {
                        z$w <- try(solve(v), silent = TRUE)
                        if(any(class(z$w) == "try-error"))
                            warning("The covariance matrix of the moment function is singular")
                    }
            } else {
                z$w <- P$weightsMatrix
            }
        for (i in 1:length(z$coefficients))
            names(z$coefficients[[i]]) <- P$namesCoef[[i]]
        z$weightsMatrix <- P$weightsMatrix
        z$infVcov <- P$vcov
        z$infWmatrix <- P$wmatrix
        z$allArg <- P$allArg
        z$met <- paste("System of Equations: ", P$type, sep="")
        z$kernel <- P$kernel
        z$coefficients <- c(z$coefficients)
        class(z) <- c("sysGmm", "gmm")
        return(z)
    }


