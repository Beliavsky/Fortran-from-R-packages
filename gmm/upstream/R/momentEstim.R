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


momentEstim <- function(object, ...)
    {
        UseMethod("momentEstim")
    }

momentEstim.sysGmm.twoStep.formula <- function(object, ...)
{
    dat <- object$x
    y <- lapply(1:length(dat), function(i) dat[[i]]$x[,1])
    y <- do.call(c, y)
    z <- lapply(1:length(dat), function(i)
        dat[[i]]$x[,(2+dat[[i]]$k):ncol(dat[[i]]$x), drop=FALSE])
    z <- .diagMatrix(z)
    x <- lapply(1:length(dat), function(i) dat[[i]]$x[,2:(dat[[i]]$k+1), drop=FALSE])
    if (attr(dat, "sysInfo")$commonCoef)
    {
        x <- do.call(rbind, x)
    } else if (!is.null(attr(dat, "sysInfo")$crossEquConst)) {
        k <- attr(dat, "k")[[1]]
        x <- .diagMatrix(x, (1:k)[-attr(dat, "sysInfo")$crossEquConst])
    } else {
        x <- .diagMatrix(x)
    }
    names(y) <- rownames(x) <- rownames(z) <- 1:length(y)
    df <- ncol(z) - ncol(x)
    k <- ncol(x)
    q <- ncol(z)
    n <- nrow(dat[[1]]$x)
    df.residuals <- n - k        
    data <- list(y=y, x=x, z=z)
    dat2 <- getDat(y~x-1, ~z-1, data=data)
    attr(dat2, "ModelType") <- "linear"        
    if (!is.null(object$weightsMatrix))
    {
        w <- object$weightsMatrix
        attr(w, "inv") <- FALSE
        tet0 <- NULL
        fsRes <- .fsResOnly(dat2)
    } else if (object$wmatrix == "ident") {
        w <- diag(q)
        attr(w, "inv") <- FALSE
        tet0 <- NULL
        fsRes <- .fsResOnly(dat2)
    } else {
        res0 <- .tetlin(dat2, 1, "2sls")
        tet0 <- .getThetaList(res0$par, dat)
        fsRes <- res0$fsRes
        w <- .weightFct_Sys(tet=tet0, dat=dat, type=object$vcov)
    }
    res <- .tetlin(dat2, w)
    par <- .getThetaList(res$par, dat)
    names(par) <- names(dat)
    z = list(coefficients = par, objective = res$value, dat=dat, k=k, q=q, df=df, df.residual=df.residual, n=n)
    z$gt <- object$g(z$coefficients, dat)
    z$initTheta <- tet0
    tmp <- lapply(1:length(dat), function(i) .residuals(z$coefficients[[i]], dat[[i]]))
    z$fitted.values <- lapply(1:length(dat), function(i) tmp[[i]]$yhat)                      
    z$residuals <- lapply(1:length(dat), function(i) tmp[[i]]$residuals)	
    z$terms <- lapply(1:length(dat), function(i) dat[[i]]$mt)
    if(object$model) z$model <- lapply(1:length(dat), function(i) dat[[i]]$mf)
    if(object$X) z$x <- lapply(1:length(dat), function(i)
        as.matrix(dat[[i]]$x[,(dat[[i]]$ny+1):(dat[[i]]$ny+dat[[i]]$k)]))
    if(object$Y) z$y <- lapply(1:length(dat), function(i) as.matrix(dat[[i]]$x[,1:dat[[i]]$ny]))
    z$gradv <- object$gradv
    z$g <- object$g
    z$WSpec <- object$WSpec
    z$w0 <- w
    colnames(z$gt) <- do.call(c, object$namesgt)
    z$fsRes <- fsRes
    class(z) <- "sysGmm.res"
    z$specMod <- object$specMod
    return(z)	        
}

momentEstim.baseGmm.eval <- function(object, ...)
    {
        P <- object
        x <- P$x
        n <- attr(x, "n")
        q <- attr(x, "q")
        k <- attr(x, "k")
        k2 <- k
        df <- q - k
        if (is.null(object$tetw))
            object$tetw <- object$t0
        if (is.null(attr(x, "weight")$w))
            w <- .weightFct(object$tetw, x, P$vcov)
        else
            w <- .weightFct(object$tetw, x, "fixed")
        fct <- .obj1(object$t0, x, w)    
        if (attr(x,"ModelType") == "linear")
            {
                z = list(coefficients = object$t0, objective = fct, dat = x,  k = k, k2 = k2, n = n, q = q, df = df, df.residual = (n-k))            
                b <- z$coefficients
                tmp <- .residuals(b, x)
                z$fitted.values <- tmp$yhat	
                z$residuals <- tmp$rediduals	
                z$terms <- x$mt
                if(P$model) z$model <- x$mf
                if(P$X) z$x <- as.matrix(x$x[,(x$ny+1):(x$ny+x$k)])
                if(P$Y) z$y <- as.matrix(x$x[,1:x$ny])
            } else {
                attr
                z = list(coefficients = object$t0, objective = fct, k=k, k2=k2, n=n, q=q, df=df, initTheta = object$t0, dat=x)	
            }
        z$gt <- .momentFct(z$coefficients, x)
        if (is.null(colnames(z$gt)))
            colnames(z$gt) <- paste("gt",1:ncol(z$gt),sep="")
        else
            colnames(z$gt) <- P$namesgt
        z$gradv <- P$gradv
        z$iid <- P$iid
        z$g <- P$g
        z$WSpec <- P$WSpec
        z$w0 <- w
        names(z$coefficients) <- P$namesCoef
        class(z) <- paste(P$TypeGmm,".res",sep="")
        z$specMod <- P$specMod
        
        return(z)	        
    }

momentEstim.baseGmm.twoStep <- function(object, ...)
    {
        P <- object
        x <- P$x
        n <- attr(x, "n")
        q <- attr(x, "q")
        k <- attr(x, "k")
        k2 <- k
        df <- q - k
        w = .weightFct(NULL, x, "ident")
        chkOptim <- any(P$optfct == c("optim", "constrOptim"))
        if (chkOptim)
            {
                if (P$gradvf)
                    {
                        gr2 <- function(thet, x,  w)
                            {
                                gt <- .momentFct(thet, x)
                                Gbar <- .DmomentFct(thet, x) 
                                gbar <- as.vector(colMeans(gt))
                                INV <- attr(w, "inv")
                                if (INV)		
                                    obj <- crossprod(Gbar, solve(w, gbar))
                                else
                                    obj <- crossprod(Gbar,w)%*%gbar
                                return(obj*2)
                            }
                    } else {
                        gr2 <- NULL
                    }
            }
        if (P$optfct == "optim")
            {                     
                argDots <- list(...)
                allArgOptim <- list(par = P$t0, fn = .obj1, gr = gr2, x = x, w = w)
                allArgOptim <- c(allArgOptim,argDots)
                res <- do.call(optim,allArgOptim)
            }
        if (P$optfct == "constrOptim")
            {
                if (!any(c("ui","ci") %in% names(list(...))))
                    stop("You must specify ui and ci when optfct is set to constrOptim")
                argDots <- list(...)
                ui <- argDots$ui
                ci <- argDots$ci
                argDots$ui <- NULL
                argDots$ci <- NULL
                allArgOptim <- list(theta = P$t0, f = .obj1, grad = gr2, ui = ui, ci = ci, x = x, w = w)
                allArgOptim <- c(allArgOptim,argDots)
                res <- do.call(constrOptim,allArgOptim)
            }
        if (P$optfct == "nlminb")
            {
                res <- nlminb(P$t0, .obj1, x = x, w = w, ...)
                res$value <- res$objective
            }
        if (P$optfct == "optimize")
            {
                res <- optimize(.obj1, P$t0, x = x, w = w, ...)
                res$par <- res$minimum
                res$value <- res$objective
            }
        if (q == k2 | P$wmatrix == "ident")
            {
                z = list(coefficients = res$par, objective = res$value, k=k, k2=k2, n=n, q=q, df=df)	
                if (chkOptim)
                    z$algoInfo <- list(convergence = res$convergence, counts = res$counts, message = res$message)
                else if(P$optfct == "nlminb")
                    z$algoInfo <- list(convergence = res$convergence, counts = res$evaluations, message = res$message)
            } else {
                w <- .weightFct(res$par, P$x, P$vcov)
                if (chkOptim)
                    {
                        allArgOptim$w <- w
                        res2 <- do.call(get(object$optfct),allArgOptim)
                    }
                if (P$optfct == "nlminb")
                    {
                        res2 <- nlminb(P$t0, .obj1, x = x, w = w, ...)
                        res2$value <- res2$objective
                    }                
                if (P$optfct == "optimize")
                    {
                        res2 <- optimize(.obj1, P$t0, x = x, w = w, ...)
                        res2$par <- res2$minimum
                        res2$value <- res2$objective
                    }	
                z = list(coefficients = res2$par, objective = res2$value, k=k, k2=k2, n=n, q=q, df=df, initTheta = res$par)	
                if (chkOptim)
                    {
                        z$algoInfo <- list(convergence = res2$convergence, counts = res2$counts, message = res2$message)
                        z$InitialAlgoInfo <- list(convergence = res$convergence, counts = res$counts, message = res$message)
                    } else if(P$optfct == "nlminb") {
                        z$algoInfo <- list(convergence = res2$convergence, counts = res2$evaluations, message = res2$message)
                        z$InitialAlgoInfo <- list(convergence = res$convergence, counts = res$evaluations, message = res$message)
                    }
            }
        z$dat <- P$x
        z$gt <- P$g(z$coefficients, P$x)
        z$gradv <- P$gradv
        z$iid <- P$iid
        z$g <- P$g
        z$WSpec <- P$WSpec
        names(z$coefficients) <- P$namesCoef
        if (is.null(colnames(z$gt)))
            colnames(z$gt) <- paste("gt",1:ncol(z$gt),sep="")
        class(z) <- paste(P$TypeGmm,".res",sep="")	
        z$specMod <- P$specMod 
        z$w0 <- w
        return(z)
    }


momentEstim.tsls.twoStep.formula <- function(object, ...)
    {
        P <- object
        g <- P$g
        dat <- P$x
        if (dat$ny > 1)
            stop("tsls is for one dimentional dependent variable")
        n <- attr(dat, "n")
        q <- attr(dat, "q")
        k <- attr(dat, "k")
        k2 <- k*dat$ny
        df <- q-k*dat$ny
        w = .weightFct(NULL, dat, "ident")
        if (q == k2)
            {
                res2 <- .tetlin(dat, w)
                z = list(coefficients = res2$par, objective = res2$value, dat = dat, k = k,
                    k2 = k2, n = n, q = q, df = df, df.residual = (n-k))
            } else {
                res2 <- .tetlin(dat, w, type="2sls")
                z = list(coefficients = res2$par, objective = res2$value, dat=dat, k=k, k2=k2,
                    n=n, q=q, df=df, df.residual = (n-k))	
            }
        z$gt <- g(z$coefficients, dat) 
        tmp <- .residuals(z$coefficients, dat)
        z$fitted.values <- tmp$yhat	
        z$residuals <- tmp$residuals	
        z$terms <- dat$mt
        if(P$model) z$model <- dat$mf
        if(P$X) z$x <- as.matrix(dat$x[,(dat$ny+1):(dat$ny+dat$k)])
        if(P$Y) z$y <- as.matrix(dat$x[,1:dat$ny])
        z$gradv <- P$gradv
        z$iid <- P$iid
        z$g <- P$g
        z$WSpec <- P$WSpec
        z$w0 <- w
        names(z$coefficients) <- P$namesCoef
        colnames(z$gt) <- P$namesgt        
        z$fsRes <- res2$fsRes        
        class(z) <- "baseGmm.res"
        z$specMod <- P$specMod
        return(z)	
    }


momentEstim.baseGmm.twoStep.formula <- function(object, ...)
    {
        P <- object
        g <- P$g
        dat <- P$x
        n <- attr(dat, "n")
        q <- attr(dat, "q")
        k <- attr(dat, "k")
        k2 <- k*dat$ny
        df <- q-k*dat$ny
        w = .weightFct(NULL, dat, "ident")
        if (q == k2 | P$wmatrix == "ident")
            {
                res2 <- .tetlin(dat, w)
                z = list(coefficients = res2$par, objective = res2$value, dat = dat, k = k,
                    k2 = k2, n = n, q = q, df = df, df.residual = (n-k))
            } else {
                res1 <- .tetlin(dat, w, type="2sls")
                initTheta <- res1$par
                w <- .weightFct(res1$par, dat, P$vcov)
                res2 <- .tetlin(dat, w)
                res2$firstStageReg <- res1$firstStageReg
                res2$fsRes <- res1$fsRes
                z = list(coefficients = res2$par, objective = res2$value, dat=dat, k=k, k2=k2,
                    n=n, q=q, df=df, initTheta = initTheta, df.residual = (n-k))
            }
        z$gt <- g(z$coefficients, dat) 
        b <- z$coefficients
        tmp <- .residuals(b, dat)
        z$fitted.values <- tmp$yhat	
        z$residuals <- tmp$residuals	
        z$terms <- dat$mt
        if(P$model) z$model <- dat$mf
        if(P$X) z$x <- as.matrix(dat$x[,(dat$ny+1):(dat$ny+dat$k)])
        if(P$Y) z$y <- as.matrix(dat$x[,1:dat$ny])
        z$gradv <- P$gradv
        z$iid <- P$iid
        z$g <- P$g
        z$WSpec <- P$WSpec
        z$w0 <- w
        names(z$coefficients) <- P$namesCoef
        colnames(z$gt) <- P$namesgt
        
        if (P$vcov == "iid" & P$wmatrix != "ident")
            z$fsRes <- res2$fsRes
        
        class(z) <- paste(P$TypeGmm,".res",sep="")
        z$specMod <- P$specMod
        return(z)	
    }

momentEstim.baseGmm.iterative.formula <- function(object, ...)
    {
        P <- object
        g <- P$g
        dat <- P$x
        n <- attr(dat, "n")
        q <- attr(dat, "q")
        k <- attr(dat, "k")
        k2 <- k*dat$ny
        df <- q-k*dat$ny
        w = .weightFct(NULL, dat, "ident")
  
        if (q == k2 | P$wmatrix == "ident")
            {
                res <- .tetlin(dat, w)
                z = list(coefficients = res$par, objective = res$value, dat = dat, k = k, k2 = k2,
                    n = n, q = q, df = df, df.residual = (n-k))
            } else {
                res <- .tetlin(dat, w, type="2sls")
                fsRes <- res$fsRes
                initTheta <- res$par
                ch <- 100000
                j <- 1
                while(ch > P$crit)
                    {
                        tet <- res$par
                        w <- .weightFct(tet, dat, P$vcov)
                        res <- .tetlin(dat, w)
                        ch <- crossprod(abs(tet- res$par)/tet)^.5
                        if (j>P$itermax)
                            {
                                cat("No convergence after ", P$itermax, " iterations")
                                ch <- P$crit
                            }
                        if(P$traceIter)
                            cat("Iter :",j,": value=",res$value,", Coef = ", res$par,"\n") 
                        j <- j+1	
                    }
                z = list(coefficients = res$par, objective = res$value, dat=dat, k=k, k2=k2,
                    n=n, q=q, df=df, initTheta=initTheta, df.residual = (n-k))	
            }
        z$gt <- g(z$coefficients, dat)
        tmp <- .residuals(z$coefficients, dat)
        z$fitted.values <- tmp$yhat	
        z$residuals <- tmp$residuals	
        z$terms <- dat$mt
        if(P$model) z$model <- dat$mf
        if(P$X) z$x <- as.matrix(dat$x[,(dat$ny+1):(dat$ny+dat$k)])
        if(P$Y) z$y <- as.matrix(dat$x[,1:dat$ny])  
        z$terms <- dat$mt
        if(P$model) z$model <- dat$mf
        z$gradv <- P$gradv
        z$iid <- P$iid
        z$g <- P$g
        z$WSpec <- P$WSpec
        z$w0 <- w

        names(z$coefficients) <- P$namesCoef
        colnames(z$gt) <- P$namesgt
        if (P$vcov == "iid" & P$wmatrix != "ident")
            z$fsRes <- fsRes
        class(z) <- paste(P$TypeGmm,".res",sep="")
        z$specMod <- P$specMod
        return(z)	
    }

momentEstim.baseGmm.iterative <- function(object, ...)
    {
        P <- object
        x <- P$x
        n <- attr(x, "n")
        q <- attr(x, "q")
        k <- attr(x, "k")
        k2 <- k
        df <- q - k
        w = .weightFct(NULL, x, "ident")
        chkOptim <- any(P$optfct == c("optim", "constrOptim"))
        if (chkOptim)
            {
                if (P$gradvf)
                    {
                        gr2 <- function(thet, x,  w)
                            {
                                gt <- .momentFct(thet, x)
                                Gbar <- .DmomentFct(thet, x) 
                                gbar <- as.vector(colMeans(gt))
                                INV <- attr(w, "inv")
                                if (INV)		
                                    obj <- crossprod(Gbar, solve(w, gbar))
                                else
                                    obj <- crossprod(Gbar,w)%*%gbar
                                return(obj*2)
                            }
                    } else {
                        gr2 <- NULL
                    }
            }
        if (P$optfct == "optim")
            {                     
                argDots <- list(...)
                allArgOptim <- list(par = P$t0, fn = .obj1, gr = gr2, x = x, w = w)
                allArgOptim <- c(allArgOptim,argDots)
                res <- do.call(optim,allArgOptim)
            }
        if (P$optfct == "constrOptim")
            {
                if (!any(c("ui","ci") %in% names(list(...))))
                    stop("You must specify ui and ci when optfct is set to constrOptim")
                argDots <- list(...)
                ui <- argDots$ui
                ci <- argDots$ci
                argDots$ui <- NULL
                argDots$ci <- NULL
                allArgOptim <- list(theta = P$t0, f = .obj1, grad = gr2, ui = ui, ci = ci, x = x, w = w)
                allArgOptim <- c(allArgOptim,argDots)
                res <- do.call(constrOptim,allArgOptim)
            }
        if (P$optfct == "nlminb")
            {
                res <- nlminb(P$t0, .obj1, x = x, w = w, ...)
                res$value <- res$objective
            }
        if (P$optfct == "optimize")
            {
                res <- optimize(.obj1, P$t0, x = x, w = w, ...)
                res$par <- res$minimum
                res$value <- res$objective
            }	
        if (q == k2 | P$wmatrix == "ident")
            {
                z <- list(coefficients = res$par, objective = res$value, k=k, k2=k2, n=n, q=q, df=df)
                if (chkOptim)
                    z$algoInfo <- list(convergence = res$convergence, counts = res$counts,
                                       message = res$message)
                else if(P$optfct == "nlminb")
                    z$algoInfo <- list(convergence = res$convergence, counts = res$evaluations,
                                       message = res$message)
            } else {                
                initTheta = res$par
                z <- list()
                if (chkOptim)
                    z$initialAlgoInfo <- list(convergence = res$convergence, counts = res$counts,
                                              message = res$message)
                else if(P$optfct == "nlminb")
                    z$initialAlgoInfo <- list(convergence = res$convergence, counts = res$evaluations,
                                              message = res$message)                
                ch <- 100000
                j <- 1
                while(ch > P$crit)
                    {
                        tet <- res$par
                        w <- .weightFct(tet, P$x, P$vcov)          
                        if (chkOptim)
                            {
                                allArgOptim$w <- w
                                allArgOptim[[1]] <- tet
                                res <- do.call(get(P$optfct),allArgOptim)
                            }
                        if (P$optfct == "nlminb")
                            {
                                res <- nlminb(tet, .obj1, x = P$x, w = w, ...)
                                res$value <- res$objective
                            }
                        if (P$optfct == "optimize")
                            {
                                res <- optimize(.obj1, P$t0, x = P$x, w = w, ...)
                                res$par <- res$minimum
                                res$value <- res$objective
                            }	
                        ch <- crossprod(tet-res$par)^.5/(1+crossprod(tet)^.5)
                        if (j>P$itermax)
                            {
                                cat("No convergence after ", P$itermax, " iterations")
                                ch <- P$crit
                            }
                        if(P$traceIter)
                            cat("Iter :",j,": value=",res$value,", Coef = ", res$par,"\n") 
                        j <- j+1	
                    }
                z2 = list(coefficients = res$par, objective = res$value,k=k, k2=k2, n=n, q=q,
                    df=df, initTheta=initTheta)
                z <- c(z, z2)
                if (chkOptim)
                    z$algoInfo <- list(convergence = res$convergence, counts = res$counts,
                                       message = res$message)
                else if(P$optfct == "nlminb")
                    z$algoInfo <- list(convergence = res$convergence, counts = res$evaluations,
                                       message = res$message)
                
            }
        z$dat <- P$x
        z$gt <- P$g(z$coefficients, P$x)
        z$gradv <- P$gradv
        z$iid <- P$iid
        z$g <- P$g
        z$WSpec <- P$WSpec
        z$w0 <- w

        names(z$coefficients) <- P$namesCoef
        if (is.null(colnames(z$gt)))
            colnames(z$gt) <- paste("gt",1:ncol(z$gt),sep="")
        z$specMod <- P$specMod
        class(z) <- paste(P$TypeGmm,".res",sep="")	
        return(z)
    }

momentEstim.baseGmm.cue.formula <- function(object, ...)
    {
        P <- object
        g <- P$g  
        dat <- P$x
        n <- attr(dat, "n")
        q <- attr(dat, "q")
        k <- attr(dat, "k")  
        k2 <- k*dat$ny
        df <- q-k*dat$ny
        w <- .weightFct(NULL, dat, "ident")                  
        
        if (q == k2 | P$wmatrix == "ident")
            {
                res <- .tetlin(dat, w)
                z = list(coefficients = res$par, objective = res$value, dat = dat, k = k, k2 = k2, n = n, q = q, df = df, df.residual = (n-k))
                P$weightMessage <- "No CUE needed because the model is just identified"
            } else {
                if (is.null(P$t0))
                    {
                        P$t0 <- .tetlin(dat, w, type="2sls")$par
                        initTheta <- P$t0
                    } else {
                        initTheta <- P$t0
                    }
                if (P$vcov == "HAC")
                    {
                        w <- .weightFct(P$t0, dat, P$vcov)
                        attr(dat, "weight")$WSpec$sandwich$bw <- attr(w,"Spec")$bw
                        P$weightMessage <- "Weights for kernel estimate of the covariance are fixed and based on the first step estimate of Theta"
                    }
                if (P$optfct == "optim")
                    res2 <- optim(P$t0,.objCue, x = dat, type = P$vcov, ...)
                if (P$optfct == "constrOptim")
                    {
                        if (!any(c("ui","ci") %in% names(list(...))))
                            stop("You must specify ui and ci when optfct is set to constrOptim")
                        argDots <- list(...)
                        ui <- argDots$ui
                        ci <- argDots$ci
                        argDots$ui <- NULL
                        argDots$ci <- NULL
                        allArgOptim <- list(theta = P$t0, f = .objCue, grad = NULL, ui = ui, ci = ci, x = dat, type = P$vcov)
                        allArgOptim <- c(allArgOptim,argDots)
                        res2 <- do.call(constrOptim,allArgOptim)
                    }
                if (P$optfct == "nlminb")
                    {
                        res2 <- nlminb(P$t0,.objCue, x = dat, type = P$vcov, ...)
                        res2$value <- res2$objective
                    }
                if (P$optfct == "optimize")
                    {
                        res2 <- optimize(.objCue,P$t0, x = dat, type = P$vcov, ...)
                        res2$par <- res2$minimum
                        res2$value <- res2$objective
                    }
                z = list(coefficients = res2$par, objective = res2$value, dat = dat, k = k, k2 = k2,
                    n = n, q = q, df = df, initTheta=initTheta, df.residual = (n-k))
                if (any(P$optfct == c("optim", "constrOptim")))
                    z$algoInfo <- list(convergence = res2$convergence, counts = res2$counts, message = res2$message)
                else if(P$optfct == "nlminb")
                    z$algoInfo <- list(convergence = res2$convergence, counts = res2$evaluations, message = res2$message)
            }

        z$gt <- g(z$coefficients, dat)
        tmp <- .residuals(z$coefficients, dat)
        z$fitted.values <- tmp$yhat	
        z$residuals <- tmp$residuals	
        z$terms <- dat$mt
        if(P$model) z$model <- dat$mf
        if(P$X) z$x <- as.matrix(dat$x[,(dat$ny+1):(dat$ny+dat$k)])
        if(P$Y) z$y <- as.matrix(dat$x[,1:dat$ny])  
        z$dat <- dat 
        z$terms <- dat$mt
        if(P$model) z$model <- dat$mf
        z$gradv <- P$gradv
        z$iid <- P$iid
        z$g <- P$g
        z$specMod <- P$specMod
        z$cue <- list(weights=P$fixedKernW,message=P$weightMessage)
        z$WSpec <- P$WSpec
        z$w0 <- .weightFct(z$coefficients, dat, P$vcov)
        names(z$coefficients) <- P$namesCoef
        colnames(z$gt) <- P$namesgt
        
        class(z) <- paste(P$TypeGmm,".res",sep="")
        return(z)	
    }

momentEstim.baseGmm.cue <- function(object, ...)
    {
        P <- object
        x <- P$x
        n <- attr(x, "n")
        q <- attr(x, "q")
        k <- attr(x, "k")
        k2 <- k
        df <- q - k

        if (q == k2 | P$wmatrix == "ident")
            {
                w = .weightFct(NULL, x, "ident")
                res <- gmm(P$allArg$g,P$allArg$x,P$t0,wmatrix="ident",optfct=P$optfct, ...)
                z <- list(coefficients = res$coef, objective = res$objective,
                          algoInfo = res$algoInfo, k=k, k2=k2, n=n, q=q, df=df,
                          initTheta=P$t0)
                P$weightMessage <- "No CUE needed because the model if just identified or you set wmatrix=identity"
            } else {
                w <- .weightFct(P$t0, x, P$vcov)
                initTheta <- P$t0
                if (P$vcov == "HAC")
                    {
                        res <- try(gmm(P$allArg$g,P$allArg$x,P$t0,wmatrix="ident",
                                       optfct=P$optfct, ...))                      
                        if(any(class(res)=="try-error"))
                            stop("Cannot get a first step estimate to compute the weights for the Kernel estimate of the covariance matrix; try different starting values")
                        w <- .weightFct(res$coefficients, x, P$vcov)
                        attr(x, "weight")$WSpec$sandwich$bw <- attr(w,"Spec")$bw
                        P$weightMessage <- "Weights for kernel estimate of the covariance are fixed and based on the first step estimate of Theta"
                    } else {
                        res <- list()
                    }
                if (P$optfct == "optim")
                    {
                        res2 <- optim(P$t0, .objCue, x = x, type = P$vcov, ...)
                    }
                if (P$optfct == "constrOptim")
                    {
                        if (!any(c("ui","ci") %in% names(list(...))))
                            stop("You must specify ui and ci when optfct is set to constrOptim")
                        argDots <- list(...)
                        ui <- argDots$ui
                        ci <- argDots$ci
                        argDots$ui <- NULL
                        argDots$ci <- NULL
                        allArgOptim <- list(theta = P$t0, f = .objCue, grad = NULL, ui = ui, ci = ci, x = x, type = P$vcov)
                        allArgOptim <- c(allArgOptim,argDots)
                        res2 <- do.call(constrOptim,allArgOptim)
                    }
                if (P$optfct == "nlminb")
                    {
                        res2 <- nlminb(P$t0, .objCue, x = x, type = P$vcov, ...)
                        res2$value <- res2$objective
                    }
                if (P$optfct == "optimize")
                    {
                        res2 <- optimize(.objCue,P$t0, x = x, type = P$vcov, ...)
                        res2$par <- res2$minimum
                        res2$value <- res2$objective
                    }
                z = list(coefficients=res2$par,objective=res2$value, k=k, k2=k2,
                    n=n, q=q, df=df, initTheta=initTheta)
                if (any(P$optfct == c("optim", "constrOptim")))
                    {
                        z$algoInfo <- list(convergence = res2$convergence, counts =
                                               res2$counts, message = res2$message)
                        z$InitialAlgoInfo <- list(convergence = res$algoInfo$convergence,
                                                  counts = res$algoInfo$counts,
                                                  message = res$algoInfo$message)
                    } else if (P$optfct == "nlminb") {
                        z$algoInfo <- list(convergence = res2$convergence, counts =
                                               res2$evaluations, message = res2$message)
                        z$InitialAlgoInfo <- list(convergence = res$algoInfo$convergence,
                                                  counts = res$algoInfo$evaluations,
                                                  message = res$algoInfo$message)
                    }
            }
        z$dat <- x
        z$gradv <- P$gradv
        z$gt <- P$g(z$coefficients, x)
        z$w0 <- .weightFct(z$coefficients, x, P$vcov)        
        z$iid <- P$iid
        z$g <- P$g
        z$cue <- list(weights=P$fixedKernW,message=P$weightMessage)
        names(z$coefficients) <- P$namesCoef
        if (is.null(colnames(z$gt)))
            colnames(z$gt) <- paste("gt",1:ncol(z$gt),sep="")
        z$WSpec <- P$WSpec        
        z$specMod <- P$specMod
        class(z) <- paste(P$TypeGmm, ".res", sep = "")	
        return(z)
    }

momentEstim.baseGel.modFormula <- function(object, ...)
    {
        P <- object
        g <- P$g
        q <- attr(P$x, "q")
        n <- attr(P$x, "n")
        k <- attr(P$x, "k")
        df <- q-k*P$x$ny
        l0Env <- new.env()
        assign("l0",rep(0,q),envir=l0Env)
        if (!P$constraint)
            {
                if (P$optfct == "optim")
                    res <- optim(P$tet0, .thetf, P = P, l0Env = l0Env, ...)
                if (P$optfct == "nlminb")
                    res <- nlminb(P$tet0, .thetf, P = P, l0Env = l0Env, ...)	
                if (P$optfct == "optimize")
                    { 
                        res <- optimize(.thetf, P$tet0, P = P, l0Env = l0Env, ...)
                        res$par <- res$minimum
                        res$convergence <- "There is no convergence code for optimize"
                    }
            }

        if(P$constraint)
            res <- constrOptim(P$tet0, .thetf, grad = NULL, P = P, l0Env = l0Env, ...)
        All <- .thetf(res$par, P, "all",l0Env = l0Env)
        gt <- All$gt
        rlamb <- All$lambda

        z <- list(coefficients = res$par, lambda = rlamb$lambda, conv_lambda = rlamb$conv, conv_par = res$convergence, dat=P$x)

        z$type <- P$type
        z$gt <- gt
        pt <- .getImpProb(z$gt, z$lambda, P$type, P$k1, P$k2)
        z$pt <- c(pt) 
        z$conv_moment <- attr(pt, "conv_moment")
        z$conv_pt <- attr(pt, "conv_pt")
        z$objective <- All$obj
        z$call <- P$call
        z$k1 <- P$k1
        z$k2 <- P$k2
        z$CGEL <- P$CGEL
        z$typeDesc <- P$typeDesc
        z$specMod <- P$specMod
        z$df <- df
        names(z$coefficients) <- object$namesCoef
        if (P$onlyCoefficients)
            return(z[c("coefficients","lambda","conv_lambda","conv_par","objective")])

        
        if (!is.null(object$namesgt))
            {
                colnames(z$gt) <- object$namesgt
            } else {
                colnames(z$gt) <- paste("g",1:ncol(z$gt), sep="")
            }
        names(z$lambda) <- paste("Lam(", colnames(z$gt), ")", sep="")

        if (!is.null(attr(P$x,"eqConst")) & P$allArg$eqConstFullVcov)
            {
                eqConst <- attr(P$x,"eqConst")$eqConst
                coef <- rep(0,length(eqConst[,1])+length(z$coefficients))
                ncoef <- rep("",length(eqConst[,1])+length(z$coefficients))
                coef[-eqConst[,1]] <- z$coefficients
                ncoef[-eqConst[,1]] <- names(z$coefficients)
                coef[eqConst[,1]] <- eqConst[,2]
                ncoef[eqConst[,1]] <- rownames(eqConst)
                names(coef) <- ncoef
                z$coefficients <- coef
                attr(P$x, "k") <- attr(P$x, "k") + nrow(eqConst)
                z$df <- z$df - nrow(eqConst)
                attr(P$x,"eqConst") <- NULL
                z$specMod <- paste(z$specMod, "** Note: Covariance matrix computed for all coefficients based on restricted values \n   Tests non-valid**\n\n")
            }
  
        if(P$gradvf)
            G <- P$gradv(z$coefficients, P$x)
        else
            G <- P$gradv(z$coefficients, P$x, z$pt)
        allVcov <- try(.vcovGel(gt, G, P$k1, P$k2, P$bwVal, z$pt),
                       silent=TRUE)
        if (any(class(allVcov) == "try-error"))
            {
                z$vcov_par <- matrix(NA, length(z$coefficients), length(z$coefficients))
                z$vcov_lambda <- matrix(NA, length(z$lambda), length(z$lambda))
                z$khat <- NULL
                warning("Cannot compute the covariance matrices")
            } else {
                z <- c(z, allVcov)
            }        
        z$weights <- P$w
        z$bwVal <- P$bwVal
        names(z$bwVal) <- "Bandwidth"
        dimnames(z$vcov_par) <- list(names(z$coefficients), names(z$coefficients))
        dimnames(z$vcov_lambda) <- list(names(z$lambda), names(z$lambda))
        tmp <- .residuals(z$coefficients, P$x)
        z$fitted.values <- tmp$yhat	
        z$residuals <- tmp$residuals
        z$dat <- P$x
        z$terms <- P$x$mt
        if(P$model) z$model <- P$x$mf
        if(P$X) z$x <- as.matrix(P$x$x[,(P$x$ny+1):(P$x$ny+P$x$k)])
        if(P$Y) z$y <- as.matrix(P$x$x[,1:P$x$ny])  
        class(z) <- paste(P$TypeGel, ".res", sep = "")
        z$allArg <- P$allArg
        return(z)
    }

momentEstim.baseGel.mod <- function(object, ...)
    {
        P <- object
        x <- P$x
        q <- attr(x, "q")
        n <- attr(x, "n")
        df <- q - attr(x, "k")
        l0Env <- new.env()
        assign("l0",rep(0,q),envir=l0Env)
        if (!P$constraint)
            {
                if (P$optfct == "optim")
                    res <- optim(P$tet0, .thetf, P = P, l0Env = l0Env, ...)
                if (P$optfct == "nlminb")
                    res <- nlminb(P$tet0, .thetf, P = P, l0Env = l0Env, ...)                
                if (P$optfct == "optimize")
                    { 
                        res <- optimize(.thetf, P$tet0, P = P, l0Env = l0Env, ...)
                        res$par <- res$minimum
                        res$convergence <- "There is no convergence code for optimize"
                    }
            }
        if(P$constraint)
            res <- constrOptim(P$tet0, .thetf, grad = NULL, P = P,l0Env = l0Env, ...)
        All <- .thetf(res$par, P, "all",l0Env = l0Env)
        gt <- All$gt
        rlamb <- All$lambda

        z <- list(coefficients = res$par, lambda = rlamb$lambda, conv_lambda = rlamb$conv, conv_par = res$convergence, dat=x)

        z$type <- P$type
        z$gt <- gt
        pt <- .getImpProb(z$gt, z$lambda, P$type, P$k1, P$k2)
        z$pt <- c(pt) 
        z$conv_moment <- attr(pt, "conv_moment")
        z$conv_pt <- attr(pt, "conv_pt")
        z$objective <- All$obj
        names(z$coefficients) <- P$namesCoef
        z$specMod <- P$specMod
        z$df <- df
        if (!is.null(object$namesgt))
            {
                colnames(z$gt) <- object$namesgt
            } else {
                colnames(z$gt) <- paste("g",1:ncol(z$gt), sep="")
            }
        names(z$lambda) <- paste("Lam(", colnames(z$gt), ")", sep="")
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
                attr(x, "k") <- attr(x, "k") +  nrow(eqConst)
                z$df <- z$df -  nrow(eqConst)
                attr(x,"eqConst") <- NULL
                z$specMod <- paste(z$specMod, "** Note: Covariance matrix computed for all coefficients based on restricted values \n   Tests non-valid**\n\n")
            }
        if (P$onlyCoefficients)
            return(z[c("coefficients", "lambda", "conv_lambda", "conv_par", "objective")])
        
        if(P$gradvf)
            G <- P$gradv(z$coefficients, x)
        else
            G <- P$gradv(z$coefficients, x, z$pt)
        z$G <- G
        allVcov <- try(.vcovGel(gt, G, P$k1, P$k2, P$bwVal, z$pt),
                         silent=TRUE)
        if (any(class(allVcov) == "try-error"))
            {
                z$vcov_par <- matrix(NA, length(z$coefficients), length(z$coefficients))
                z$vcov_lambda <- matrix(NA, length(z$lambda), length(z$lambda))
                z$khat <- NULL
                warning("Cannot compute the covariance matrices")
            } else {
                z <- c(z, allVcov)
            } 
        z$weights <- P$w
        z$bwVal <- P$bwVal
        names(z$bwVal) <- "Bandwidth"
        dimnames(z$vcov_par) <- list(names(z$coefficients), names(z$coefficients))
        dimnames(z$vcov_lambda) <- list(names(z$lambda), names(z$lambda))
        if(P$X) z$x <- x
        z$call <- P$call
        z$k1 <- P$k1
        z$k2 <- P$k2
        z$CGEL <- P$CGEL
        z$typeDesc <- P$typeDesc
        z$allArg <- P$allArg        
        class(z) <- paste(P$TypeGel, ".res", sep = "")
        return(z)
    }

momentEstim.fixedW.formula <- function(object, ...)
    {
        P <- object
        g <- P$g
        dat <- P$x
        n <- attr(dat, "n")
        q <- attr(dat, "q")
        k <- attr(dat, "k")
        k2 <- k*dat$ny
        df <- q-k*dat$ny
        w = .weightFct(NULL, dat, "fixed")
        if(!all(dim(w) == c(q,q)))
            stop("The matrix of weights must be qxq")
        eigenW <- svd(w)$d
        if(!is.double(eigenW))
            warning("The matrix of weights is not strictly positive definite")
        if(is.double(eigenW))
            {
                if(any(eigenW<=0))
                    warning("The matrix of weights is not strictly positive definite")
            }  
        res2 <- .tetlin(dat, w)
        
        z = list(coefficients = res2$par, objective = res2$value, dat=dat, k=k, k2=k2, n=n, q=q, df=df, df.residual = (n-k))	

        z$gt <- g(z$coefficients, dat)
        tmp <- .residuals(z$coefficients, dat)
        z$fitted.values <- tmp$yhat	
        z$residuals <- tmp$rediduals	
        z$terms <- dat$mt
        if(P$model) z$model <- dat$mf
        if(P$X) z$x <- as.matrix(dat$x[,(dat$ny+1):(dat$ny+dat$k)])
        if(P$Y) z$y <- as.matrix(dat$x[,1:dat$ny])  
        z$dat <- dat 
        z$terms <- dat$mt
        if(P$model) z$model <- dat$mf
        z$gradv <- P$gradv
        z$iid <- P$iid
        z$g <- P$g
        z$WSpec <- P$WSpec
        
        names(z$coefficients) <- P$namesCoef
        colnames(z$gt) <- P$namesgt
        
        z$specMod <- P$specMod
        class(z) <- paste(P$TypeGmm,".res",sep="")
        return(z)	
    }

momentEstim.fixedW <- function(object, ...)
    {
        P <- object
        x <- P$x
        n <- attr(x, "n")
        q <- attr(x, "q")
        k <- attr(x, "k")
        k2 <- k
        df <- q - k
        w = .weightFct(NULL, x, "fixed")
        if(!all(dim(w) == c(q,q)))
            stop("The matrix of weights must be qxq")
        eigenW <- svd(w)$d
        if(!is.double(eigenW))
            warning("The matrix of weights is not strictly positive definite")
        if(is.double(eigenW))
            {
                if(any(eigenW<=0))
                    warning("The matrix of weights is not strictly positive definite")
            }
        chkOptim <- any(P$optfct == c("optim", "constrOptim"))
        if (chkOptim)
            {
                if (P$gradvf)
                    {
                        gr2 <- function(thet, x,  w)
                            {
                                gt <- .momentFct(thet, x)
                                Gbar <- .DmomentFct(thet, x) 
                                gbar <- as.vector(colMeans(gt))
                                INV <- attr(w, "inv")
                                if (INV)		
                                    obj <- crossprod(Gbar, solve(w, gbar))
                                else
                                    obj <- crossprod(Gbar,w)%*%gbar
                                return(obj*2)
                            }
                    } else {
                        gr2 <- NULL
                    }
            }
        if (P$optfct == "optim")
            {
                argDots <- list(...)
                allArgOptim <- list(par = P$t0, fn = .obj1, gr = gr2, x = x, w = w)
                argDots$gr <- NULL
                allArgOptim <- c(allArgOptim,argDots)
                res2 <- do.call(optim,allArgOptim)
            }    
        if (P$optfct == "constrOptim")
            {
                if (!any(c("ui","ci") %in% names(list(...))))
                    stop("You must specify ui and ci when optfct is set to constrOptim")
                argDots <- list(...)
                ui <- argDots$ui
                ci <- argDots$ci
                argDots$ui <- NULL
                argDots$ci <- NULL
                allArgOptim <- list(theta = P$t0, f = .obj1, grad = gr2, ui = ui, ci = ci, x = x, w = w)
                allArgOptim <- c(allArgOptim,argDots)
                res2 <- do.call(constrOptim,allArgOptim)
            }
        if (P$optfct == "nlminb")
            {
                res2 <- nlminb(P$t0, .obj1, x = P$x, w = w, ...)
                res2$value <- res2$objective
            }
        if (P$optfct == "optimize")
            {
                res2 <- optimize(.obj1, P$t0, x = P$x, w = w, ...)
                res2$par <- res2$minimum
                res2$value <- res2$objective
            }	
        z = list(coefficients = res2$par, objective = res2$value, k=k, k2=k2, n=n, q=q, df=df)	
        if (chkOptim)
            z$algoInfo <- list(convergence = res2$convergence, counts = res2$counts, message = res2$message)
        else if(P$optfct == "nlminb")
            z$algoInfo <- list(convergence = res2$convergence, counts = res2$evaluations, message = res2$message)

        z$dat <- P$x
        z$gt <- P$g(z$coefficients, P$x)
        z$gradv <- P$gradv
        z$iid <- P$iid
        z$g <- P$g
        z$WSpec <- P$WSpec

        names(z$coefficients) <- P$namesCoef
        if (is.null(colnames(z$gt)))
            colnames(z$gt) <- paste("gt",1:ncol(z$gt),sep="") 
        z$specMod <- P$specMod
        class(z) <- paste(P$TypeGmm,".res",sep="")	
        return(z)
    }

momentEstim.baseGel.eval <- function(object, ...)
    {
        P <- object
        q <- attr(P$x, "q")
        n <- attr(P$x, "n")
        l0Env <- new.env()
        assign("l0",rep(0,q),envir=l0Env)
        All <- .thetf(P$tet0, P, "all",l0Env = l0Env)
        gt <- All$gt
        rlamb <- All$lambda
        z <- list(coefficients = P$tet0, lambda = rlamb$lambda, conv_lambda = rlamb$conv, conv_par = NULL, dat=P$x)

        z$type <- P$type
        z$gt <- gt
        pt <- .getImpProb(z$gt, z$lambda, P$type, P$k1, P$k2)
        z$pt <- c(pt) 
        z$conv_moment <- attr(pt, "conv_moment")
        z$conv_pt <- attr(pt, "conv_pt")
        z$objective <- All$obj
        z$call <- P$call
        z$k1 <- P$k1
        z$k2 <- P$k2
        z$CGEL <- P$CGEL
        z$typeDesc <- paste(P$typeDesc, " (Eval only, tests non-valid) ", sep="")
        z$specMod <- P$specMod
        names(z$coefficients) <- P$namesCoef
        z$df <- length(z$lambda) - length(z$coefficients) 
        if (!is.null(object$namesgt))
            {
                colnames(z$gt) <- object$namesgt
            } else {
                colnames(z$gt) <- paste("g",1:ncol(z$gt), sep="")
            }
        names(z$lambda) <- paste("Lam(", colnames(z$gt), ")", sep="")
        if(P$gradvf)
            G <- P$gradv(z$coefficients, P$x)
        else
            G <- P$gradv(z$coefficients, P$x, z$pt)
        z$G <- G
        allVcov <- .vcovGel(gt, G, P$k1, P$k2, P$bwVal, z$pt)
        z <- c(z, allVcov)       
        z$weights <- P$w
        z$bwVal <- P$bwVal
        names(z$bwVal) <- "Bandwidth"
        dimnames(z$vcov_par) <- list(names(z$coefficients), names(z$coefficients))
        dimnames(z$vcov_lambda) <- list(names(z$lambda), names(z$lambda))
        if (attr(P$x,"ModelType") == "linear")
            {
                tmp <- .residuals(z$coefficients, P$x)
                z$fitted.values <- tmp$yhat	
                z$residuals <- tmp$residuals
                z$terms <- P$x$mt
                if(P$model) z$model <- P$x$mf
                if(P$X) z$x <- as.matrix(P$x$x[,(P$x$ny+1):(P$x$ny+P$x$k)])
                if(P$Y) z$y <- as.matrix(P$x$x[,1:P$x$ny])
            } 
        return(z)
    }


