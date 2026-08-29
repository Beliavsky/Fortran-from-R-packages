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

gmm <- function(g,x,t0=NULL,gradv=NULL, type=c("twoStep","cue","iterative"), wmatrix = c("optimal","ident"),  vcov=c("HAC","MDS","iid","TrueFixed"), 
                kernel=c("Quadratic Spectral","Truncated", "Bartlett", "Parzen", "Tukey-Hanning"),crit=10e-7,bw = bwAndrews, 
                prewhite = 1, ar.method = "ols", approx="AR(1)",tol = 1e-7, itermax=100,optfct=c("optim","optimize","nlminb", "constrOptim"),
                model=TRUE, X=FALSE, Y=FALSE, TypeGmm = "baseGmm", centeredVcov = TRUE, weightsMatrix = NULL, traceIter = FALSE, data, eqConst = NULL, 
                eqConstFullVcov = FALSE, mustar = NULL, onlyCoefficients=FALSE, ...)
    {
        type <- match.arg(type)
        kernel <- match.arg(kernel)
        vcov <- match.arg(vcov)
        wmatrix <- match.arg(wmatrix)
        optfct <- match.arg(optfct)

        if (!is.null(eqConst))
                TypeGmm <- "constGmm"
        
        if(vcov=="TrueFixed" & is.null(weightsMatrix))
            stop("TrueFixed vcov only for fixed weighting matrix")
        if(!is.null(weightsMatrix))
            wmatrix <- "optimal"

        if(missing(data))
            data<-NULL
        all_args<-list(data = data, g = g, x = x, t0 = t0, gradv = gradv, type = type, wmatrix = wmatrix, vcov = vcov, kernel = kernel,
                       crit = crit, bw = bw, prewhite = prewhite, ar.method = ar.method, approx = approx, 
                       weightsMatrix = weightsMatrix, centeredVcov = centeredVcov, tol = tol, itermax = itermax, 
                       optfct = optfct, model = model, X = X, Y = Y, call = match.call(), traceIter = traceIter, 
                       eqConst = eqConst, eqConstFullVcov = eqConstFullVcov, mustar = mustar)
        class(all_args)<-TypeGmm
        Model_info<-getModel(all_args, ...)
        z <- momentEstim(Model_info, ...)
        if (onlyCoefficients)
            return(z[c("coefficients","objective")])
        z <- FinRes(z, Model_info)
        z
    }

evalGmm <- function(g, x, t0, tetw=NULL, gradv=NULL, wmatrix = c("optimal","ident"),
                    vcov=c("HAC","iid","TrueFixed"), 
                    kernel=c("Quadratic Spectral","Truncated", "Bartlett", "Parzen",
                        "Tukey-Hanning"),crit=10e-7,bw = bwAndrews, 
                    prewhite = FALSE, ar.method = "ols", approx="AR(1)",tol = 1e-7,
                    model=TRUE, X=FALSE, Y=FALSE,  centeredVcov = TRUE,
                    weightsMatrix = NULL, data, mustar = NULL)
    {
        TypeGmm = "baseGmm"
        type <- "eval"    
        kernel <- match.arg(kernel)
        vcov <- match.arg(vcov)
        wmatrix <- match.arg(wmatrix)
        if (is.null(tetw) & is.null(weightsMatrix))
            stop("If the weighting matrix is not provided, you need to provide the vector of parameters tetw")
        if(vcov=="TrueFixed" & is.null(weightsMatrix))
            stop("TrueFixed vcov only for fixed weighting matrix")
        if(!is.null(weightsMatrix))
            wmatrix <- "optimal"
        if(missing(data))
            data<-NULL
        all_args<-list(data = data, g = g, x = x, t0 = t0, tetw = tetw, gradv = gradv,
                       type = type, wmatrix = wmatrix, vcov = vcov, kernel = kernel,
                       crit = crit, bw = bw, prewhite = prewhite, ar.method = ar.method,
                       approx = approx, weightsMatrix = weightsMatrix,
                       centeredVcov = centeredVcov, tol = tol, itermax = 100, 
                       model = model, X = X, Y = Y, call = match.call(),
                       traceIter = NULL, optfct="optim",
                       eqConst = NULL, eqConstFullVcov = FALSE, mustar = mustar)
        class(all_args)<-TypeGmm
        Model_info<-getModel(all_args)
        class(Model_info) <- "baseGmm.eval"
        z <- momentEstim(Model_info)
        z <- FinRes(z, Model_info)
        z
    }

tsls <- function(g,x,data)
    {
        if(class(g)[1] != "formula")
            stop("2SLS is for linear models expressed as formula only")
        ans <- gmm(g,x,data=data,vcov="iid", TypeGmm="tsls")
        ans$met <- "Two Stage Least Squares"
        ans$call <- match.call()
        class(ans) <- c("tsls","gmm")
        ans$vcov <- vcov(ans, type="Classical")
        return(ans)
    }

.myKernHAC <- function(gmat, obj)
    {
        gmat <- as.matrix(gmat)
        if(obj$centeredVcov) 
            gmat <- scale(gmat, scale=FALSE)
        class(gmat) <- "gmmFct"
	AllArg <- obj$WSpec$sandwich
        AllArg$x <- gmat
	if (is.function(AllArg$bw))
            {
                if (identical(AllArg$bw, bwWilhelm))
                    {
                        G <- .DmomentFct(obj$tet, obj$dat)
                        obj <- list(gt=gmat, G=G)
                        class(obj) <- "gmm"
                    } else {
                        obj <- gmat
                    }
		AllArg$bw <- AllArg$bw(obj, order.by = AllArg$order.by,
                                       kernel = AllArg$kernel, 
                                       prewhite = AllArg$prewhite,
                                       ar.method = AllArg$ar.method,
                                       approx=AllArg$approx)
            }
	weights <- weightsAndrews(x=gmat, bw=AllArg$bw, kernel=AllArg$kernel,
                                  prewhite=AllArg$prewhite, tol=AllArg$tol)
	w <- vcovHAC(x=gmat, order.by=AllArg$order.by, weights=weights,
                     prewhite=AllArg$prewhite, sandwich=FALSE,
                     ar.method=AllArg$ar.method, adjust=FALSE)
	attr(w,"Spec") <- list(weights = weights, bw = AllArg$bw,
                               kernel = AllArg$kernel)
	w
    }

getDat <- function (formula, h, data, error=TRUE) 
{
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf$na.action <- "na.pass"	
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    y <- as.matrix(model.response(mf, "numeric"))
    namey <- as.character(formula)[2]
    if (ncol(y)>1)
        namey <- paste(namey, ".", 1:ncol(y), sep="")
    xt <- as.matrix(model.matrix(mt, mf, NULL))
    n <- NROW(y)
    if (inherits(h,'formula'))
        {
            tmp <- as.character(formula)
            termsh <- terms(h)
            h <- paste(tmp[2], "~", as.character(h)[2], sep="")
            h <- as.formula(h)
            mfh <- match.call(expand.dots = FALSE)
            mh <- match(c("h", "data"), names(mfh), 0L)
            mfh <- mfh[c(1L, mh)]
            mfh$formula <- h
            mfh$h <- NULL
            mfh$drop.unused.levels <- TRUE
            mfh$na.action <- "na.pass"
            mfh[[1L]] <- as.name("model.frame")
            mfh <- eval(mfh, parent.frame())
            mth <- attr(mfh, "terms")
            h <- as.matrix(model.matrix(mth, mfh, NULL))
        }
    else
        {
            h <- as.matrix(h)
            chkInt <- sapply(1:ncol(h), function(i) all(h[,i]/mean(h[,i]) == 1))
            if (sum(chkInt) > 1)
                stop("Too many intercept in the matrix h")
            if (any(chkInt))
                {
                    h <- h[,!chkInt, drop=FALSE]
                    h <- cbind(1,h)
                    intercept_h <- TRUE
                } else {
                    if (attr(mt,"intercept")==1)
                        {
                            h <- cbind(1, h)
                            intercept_h <- TRUE
                        } else {
                            intercept_h <- FALSE
                        }
                }
            if(is.null(colnames(h)))
                {
                    if (intercept_h)
                        colnames(h) <- c("(Intercept)",paste("h",1:(ncol(h)-1),sep=""))
                    else
                        colnames(h) <- paste("h",1:ncol(h),sep="")
                } else {
                    if (intercept_h)
                        colnames(h)[1] <- "(Intercept)"
                    coln_h <- colnames(h)
                    coln_h <- gsub(" ", "", coln_h)
                    chk <- which(coln_h == "")
                    if (length(chk) >0)
                        coln_h[chk] <- paste("h", 1:length(chk), sep="")
                    if (any(duplicated(coln_h)))
                        stop("colnames of the matrix h must be unique")
                    colnames(h) <-  coln_h                    
                }
            if (!intercept_h)
                {
                    hFormula <- paste(colnames(h), collapse="+")
                    hFormula <- as.formula(paste("~", hFormula, "-1", sep=""))
                } else {
                    hFormula <- paste(colnames(h)[-1], collapse="+")
                    hFormula <- as.formula(paste("~", hFormula, sep=""))
                }
            termsh <- terms(hFormula)            
        }
    ny <- ncol(y)
    k <- ncol(xt)
    nh <- ncol(h)
    if (nh<k)
        {
            if (error)
                stop("The number of moment conditions must be at least equal to the number of coefficients to estimate")
        }
    if (is.null(colnames(y)))
        colnames(y) <- namey
    rownames(xt) <- rownames(y)
    rownames(h) <- rownames(y)
    x <- cbind(y,xt,h)
    if(any(is.na(x)))
        {
            warning("There are missing values. Associated observations have been removed")
            x <- na.omit(x)
            if (nrow(x)<=k)
                {
                    if (error)
                        stop("The number of observations must be greater than the number of coefficients")
                }
        }
    colnames(x)<-c(colnames(y),colnames(xt),colnames(h))
    return(list(x=x,nh=nh,ny=ny,k=k,mf=mf,mt=mt,cl=cl,termsh=termsh,termsx=mt))
}

.fsResOnly <- function(dat)
{
    x <- dat$x
    g <- .momentFct
    gradv <- .DmomentFct
    ny <- dat$ny
    nh <- dat$nh
    k <- dat$k
    n <- nrow(x)
    ym <- as.matrix(x[,1:ny])
    xm <- as.matrix(x[,(ny+1):(ny+k)])
    hm <- as.matrix(x[,(ny+k+1):(ny+k+nh)])
    if (!is.null(attr(dat, "eqConst")))
    {
        resTet <- attr(dat,"eqConst")$eqConst
        y2 <- xm[, resTet[,1],drop=FALSE]%*%resTet[,2]
        ym <- ym-c(y2)
        xm <- xm[,-resTet[,1],drop=FALSE]
        k <- ncol(xm)
    }
    includeExo <- which(colnames(xm)%in%colnames(hm))
    if (length(includeExo) > 0)
    {                        
        endo <- xm[, -includeExo, drop = FALSE]
        endoName <- colnames(endo)
        if (ncol(endo) != 0)
        {
            if (attr(dat$termsh, "intercept") == 1)
                restsls <- lm(endo~hm[,-1])
            else
                restsls <- lm(endo~hm-1)
            fsls <- xm
            fsls[, -includeExo] <- restsls$fitted
        } else {
            fsls <- xm
            restsls <- NULL
        }
    } else {
        if (attr(dat$termsh, "intercept") == 1)
            restsls <- lm(xm~hm[,-1])
        else
            restsls <- lm(xm~hm-1)
        endoName <- colnames(xm)
    }                
    if (!is.null(restsls))
    {
        chk <- .chkPerfectFit(restsls)
        fsRes <- suppressWarnings(summary(restsls))[!chk]
        attr(fsRes, "Endo") <- endoName[!chk]
    }
    return(fsRes)
}
    
    
.tetlin <- function(dat, w, type=NULL)
    {
        x <- dat$x
        g <- .momentFct
        gradv <- .DmomentFct
        ny <- dat$ny
        nh <- dat$nh
        k <- dat$k
        n <- nrow(x)
        ym <- as.matrix(x[,1:ny])
        xm <- as.matrix(x[,(ny+1):(ny+k)])
        hm <- as.matrix(x[,(ny+k+1):(ny+k+nh)])
        if (!is.null(attr(dat, "eqConst")))
            {
                resTet <- attr(dat,"eqConst")$eqConst
                y2 <- xm[, resTet[,1],drop=FALSE]%*%resTet[,2]
                ym <- ym-c(y2)
                xm <- xm[,-resTet[,1],drop=FALSE]
                k <- ncol(xm)
            }
        includeExo <- which(colnames(xm)%in%colnames(hm))
        inv <- attr(w, "inv")
        mustar <- attr(dat, "mustar")
        if (!is.null(type))
            {            
                if(type=="2sls")
                    {
                        if (length(includeExo) > 0)
                            {                        
                                endo <- xm[, -includeExo, drop = FALSE]
                                endoName <- colnames(endo)
                                if (ncol(endo) != 0)
                                    {
                                        if (attr(dat$termsh, "intercept") == 1)
                                            restsls <- lm(endo~hm[,-1])
                                        else
                                            restsls <- lm(endo~hm-1)
                                        fsls <- xm
                                        fsls[, -includeExo] <- restsls$fitted
                                    } else {
                                        fsls <- xm
                                        restsls <- NULL
                                    }
                            } else {
                                if (attr(dat$termsh, "intercept") == 1)
                                    restsls <- lm(xm~hm[,-1])
                                else
                                    restsls <- lm(xm~hm-1)
                                fsls <- restsls$fitted
                                endoName <- colnames(xm)
                            }                
                        par <- lm.fit(as.matrix(fsls), ym)$coefficients
                        if (ny == 1)
                            {                                    
                                e2sls <- ym-xm%*%par
                                v2sls <- lm.fit(as.matrix(hm), e2sls)$fitted
                                value <- sum(v2sls^2)/sum(e2sls^2)                
                            } else {
                                par <- c(t(par))	
                                g2sls <- g(par, dat)
                                w <- crossprod(g2sls)/n
                                gb <- matrix(colMeans(g2sls), ncol = 1)
                                value <- crossprod(gb, solve(w, gb)) 
                            }
                    }
            } else {            
                if (ny>1)
                    {
                        if (inv) 
                            {
                                whx <- solve(w, (crossprod(hm, xm) %x% diag(ny)))
                                wvecyh <- solve(w, matrix(crossprod(ym, hm), ncol = 1))	
                            } else {
                                whx <- w%*% (crossprod(hm, xm) %x% diag(ny))
                                wvecyh <- w%*%matrix(crossprod(ym, hm), ncol = 1)
                            }
                        dg <- gradv(NULL, dat)
                        xx <- crossprod(dg, whx)
                        par <- solve(xx, crossprod(dg, wvecyh))
                    } else {   
                        if (nh>k)
                            {
                                if(inv)
                                    xzwz <- crossprod(xm,hm)%*%solve(w,t(hm))	
                                else
                                    xzwz <- crossprod(xm,hm)%*%w%*%t(hm)
                                par <- solve(xzwz%*%xm,xzwz%*%ym)
                            } else {
                                par <- solve(crossprod(hm,xm),crossprod(hm,ym))
                            }
                    }
                gb <- matrix(colSums(g(par, dat))/n, ncol = 1)
                if(inv)
                    value <- crossprod(gb, solve(w, gb)) 
                else
                    value <- crossprod(gb, w%*%gb) 
            }        
        res <- list(par = par, value = value)
        if (!is.null(mustar))
            {
                if (!is.null(type))
                    {
                        w <- crossprod(hm)/NROW(hm)
                        attr(w, "inv") <- TRUE
                    }        
                res <- .mustarLin(par, xm, hm, w, dat, mustar)
            }
        if (!is.null(type))
            {    
                if (type == "2sls")
                    res$firstStageReg <- restsls
                if (!is.null(restsls))
                    {
                        chk <- .chkPerfectFit(restsls)
                        res$fsRes <- suppressWarnings(summary(restsls))[!chk]
                        attr(res$fsRes, "Endo") <- endoName[!chk]
                    }
            }
        return(res)
    }

.mustarLin <- function(par, xm, hm, w, dat, mustar)
    {
        if (ncol(xm) == ncol(hm))
            {
                par <- par-solve(crossprod(hm,xm),mustar)
            } else {
                hmxm <- crossprod(hm,xm)
                if (attr(w, "inv"))
                    T1 <- solve(w, hmxm)
                else
                    T1 <- w%*%hmxm
                par <- par-solve(crossprod(hmxm, T1), crossprod(T1, mustar))
            }
        gb <- colSums(.momentFct(par, dat))/NROW(xm)
        if(attr(w, "inv"))
            value <- crossprod(gb, solve(w, gb)) 
        else
            value <- crossprod(gb, w%*%gb)
        list(value=value, par=par)
    }

.obj1 <- function(thet, x, w)
    {
        gt <- .momentFct(thet, x)
        gbar <- as.vector(colMeans(gt))
        INV <- attr(w, "inv")
        if (INV)		
            obj <- crossprod(gbar, solve(w, gbar))
        else
            obj <- crossprod(gbar,w)%*%gbar
        return(obj)
    }

.Gf <- function(thet, x, pt = NULL)
    {
        myenv <- new.env()
        assign('x', x, envir = myenv)
        assign('thet', thet, envir = myenv)
        barg <- function(x, thet)
            {
                gt <- .momentFct(thet, x)
                if (is.null(pt))
                    gbar <- as.vector(colMeans(gt))
                else
                    gbar <- as.vector(colSums(c(pt)*gt))

                return(gbar)
            }
        G <- attr(numericDeriv(quote(barg(x, thet)), "thet", myenv), "gradient")
        return(G)
    }

.objCue <- function(thet, x, type=c("HAC", "MDS", "iid", "ident", "fct", "fixed"))
    {
        type <- match.arg(type)
        gt <- .momentFct(thet, x)
        gbar <- as.vector(colMeans(gt))
        w <- .weightFct(thet, x, type)
        if (attr(w, "inv"))
            obj <- crossprod(gbar,solve(w,gbar))
        else
            obj <- crossprod(gbar,w%*%gbar)
        return(obj)
    }	

  
.momentFct <- function(tet, dat)
    {
        if (!is.null(attr(dat, "eqConst")))
            {
                resTet <- attr(dat,"eqConst")$eqConst
                tet2 <- vector(length=length(tet)+nrow(resTet))                
                tet2[resTet[,1]] <- resTet[,2]
                tet2[-resTet[,1]] <- tet
            } else {
                tet2 <- tet
            }
        if (attr(dat, "ModelType") == "linear")
            {
                x <- dat$x
                ny <- dat$ny
                nh <- dat$nh
                k <- dat$k
                tet2 <- matrix(tet2, ncol = k)
                e <- x[,1:ny] - x[,(ny+1):(ny+k)] %*% t(tet2)
                gt <- e * x[, ny+k+1]
                if(nh > 1)
                    for (i in 2:nh)
                        gt <- cbind(gt, e*x[, (ny+k+i)])
            } else {
                gt <- attr(dat,"momentfct")(tet2, dat)
            }
        if (!is.null(attr(dat, "smooth")))
            {
                bw <- attr(dat, "smooth")$bw
                w <- attr(dat, "smooth")$w
                gt <- smoothG(gt, bw = bw, weights = w)$smoothx
            }
        gt <- as.matrix(gt)
        if (!is.null(attr(dat, "mustar")))
            {
                if (length(attr(dat, "mustar")) != ncol(gt))
                    stop("dim of mustar must match the number of moment conditions")
                gt <- sweep(gt, 2, attr(dat, "mustar"), "-")
            }
        return(gt)
    }

.DmomentFct <- function(tet, dat, pt=NULL)
    {
        if (!is.null(attr(dat, "eqConst")))
            {
                resTet <- attr(dat,"eqConst")$eqConst
                tet2 <- vector(length=length(tet)+nrow(resTet))
                tet2[resTet[,1]] <- resTet[,2]
                tet2[-resTet[,1]] <- tet
            } else {
                tet2 <- tet
            }
        if ((attr(dat,"ModelType") == "linear") & (is.null(attr(dat, "smooth"))))
            {
                x <- dat$x
                ny <- dat$ny
                nh <- dat$nh
                k <- dat$k
                dgb <- -(t(x[,(ny+k+1):(ny+k+nh)]) %*% x[,(ny+1):(ny+k)]) %x% diag(rep(1,ny))/nrow(x)
                if (!is.null(attr(dat, "eqConst")))
                    dgb <- dgb[,-attr(dat,"eqConst")$eqConst[,1], drop=FALSE]
            } else {
                if (is.null(attr(dat,"gradv")))
                    {
                        dgb <- .Gf(tet, dat, pt)
                    } else {
                        dgb <- attr(dat,"gradv")(tet2, dat)
                        if (!is.null(attr(dat, "eqConst")))
                            dgb <- dgb[,-attr(dat,"eqConst")$eqConst[,1], drop=FALSE]
                    }
            }
        return(as.matrix(dgb))
    }

.weightFct <- function(tet, dat, type=c("HAC", "MDS", "iid", "ident", "fct", "fixed")) 
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
                gt <- .momentFct(tet,dat)
                if (!is.null(attr(dat, "namesgt")))
                    colnames(gt) <- attr(dat, "namesgt")
                if(attr(dat, "weight")$centeredVcov)
                    gt <- scale(gt, scale=FALSE)
                n <- NROW(gt)
            }
        if (type == "HAC")
            {
                obj <- attr(dat, "weight")
                obj$centeredVcov <- FALSE
                obj$tet <- tet
                obj$dat <- dat
                w <- .myKernHAC(gt, obj)
                attr(w, "inv") <- TRUE
            }
        if (type == "MDS")
            {
                w <- crossprod(gt)/n
                attr(w, "inv") <- TRUE
            }
        if (type == "iid")
            {
                if (attr(dat, "ModelType") == "linear")
                    {
                        if (dat$ny == 1)
                            {
                                e <- .residuals(tet, dat)$residuals
                                sig <- mean(scale(e,scale=FALSE)^2)
                                z <- dat$x[,(1+dat$ny+dat$k):ncol(dat$x)]
                                w <- sig*crossprod(z)/length(e)
                            } else {
                                w <- crossprod(gt)/n
                            }
                    } else {
                        w <- crossprod(gt)/n
                    }
                attr(w, "inv") <- TRUE
            }        
        if (type == "fct")
            {
                w <- attr(dat, "weight")$fct(gt, attr(dat, "weight")$fctArg)
                attr(w, "inv") <- TRUE
            }
        return(w)
    }

.residuals <- function(tet, dat)
    {
        if (!is.null(attr(dat, "eqConst")))
            {
                resTet <- attr(dat,"eqConst")$eqConst
                tet2 <- vector(length=length(tet)+nrow(resTet))
                tet2[resTet[,1]] <- resTet[,2]
                tet2[-resTet[,1]] <- tet
            } else {
                tet2 <- tet
            }
        tet2 <- t(matrix(tet2, nrow = dat$ny))
        y <- as.matrix(dat$x[,1:dat$ny])
        x <- as.matrix(dat$x[,(dat$ny+1):(dat$ny+dat$k)])
        yhat <- x%*%tet2
        e <- y-yhat
        return(list(residuals=e, yhat=yhat))
    }

