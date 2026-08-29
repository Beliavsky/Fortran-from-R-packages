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

.rho <- function(x, lamb, derive = 0, type = c("EL", "ET", "CUE", "HD", "ETEL", "ETHD",
                                          "RCUE"), k = 1)
    {
	type <- match.arg(type)
        if (type == "RCUE")
            type <- "CUE"
	gml <- x%*%c(lamb)*k
	if (derive == 0)
            {
		if (type == "EL")
                    {
			if (any(gml>=1))
                            stop("Computation of Lambda fails because NAs produced by log(1-gt*l)")
			rhomat <- log(1 - gml) 
                    }
		if (type == "ET")
                    rhomat <- -exp(gml)
		if (type == "CUE")
                    rhomat <- -gml -0.5*gml^2
                if (type == "HD")
                    rhomat <- -2/(1-gml/2)
                if (type == "ETEL")
                    {
                        w <- -exp(gml)
                        w <- w/sum(w)
                        rhomat <- -log(w*NROW(gml))
                    }
                if (type == "ETHD")
                    {
                        w <- -exp(gml)
                        w <- w/sum(w)
                        rhomat <- (sqrt(w)-1/sqrt(NROW(gml)))^2
                    }
            }
	if (derive==1)
            {
		if (type == "EL")
                    rhomat <- -1/(1 - gml) 
		if (type == "ET")
                    rhomat <- -exp(gml)
		if (type == "CUE")
                    rhomat <- -1 - gml
                if (type == "HD")
                    rhomat <- -1/((1-gml/2)^2)
                if (any(type == c("ETEL","ETHD")))
                    rhomat <- NULL
            }
	if (derive==2)
            {
		if (type == "EL")
                    rhomat <- -1/(1 - gml)^2 			
		if (type == "ET")
                    rhomat <- -exp(gml)		
		if (type == "CUE")
                    rhomat <- -rep(1,nrow(x))
                if (type == "HD")
                    rhomat <- -1/((1-gml/2)^3)
                if (any(type == c("ETEL","ETHD")))
                    rhomat <- NULL                
            }
	return(c(rhomat))
    }

.getCgelLam <- function(gt, l0, type = c('EL', 'ET', 'CUE', "HD"),
                        method = c("nlminb", "optim", "constrOptim"),
                        control=list(), k = 1, alpha = 0.01)
    {
	type <- match.arg(type)
	method <- match.arg(method)
	fct <- function(l, X)
            {
		r1 <- colMeans(.rho(gt,l,derive=1,type=type,k=k)*X)
		crossprod(r1) + alpha*crossprod(l)
            }
	Dfct <- function(l, X)
            {
		r2 <- .rho(X,l,derive=2,type=type,k=k)
		r1 <- .rho(X,l,derive=1,type=type,k=k)
		H <- t(X*r2)%*%X/nrow(X)
		2*H%*%colMeans(r1*X) + 2*alpha*l
            }

	if (method == "nlminb")
            res <- nlminb(l0, fct, Dfct, X = gt, control = control) 
	if (method == "optim")
            res <- optim(l0, fct, Dfct, X = gt, method="BFGS", control = control) 
	if (method == "constrOptim")
            {
		ci <- -rep(1,nrow(gt))
		res <- constrOptim(rep(0,ncol(gt)),fct,Dfct,-gt,ci,control=control,X=gt)
            }
	if (method == "optim")
            {
		conv <- list(convergence = res$convergence,
                             counts = res$counts, message = res$message)
            } else {
		conv <- list(convergence = res$convergence, counts = res$evaluations,
                             message = res$message)
            }
        
	return(list(lambda = res$par, convergence = conv, 
                    obj = mean(.rho(gt,res$par, derive=0,type=type,k=k)) -
                        .rho(1,0, derive=0,type=type,k=k)))
    }

.Wu <- function(gt, tol_lam = 1e-8, maxiterlam = 50, K=1)
    {
        u <- as.matrix(gt)
        n=length(nrow(u))
        M=rep(0,ncol(u))
        dif=1
        tol=tol_lam
        k=0
        while(dif>tol & k<=maxiterlam)
            {
                D1=t(u)%*%((1/(1+u%*%M))*rep(1,n))
                DD=-t(u)%*%(c((1/(1+u%*%M)^2))*u)
                D2=solve(DD,D1,tol=1e-40)
                dif=max(abs(D2))
                rule=1
                while(rule>0)
                    {
                        rule=0
                        if(min(1+t(M-D2)%*%t(u))<=0) rule=rule+1
                        if(rule>0) D2=D2/2
                    }
                M=M-D2
                k=k+1
            }
        if(k>=maxiterlam)
            {
                M=rep(0,ncol(u))
                conv = list(convergence=1)
            } else {
                conv = list(convergence=0)
            }
	return(list(lambda = c(-M), convergence = conv, obj =
                        mean(.rho(gt,-M,derive=0,type="EL",k=K))))        
    }

.Wu2 <- function(gt, tol_lam = 1e-8, maxiter = 50, K = 1)
    {
        gt <- as.matrix(gt)
        res <- .Fortran(F_wu, as.double(gt), as.double(tol_lam),
                        as.integer(maxiter), as.integer(nrow(gt)),
                        as.integer(ncol(gt)), as.double(K),
                        conv=integer(1), obj=double(1),
                        lambda=double(ncol(gt)))
        list(lambda=res$lambda, convergence=list(convergence=res$conv),
             obj = res$obj)
    }

.CUE_lam <- function(gt, K=1)
    {
        q <- qr(gt)
        n <- nrow(gt)
        l0 <- -qr.coef(q, rep(1,n))
        conv <- list(convergence=0)
        list(lambda = l0, convergence = conv, obj =
                 mean(.rho(gt,l0,derive=0,type="CUE",k=K)))
    }

.CUE_lamPos <- function(gt, K=1, control=list())
    {
        getpt <- function(gt,lam)
            {
                gl <- c(gt%*%lam)
                pt <- 1 + gl
                pt/sum(pt)
            }
        maxit <- ifelse("maxit" %in% names(control),
                        control$maxit, 50)
        res <-.CUE_lam(gt, K)
        n <- nrow(gt)
        i <- 1
        pt <- getpt(gt, res$lambda)
        w <- pt<0        
        while (TRUE)
            {
                gt2 <- gt[!w,]
                n1 <- nrow(gt2)
                if (n1 == n)
                    break
                res <-  try(.CUE_lam(gt2), silent=TRUE)
                if (i > maxit)
                    return(list(lambda=rep(0,ncol(gt)), obj=0, pt=rep(1/n,n),
                                convergence=list(convergence=1)))
                if (any(class(res) == "try-error"))
                    return(list(lambda=rep(0,ncol(gt)), obj=0, pt=rep(1/n,n),
                                convergence=list(convergence=2)))
                pt[!w] <- getpt(gt2, res$lambda)
                pt[w] <- 0
                if (all(pt>=0))
                    break
                i <- i+1
                w[!w] <- pt[!w]<0
            }
        n0 <- n-n1
        res$obj <- res$obj*n1/n + n0/(2*n)
        res
    }

.CUE_lamPos2 <- function(gt, K=1, control=list())
    {
        gt <- as.matrix(gt)
        n <- nrow(gt)
        q <- ncol(gt)
        maxit <- ifelse("maxit" %in% names(control),
                        control$maxit, 50)        
        res <- try(.Fortran(F_lamcuep, as.double(gt),
                        as.integer(n), as.integer(q), as.double(K),
                        as.integer(maxit),conv=integer(1),
                        lam=double(q),pt=double(n),
                        obj=double(1)
                            ), silent=TRUE)
        if (any(class(res) == "try-error"))
            return(list(lambda=rep(0,q), obj=0, pt=rep(1/n,n),
                        convergence=list(convergence=3)))
        list(lambda=res$lam, obj=res$obj, pt=res$pt,
             convergence=list(convergence=res$conv))
    }

getLamb <- function(gt, l0, type = c('EL', 'ET', 'CUE', "ETEL", "HD", "ETHD", "RCUE"),
                    tol_lam = 1e-7, maxiterlam = 100, tol_obj = 1e-7, 
                    k = 1, method = c("nlminb", "optim", "iter", "Wu"),
                    control=list())
    {
	method <- match.arg(method)
        type <- match.arg(type)
        gt <- as.matrix(gt)
        if (method == "Wu" & type != "EL")
            stop("Wu (2005) method to compute Lambda is for EL only")
        if (method == "Wu")
            return(.Wu2(gt, tol_lam, maxiterlam, k))
        if (type == "CUE")
            return(.CUE_lam(gt, k))
        if (type == "RCUE")
            return(.CUE_lamPos2(gt, k, control))
	if (method == "iter")
            {
		if ((type == "ETEL") | (type == "ETHD"))
                    type = "ET"
		for (i in 1:maxiterlam)
                    {
			r1 <- .rho(gt,l0,derive=1,type=type,k=k)
			r2 <- .rho(gt,l0,derive=2,type=type,k=k)
			F <- -colMeans(r1*gt)
			J <- crossprod(r2*gt,gt)
			if (sum(abs(F))<tol_obj)
                            {
				conv <- list(convergence="Tolerance for the FOC reached")
				break
                            }
			P <- solve(J,F)
			if (sum(abs(P))<tol_lam)
                            {
				conv <- list(convergence="Tolerance on lambda reached")	
				break
                            }
			l0 <- l0 + P
			conv <- list(convergence="maxiterlam reached")
                    }
            } else {
		fct <- function(l,X,type,k)
                    {
			r0 <- .rho(X,l,derive=0,type=type,k=k)
			-mean(r0)
                    }
		Dfct <- function(l,X,type,k)
                    {
			r1 <- .rho(X,l,derive=1,type=type,k=k)
		        -colMeans(r1*X)
                    }
		DDfct <- function(l,X,type,k)
                    {
			r2 <- .rho(X,l,derive=2,type=type,k=k)
			-t(X*r2)%*%X/nrow(X)
                    }
		if ((type == "ETEL")|(type=="ETHD"))                
                    type = "ET"
		if (method=="optim")
                    {
                        if (type != "EL")
                            {
                                res <- optim(rep(0,ncol(gt)),fct,gr=Dfct,X=gt,type=type,
                                             k=k,method="BFGS",control=control)
                            } else {	
                                ci <- -rep(1,nrow(gt))
                                res <- constrOptim(rep(0,ncol(gt)),fct,Dfct,-gt,ci,
                                                   control=control,X=gt,type=type,k=k)
                            }
                    } else {
                        res <- nlminb(rep(0,ncol(gt)), fct, gradient = Dfct,
                                      hessian = DDfct, X = gt, type=type, k=k,
                                      control = control)
                    }
		l0 <- res$par
		if (method == "optim" | method == "constrOptim")
                    conv <- list(convergence = res$convergence, counts = res$counts,
                                 message = res$message)
		if(method == "nlminb")
                    conv <- list(convergence = res$convergence, counts =
                                     res$evaluations, message = res$message)
            }
	return(list(lambda = l0, convergence = conv, obj =
                        mean(.rho(gt,l0,derive=0,type=type,k=k))-
                            .rho(1, 0, type = type, derive = 0, k = k)))
    }

smoothG <- function (x, bw = bwAndrews, prewhite = 1, ar.method = "ols",
                     weights = weightsAndrews,
                     kernel = c("Bartlett", "Parzen", "Truncated", "Tukey-Hanning"),
                     approx = c("AR(1)", "ARMA(1,1)"),
                     tol = 1e-7) 
    {
	kernel <- match.arg(kernel)
	approx <- match.arg(approx)
        
	n <- nrow(x)
	if (is.function(weights))
            {
                class(x) <- "gmmFct"
                w <- weights(x, bw = bw, kernel = kernel,  
                             prewhite = prewhite, ar.method = ar.method, tol = tol, 
                             verbose = FALSE, approx = approx)
            } else {
                w <- weights
            }
        if (is.numeric(w))
            {
                rt <- length(w)
                if (rt >= 2)
                    {
                        w <- c(w[rt:2], w)
                        w <- w / sum(w)
                        w <- kernel(w[rt:length(w)])
                    } else {
                        w <- kernel(1)
                    }
            } else {
                if (class(w)[1] != "tskernel")                   
                    stop("Provided weights must be a numeric vector or an object of class 'tskernel'")
            }
        if (length(w$coef)>1)
            x <- kernapply(x, w)
        sx <- list("smoothx" = x, "kern_weights" = w, bw = bw)
        return(sx)		
    }
                      
gel <- function(g, x, tet0 = NULL, gradv = NULL, smooth = FALSE,
                type = c("EL", "ET", "CUE", "ETEL", "HD", "ETHD","RCUE"), 
                kernel = c("Truncated", "Bartlett"), bw = bwAndrews,
                approx = c("AR(1)", "ARMA(1,1)"), prewhite = 1, ar.method = "ols",
                tol_weights = 1e-7, tol_lam = 1e-9, tol_obj = 1e-9, 
		tol_mom = 1e-9, maxiterlam = 100, constraint = FALSE,
                optfct = c("optim", "optimize", "nlminb"), 
                optlam = c("nlminb", "optim", "iter", "Wu"), data,
                Lambdacontrol = list(),
                model = TRUE, X = FALSE, Y = FALSE, TypeGel = "baseGel", alpha = NULL,
                eqConst = NULL, eqConstFullVcov = FALSE, onlyCoefficients=FALSE, ...)
    {
	type <- match.arg(type)
	optfct <- match.arg(optfct)
	optlam <- match.arg(optlam)
	weights <- weightsAndrews
	approx <- match.arg(approx)
	kernel <- match.arg(kernel)
        if (!is.null(eqConst))
            TypeGel <- "constGel"
        
	if(missing(data))
            data<-NULL
	all_args <- list(g = g, x = x, tet0 = tet0, gradv = gradv, smooth = smooth,
                         type = type, kernel = kernel, bw = bw, approx = approx,
                         prewhite = prewhite, ar.method = ar.method,
                         tol_weights = tol_weights, tol_lam = tol_lam,
                         tol_obj = tol_obj, tol_mom = tol_mom, maxiterlam = maxiterlam,
                         constraint = constraint, optfct = optfct, weights = weights,
                         optlam = optlam, model = model, X = X, Y = Y,
                         TypeGel = TypeGel, call = match.call(),
                         Lambdacontrol = Lambdacontrol, alpha = alpha, data = data,
                         eqConst = eqConst, eqConstFullVcov = eqConstFullVcov,
                         onlyCoefficients=onlyCoefficients)
	class(all_args)<-TypeGel
	Model_info<-getModel(all_args)
	z <- momentEstim(Model_info, ...)
        if (!onlyCoefficients)
            class(z) <- "gel"
	return(z)
	}

evalGel <- function(g, x, tet0, gradv = NULL, smooth = FALSE,
                    type = c("EL", "ET", "CUE", "ETEL", "HD", "ETHD","RCUE"), 
                    kernel = c("Truncated", "Bartlett"), bw = bwAndrews,
                    approx = c("AR(1)", "ARMA(1,1)"), prewhite = 1,
                    ar.method = "ols", tol_weights = 1e-7, tol_lam = 1e-9,
                    tol_obj = 1e-9, tol_mom = 1e-9, maxiterlam = 100,
                    optlam = c("nlminb", "optim", "iter", "Wu"), data,
                    Lambdacontrol = list(),
                    model = TRUE, X = FALSE, Y = FALSE, alpha = NULL, ...)
    {
	type <- match.arg(type)
	optlam <- match.arg(optlam)
	weights <- weightsAndrews
	approx <- match.arg(approx)
	kernel <- match.arg(kernel)
        TypeGel <- "baseGel"
        
	if(missing(data))
            data<-NULL
	all_args <- list(g = g, x = x, tet0 = tet0, gradv = gradv, smooth = smooth,
                         type = type, kernel = kernel, bw = bw, approx = approx,
                         prewhite = prewhite, ar.method = ar.method, 
                         tol_weights = tol_weights, tol_lam = tol_lam,
                         tol_obj = tol_obj, tol_mom = tol_mom, maxiterlam = maxiterlam,
                         weights = weights, optlam = optlam, model = model, X = X,
                         Y = Y, call = match.call(), Lambdacontrol = Lambdacontrol,
                         alpha = alpha, data = data, optfct="optim")
	class(all_args)<-TypeGel
	Model_info<-getModel(all_args)
        class(Model_info) <- "baseGel.eval"
	z <- momentEstim(Model_info, ...)
	class(z) <- "gel"
	return(z)
      }

.thetf <- function(tet, P, output=c("obj","all"), l0Env)
    {
        output <- match.arg(output)
        gt <- P$g(tet, P$x)
        l0 <- get("l0",envir=l0Env)
        if (((P$type=="ETEL")|(P$type=="ETHD"))&(!is.null(P$CGEL)))
            {
                P$CGEL <- NULL
                warning("CGEL not implemented for ETEL or for ETHD")
            }
        if (is.null(P$CGEL))
            {
                if (P$optlam != "optim" & P$type == "EL") 
                    {
                        lamb <- try(getLamb(gt, l0, type = P$type, tol_lam = P$tol_lam,
                                            maxiterlam = P$maxiterlam, 
                                            tol_obj = P$tol_obj, k = P$k1/P$k2,
                                            control = P$Lambdacontrol, 
                                            method = P$optlam), silent = TRUE)
                        if(any(class(lamb) == "try-error"))
                            lamb <- getLamb(gt, l0, type = P$type, tol_lam = P$tol_lam,
                                            maxiterlam = P$maxiterlam, 
                                            tol_obj = P$tol_obj, k = P$k1/P$k2,
                                            control = P$Lambdacontrol, method = "optim")
                    } else {
                        lamb <- getLamb(gt, l0, type = P$type,  tol_lam = P$tol_lam,
                                        maxiterlam = P$maxiterlam, tol_obj = P$tol_obj,
                                        k = P$k1/P$k2, control = P$Lambdacontrol,
                                        method = P$optlam)
                    }
            } else {
                lamb <- try(.getCgelLam(gt, l0, type = P$type, method = "nlminb",
                                        control=P$Lambdacontrol, k = P$k1/P$k2,
                                        alpha = P$CGEL),silent=TRUE)
                if (any(class(lamb) == "try-error"))
                    lamb <- try(.getCgelLam(gt, l0, type = P$type,
                                            method = "constrOptim",
                                            control=P$Lambdacontrol, 
                                            k = P$k1/P$k2, alpha = P$CGEL),silent=TRUE)
            }
            if (P$type == "ETEL")
                obj <- mean(.rho(gt, lamb$lambda, type = P$type, derive = 0,
                                 k = P$k1/P$k2) - .rho(1, 0, type = P$type,
                                     derive = 0, k = P$k1/P$k2))
            else if (P$type == "ETHD")
                obj <- sum(.rho(gt, lamb$lambda, type = P$type, derive = 0,
                                k = P$k1/P$k2) - .rho(1, 0, type = P$type, derive = 0,
                                    k = P$k1/P$k2))
            else
                obj <- lamb$obj
        assign("l0",lamb$lambda,envir=l0Env)
        if(output == "obj")
	    return(obj)
        else
	    return(list(obj = obj, lambda = lamb, gt = gt))
    }

.getImpProb <- function(gt, lambda, type, k1, k2)
    {
        if ((type == "ETEL")|(type=="ETHD"))
            type <- "ET"
        n <- NROW(gt)
        pt <- -.rho(gt, lambda, type = type, derive = 1, k = k1/k2)/n
        # Making sure pt>0
        if (type=="CUE")
            {
                eps <- -length(pt)*min(min(pt),0)
                pt <- (pt+eps/length(pt))/(1+eps)
            }
        if (type=="RCUE")
            pt[pt<0] <- 0
        ###################
        conv_moment <- colSums(pt*gt)
        conv_pt <- sum(as.numeric(pt))
        pt <- pt/sum(pt)
        attr(pt, "conv_moment") <- conv_moment
        attr(pt, "conv_pt") <- conv_pt
        pt
    }

.vcovGel <- function(gt, G, k1, k2, bw, pt=NULL,tol=1e-16)
    {
        q <- NCOL(gt)
        n <- NROW(gt)
        if (is.null(pt))
            pt <- 1/n
        G <- G/k1
        gt <- gt*sqrt(pt*bw/k2)
        qrGt <- qr(gt)
        piv <- sort.int(qrGt$pivot, index.return=TRUE)$ix
        R <- qr.R(qrGt)[,piv]
        X <- forwardsolve(t(R), G)
        Y <- forwardsolve(t(R), diag(q))
        res <- lm.fit(X,Y)
        u <- res$residuals
        Sigma <- chol2inv(res$qr$qr)/n
        diag(Sigma)[diag(Sigma)<0] <- tol
        if (q==ncol(G))
            {
                SigmaLam <- matrix(0, q, q)
            } else {
                SigmaLam <- backsolve(R, u)/n*bw^2
                diag(SigmaLam)[diag(SigmaLam)<0] <- tol
            }
        khat <- crossprod(R)
        list(vcov_par=Sigma, vcov_lambda=SigmaLam,khat=khat)
    }
