ATEgel <- function(g, balm, w=NULL, y=NULL, treat=NULL, tet0=NULL,
                   momType=c("bal","balSample","ATT"),
                   popMom = NULL, family=c("linear","logit", "probit"),
                   type = c("EL", "ET", "CUE", "ETEL", "HD", "ETHD","RCUE"),
                   tol_lam = 1e-9, tol_obj = 1e-9, tol_mom = 1e-9, maxiterlam = 100,
                   optfct = c("optim", "nlminb"), 
                   optlam = c("nlminb", "optim", "iter", "Wu"),
                   data=NULL, Lambdacontrol = list(),
                   model = TRUE, X = FALSE, Y = FALSE, ...)
    {
        type <- match.arg(type)
	optfct <- match.arg(optfct)
	optlam <- match.arg(optlam)
        momType <- match.arg(momType)
        family <- match.arg(family)
        TypeGel <- "ateGel"

	all_args <- list(g = g, x = balm, w=w, y=y, treat=treat, tet0 = tet0, type = type,
                         tol_lam = tol_lam, tol_obj = tol_obj, tol_mom = tol_mom,
                         maxiterlam = maxiterlam, optfct = optfct, 
                         optlam = optlam, model = model, X = X, Y = Y,
                         TypeGel = TypeGel, call = match.call(),
                         Lambdacontrol = Lambdacontrol, data = data,
                         constraint=FALSE, kernel = "Truncated", bw = bwAndrews,
                         approx = "AR(1)", prewhite = 1, ar.method = "ols",
                         tol_weights = 	1e-7, alpha = NULL, eqConst = NULL,
                         eqConstFullVcov = FALSE, momType=momType, popMom=popMom,
                         family=family, onlyCoefficients=FALSE)
	class(all_args)<-TypeGel
	Model_info<-getModel(all_args)
	z <- momentEstim(Model_info, ...)        
        class(z) <- c("ategel", "gel")
        res <- try(.vcovate(z), silent=TRUE)
        if (any(class(res)=="try-error"))
            {
                warning("Could not compute the robust-to misspecification standard errors")
                z$robVcov <- FALSE
            } else {
                z$vcov_par <- res$vcov_par
                z$vcov_lambda <- res$vcov_lambda
                z$robVcov <- TRUE
            }
	return(z)
    }

.momentFctATE <- function(tet, dat)
    {
        x <- dat$x
        k <- dat$k
        if (is.null(dat$w))
            {
                Z <- x[,2:(k+1),drop=FALSE]
            } else {
                Z <- cbind(x[,2:(k+1)], dat$w)
            }
        tetz <- tet[1:ncol(Z)]
        tetb <- tail(tet, k-1)
        ZT <- c(Z%*%tetz)
        if (is.null(attr(dat, "family")))
            e <- x[,1] - ZT
        else
            e <- x[,1] - attr(dat, "family")$linkinv(ZT)
        gt1 <- e * Z
        gt2 <- sweep(x[,3:(k+1),drop=FALSE],2,tetb,"-")
        gt3 <- lapply(1:(k-1), function(i) gt2[,i]*x[,-(1:(k+1))])
        gt3 <- do.call(cbind, gt3)
        gt <- cbind(gt1,gt2,gt3)
        if (is.null(attr(dat, "popMom")))
            {
                if (attr(dat, "momType") == "balSample")
                    {
                        momB <- scale(x[,-(1:(k+1)),drop=FALSE], scale=FALSE)
                        gt <- cbind(gt, momB)
                    }
                if (attr(dat, "momType") == "ATT")
                    {                        
                        momB <- sweep(x[,-(1:3), drop=FALSE], 2,
                                      colMeans(x[x[,3]==1,-(1:3), drop=FALSE]),
                                      FUN="-")
                        gt <- cbind(gt, momB)
                    }
            } else {
                momB <- sweep(x[,-(1:(k+1)), drop=FALSE], 2, attr(dat, "popMom"), "-")
                gt <- cbind(gt, momB)
            }
        return(as.matrix(gt))
    }

.DmomentFctATE <- function(tet, dat, pt=NULL)
    {
        if (is.null(pt))
            pt <- rep(1/nrow(dat$x), nrow(dat$x))
        x <- dat$x
        k <- dat$k
        q <- dat$nh*(k-1)+2*k-1
        if (is.null(dat$w))
            {
                Z <- x[,2:(k+1),drop=FALSE]
            } else {
                Z <- cbind(x[,2:(k+1)], dat$w)
                q <- q+ncol(dat$w)
            }
        l <- ncol(Z)
        ntet <- length(tet)
        ZT <- c(Z%*%tet[1:l])
        G <- matrix(0, q, ntet)
        if (is.null(attr(dat, "family")))
            {
                tau <- rep(1, nrow(x))
            } else {                
                tau <- attr(dat, "family")$mu.eta(ZT)
            }
        G11 <- lapply(1:l, function(i) -colSums(pt*Z[,i]*tau*Z))
        G[1:l, 1:l] <- do.call(rbind, G11)
        G[(l+1):ntet, (l+1):ntet] <- -sum(pt)*diag(k-1)
        uK <- colSums(pt*x[,-(1:(k+1)),drop=FALSE])
        G[(l+k):q, (l+1):ntet] <- -kronecker(diag(k-1), uK)
        if (attr(dat, "momType") != "bal" | !is.null(attr(dat, "popMom")))
            G <- rbind(G, matrix(0, dat$nh, ntet))
        return(G)
    }


.DmomentFctATE2 <- function(tet, dat, pt=NULL)
    {
        G <- .DmomentFctATE(tet, dat, pt)
        #k <- attr(dat, "k")
        k <- dat$k
        q <- nrow(G)-dat$nh
        if (is.null(pt))
            pt <- rep(1/nrow(dat$x), nrow(dat$x))
        if (attr(dat, "momType") != "bal" & is.null(attr(dat, "popMom")))
            G <- cbind(G, rbind(matrix(0,q, dat$nh), -sum(pt)*diag(dat$nh)))
        return(G)
    }

.psiGam <- function(object)
    {
        n <- nrow(object$dat$x)
        nh <- object$dat$nh
        lam <- object$lambda
        q <- length(lam)
        k <- attr(object$dat, "k")
        theta <- object$coefficients
        gt <- object$gt
        rho1 <- .rho(x=gt, lamb=lam, derive=1, type=object$type)
        rho2 <- .rho(x=gt, lamb=lam, derive=2, type=object$type)
        if (is.null(object$dat$w))
            {
                Z <- object$dat$x[,2:(k+1)]
            } else {
                Z <- cbind(object$dat$x[,2:(k+1)], object$dat$w)
            }
        l <- ncol(Z)
        ZT <- c(Z%*%theta[1:l])
        X <- object$dat$x[,-(1:(k+1)), drop=FALSE]
        family <- attr(object$dat, "family")
        momType <- attr(object$dat, "momType")
        popMom <-  attr(object$dat, "popMom")
        if (is.null(family))
            {
                tau1 <- rep(1, n)
            } else {
                tau1 <- family$mu.eta(ZT)
                tau2 <- family$mu.eta2(ZT, family)
            }
        lG1 <- sapply(1:l, function(i) -(tau1*Z[,i]*Z)%*%lam[1:l])
        q2 <- nh*(k-1)+l+k-1
        lamM <- matrix(lam[(l+k):q2], ncol=(k-1))
        lG2 <- sapply(1:(k-1), function(i) -lam[l+i]-X%*%lamM[,i])
        lG <- cbind(lG1, lG2)
        G <- .DmomentFctATE2(theta, object$dat, rho1)
        G22 <- crossprod(rho2*gt, gt)/n
        if (momType == "bal" | !is.null(popMom))
            {
                Psi <- cbind(rho1*lG, rho1*gt)
                G11 <- crossprod(rho2*lG, lG)/n
                G12 <- t(G)/n + crossprod(rho2*lG, gt)/n
                if (!is.null(family))
                    {
                        G12tmp <- lapply(1:l, function(i)
                            colSums(-rho1*tau2*Z[,i]*c(Z%*%lam[1:l])*Z))
                        G12.2 <- matrix(0, nrow(G12), ncol(G12))
                        G12.2[1:l,1:l] <- do.call(rbind, G12tmp)
                        G12 <- G12 + G12.2/n
                    }
                Gamma <- rbind(cbind(G11, G12),
                               cbind(t(G12), G22))
                addPar <- 0
            } else {
                lG <- cbind(lG, matrix(-tail(lam, nh), n, nh, byrow=TRUE))
                G11 <- crossprod(rho2*lG, lG)/n
                G12 <- t(G)/n + crossprod(rho2*lG, gt)/n
                if (!is.null(family))
                    {
                        G12tmp <- lapply(1:l, function(i)
                            colSums(-rho1*tau2*Z[,i]*c(Z%*%lam[1:l])*Z))
                        G12.2 <- matrix(0, nrow(G12), ncol(G12))
                        G12.2[1:l,1:l] <- do.call(rbind, G12tmp)
                        G12 <- G12 + G12.2/n
                    }                
                if (momType == "balSample")
                    Xi <- rep(1,n)
                else
                    Xi <- Z[,(2:k)]
                nj <- sum(Xi)
                lam2 <- -sum(rho1)*tail(lam,nh)/nj
                theta4 <- colSums(Xi*X)/nj
                G13 <- rbind(matrix(0, l+k-1, nh), -nj/n*diag(nh))
                G23 <- matrix(0,q, nh)
                G33 <- matrix(0, nh, nh)
                Psi <- cbind(rho1*lG, rho1*gt,
                             Xi*sweep(X, 2, theta4, "-"))
                Psi[,(l+k):(l+k+nh-1)] <- Psi[,(l+k):(l+k+nh-1)]-Xi%*%t(lam2)
                Gamma <- rbind(cbind(G11, G12, G13),
                               cbind(t(G12), G22, G23),
                               cbind(t(G13), t(G23), G33))
                addPar <- nh
            }
        list(Psi=Psi, Gamma=Gamma, k=length(theta), q=q, addPar=addPar, n=n)
    }

.vcovate <- function (object) 
    {
        res <- .psiGam(object)
        k <- res$k
        q <- res$q
        addPar <- res$addPar
        qrPsi <- qr(res$Psi/sqrt(res$n))
        piv <- sort.int(qrPsi$pivot, index.return=TRUE)$ix
        R <- qr.R(qrPsi)[,piv]
        T1 <- solve(res$Gamma, t(R))
        V <- T1%*%t(T1)/res$n
        allV <- list()
        allV$vcov_par <-  V[1:k, 1:k]
        allV$vcov_lambda <- V[(k+addPar+1):(k+addPar+q), (k+addPar+1):(k+addPar+q)]
        if (addPar > 0)
            {
                allV$vcov_Allpar <-  V[1:(k+addPar), 1:(k+addPar)]
                allV$vcov_Alllambda <- V[-(1:(k+addPar)), -(1:(k+addPar))]
            }
        allV
    }


vcov.ategel <- function(object, lambda = FALSE, robToMiss=TRUE, ...)
    {
        if (robToMiss)
            {
                return(vcov.gel(object, lambda))
            } else {
                object$lambda <- rep(0, length(object$lambda))
                res <- .vcovate(object)
                object$vcov_par <- res$vcov_par
                object$vcov_lambda <- res$vcov_lambda
                return(vcov.gel(object, lambda))
            }
    }

summary.ategel <- function(object, robToMiss=TRUE, ...)
    {
        if (robToMiss)
            {
                ans <- summary.gel(object)
                ans$typeDesc = paste(ans$typeDesc,
                    "\n(S.E. are robust to misspecification)", sep="")
            } else {
                object$vcov_par <- vcov(object, robToMiss=FALSE)
                object$vcov_lambda <- vcov(object, TRUE, robToMiss=FALSE)
                ans <- summary.gel(object)
                ans$typeDesc = paste(ans$typeDesc,
                    "\n(S.E. are not robust to misspecification", sep="")
            }
        ans
    }

confint.ategel <- function (object, parm, level = 0.95, lambda = FALSE,
                            type = c("Wald", "invLR", "invLM", "invJ"), fact = 3,
                            corr = NULL, robToMiss=TRUE, ...)
    {
        type <- match.arg(type)
        if (type=="Wald")
            {
                if (!robToMiss)
                    {
                        object$vcov_par <- vcov(object, robToMiss=FALSE)
                        object$vcov_lambda <- vcov(object, TRUE, robToMiss=FALSE)
                    }
                return(confint.gel(object, parm, level, lambda, type, fact, corr, ...))
            }
        object$allArg$g <- .momentFctATE
        object$allArg$y <- NULL
        object$allArg$w <- NULL
        object$allArg$treat <- NULL
        object$allArg$popMom <- NULL
        object$allArg$momType <- NULL
        object$allArg$family <- NULL
        object$allArg$x <- object$dat        
        return(confint.gel(object, parm, level, lambda, type, fact, corr, ...))
    }


marginal <- function(object, ...)
     UseMethod("marginal")

marginal.ategel <- function(object, ...)
    {
        family <- attr(object$dat, "family")
        if (is.null(family))
            return(summary(object)$coef)
        k <- attr(object$dat, "k")
        p0 <- family$linkinv(object$coef[1])
        p1 <- family$linkinv(object$coef[1]+object$coef[2:k])
        p01 <- family$mu.eta(object$coef[1])
        p11 <- family$mu.eta(object$coef[1]+object$coef[2:k])
        A <- cbind(p11-p01, p11)
        V <- vcov(object, ...)[1:k,1:k]
        sd0 <- p01*sqrt(V[1,1])
        sdd <- sapply(1:(k-1), function(i)
            sqrt(c(t(A[i,])%*%V[c(1,i+1), c(1, i+1)]%*%A[i,])))
        coef <- cbind(c(p0,p1-p0), c(sd0, sdd))
        coef <- cbind(coef, coef[,1]/coef[,2])
        coef <- cbind(coef, 2*pnorm(-abs(coef[,3])))
        colnames(coef) <-   c("Estimate", "Std. Error", "t value", "Pr(>|t|)")            
        rownames(coef) <- c("Control",
                            paste("Treat", 1:(k-1) , " versus Control", sep=""))
        coef
    }

checkConv <- function(obj, tolConv=1e-4,verbose=TRUE, ...)
    {
        if (!any(class(obj)=="ategel"))
            stop("The function is for ategel objects produced by ATEgel()")
        momType <- obj$allArg$momType
        popMom <-  obj$allArg$popMom
        conv <- c(Lambda=obj$conv_lambda$convergence==0, Coef= obj$conv_par == 0)
        
        dat <- obj$dat$x
        nZ <- attr(obj$dat, "k")-1
        z <- dat[,3:(2+nZ),drop=FALSE]
        x <- dat[,-(1:(2+nZ)),drop=FALSE]
        pt <- getImpProb(obj, ...)
        pt1 <- lapply(1:nZ, function(i) pt[z[,i]==1]/sum(pt[z[,i]==1]))
        pt0 <- pt[rowSums(z)==0]/sum(pt[rowSums(z)==0])
        m0 <- colSums(x[rowSums(z)==0,,drop=FALSE]*pt0)
        m1 <- sapply(1:nZ, function(i) colSums(x[z[,i]==1,,drop=FALSE]*pt1[[i]]))
        mAll <- cbind(m0, m1)
        n0 <- paste(paste(colnames(z),collapse="=", sep=""),"=0",sep="")
        colnames(mAll) <- c(n0, paste(colnames(z),"=1",sep=""))
        if (!is.null(popMom))
            {
                m <- popMom
            } else {
                m <- switch(momType,
                            bal=m0,
                            balSample=colMeans(x),
                            ATT=c(m1))
            }
        chk <- all(abs(mAll-m)<tolConv)
        conv <- c(conv, Balance=all(chk))        
        if (verbose)
            {
                cat("Convergence details of the ATEgel estimation\n")
                cat("********************************************\n")
                cat(obj$typeDesc,"\n\n")
                cat("Convergence of the Lambdas: ", conv["Lambda"], "\n",sep="")
                cat("Convergence of the Coefficients: ", conv["Coef"], "\n",sep="")
                cat("Achieved moment balancing: ", conv["Balance"], "\n\n",sep="")
                cat("Moments for each group:\n")
                print.default(mAll, quote=FALSE, right=TRUE)
            }
        return(list(conv=conv, moments=mAll))
    }
