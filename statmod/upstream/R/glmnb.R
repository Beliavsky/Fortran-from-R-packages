glmnb.fit <- function(X,y,dispersion,weights=NULL,offset=0,coef.start=NULL,start.method="mean",tol=1e-6,maxit=50,trace=FALSE)
#  Fit negative binomial generalized linear model with log link
#  by Fisher scoring with Levenberg-style damped
#  Gordon Smyth and Yunshun Chen
#  Created 2 November 2010.  Last modified 29 May 2019.
{
#	Check input values for y
	y <- as.vector(y)
	if(any(y < 0)) stop("y must be non-negative")
	if(!all(is.finite(y))) stop("All y values must be finite and non-missing")
	ymax <- max(y)
	n <- length(y)

#	Handle zero length y as special case
	if(n == 0) stop("y has length zero")

#	Check input values for X
	X <- as.matrix(X)
	if(n != nrow(X)) stop("length(y) not equal to nrow(X)")
	if(!all(is.finite(X))) stop("All X values must be finite and non-missing")
	p <- ncol(X)
	if(p > n) stop("More columns than rows in X")
	if(is.null(colnames(X))) colnames(X) <- paste0("x",1:p)

#	Check input values for dispersion
	if(any(dispersion<0)) stop("dispersion values must be non-negative")
	phi <- rep_len(dispersion,n)

#	Check input values for offset
	if(!all(is.finite(offset))) stop("All offset values must be finite and non-missing")
	offset <- rep_len(offset,n)

#	Check input values for weights
	if(is.null(weights)) weights <- rep_len(1,n)
	if(any(weights <= 0)) stop("All weights must be positive")

#	Handle y all zero as special case
	if(ymax==0) {
#		Does X include an intercept term?
		if(colnames(X)[1]=="(Intercept)") {
			beta <- rep_len(0,p)
			names(beta) <- colnames(X)
			beta[1] <- -Inf
			mu <- rep.int(0,n)
			names(mu) <- rownames(X)
			return(list(coefficients=beta,fitted.values=mu,deviance=0,iter=0L,convergence="converged"))
		}
#		Does X span the intercept term, at least closely enough to preserve signs?
		One <- rep_len(1,n)
		fit <- .lm.fit(X,One)
		if(max(abs(fit$residuals)) < 1) {
			beta <- -1e10 * fit$coefficients
			names(beta) <- colnames(X)
			mu <- rep_len(0,n)
			names(mu) <- rownames(X)
			return(list(coefficients=beta,fitted.values=mu,deviance=0,iter=0L,convergence="converged"))
		}
#		If X is far from spanning the intercept term, then
#		initialize the iteration by trying to cancel out the offsets
		if(is.null(coef.start)) {
			fit <- lm.wfit(x=X,y=offset,w=weights)
			coef.start <- -fit$coefficients
		}
	}

#	Starting values
	delta <- 1/6
	y1 <- pmax(y,delta)
	if(is.null(coef.start)) {
		start.method <- match.arg(start.method,c("log(y)","mean"))
		if(start.method=="log(y)") {
			fit <- lm.wfit(X,log(y1)-offset,weights)
			beta <- fit$coefficients
			mu <- exp(fit$fitted.values+offset)
		} else {
			N <- exp(offset)
			rate <- y/N
			w <- weights*N/(1+phi*N)
			beta.mean <- log(sum(w*rate)/sum(w))
			beta <- qr.coef(qr(X),rep_len(beta.mean,n))
			mu <- drop(exp(X %*% beta + offset))
		}
	} else {
		beta <- coef.start
		mu <- drop(exp(X %*% beta + offset))
	}

	unit.dev.poissonlimit <- function(y,mu,phi) {
		b <- y-mu
		b2 <- 0.5*b^2*phi*(1+phi*(2/3*b-y))
        2 * ( y*log(y/mu) - b - b2 )
	}
	unit.dev.gamma <- function(y,mu) {
		2 * ( (y-mu)/mu - log(y/mu))
	}
	unit.dev.negbin <- function(y,mu,phi) {
		2 * ( y*log(y/mu) - (y+1/phi)*log((1+y*phi)/(1+mu*phi)) )
	}
	total.deviance <- function(y,mu,phi,w) {
		if(any(is.infinite(mu))) return(Inf)
		poisson.like <- (phi < 1e-4)
		gamma.like <- (phi*mu > 1e6)
		negbin <- !(poisson.like | gamma.like)
		y <- y+1e-8
		mu <- mu+1e-8
		unit.dev <- y
		if(any(poisson.like)) unit.dev[poisson.like] <- unit.dev.poissonlimit(y[poisson.like],mu[poisson.like],phi[poisson.like])
		if(any(gamma.like)) {
			m <- mu[gamma.like]
			alpha <- m/(1+phi[gamma.like]*m)
			unit.dev[gamma.like] <- unit.dev.gamma(y[gamma.like],m)*alpha
		}
		if(any(negbin)) unit.dev[negbin] <- unit.dev.negbin(y[negbin],mu[negbin],phi[negbin])
		sum(w*unit.dev)
	}

	dev <- total.deviance(y,mu,phi,weights)

#	Scoring iteration with Levenberg damping
	iter <- 0
	if(trace) cat("Iter =",iter,", Dev =",dev," Beta",beta,"\n")
	repeat {
		iter <- iter+1

#		test for iteration limit
		if(iter > maxit) break

#		information matrix
		v.div.mu <- 1+phi*mu
		XVX <- crossprod(X,(weights*mu/v.div.mu)*X)
		maxinfo <- max(diag(XVX))
		if(iter==1) {
			lambda <- maxinfo * 1e-6
			lambda <- max(lambda,1e-13)
			lambdaceiling <- maxinfo * 1e13
			lambdabig <- FALSE
			I <- diag(p)
		}

#		score vector
		dl <- crossprod(X,weights*(y-mu)/v.div.mu)

#		Levenberg damping
		dbeta <- beta
		lev <- 0
		repeat {
			lev <- lev+1

#			trial step
			R <- chol(XVX + lambda*I, pivot=TRUE)
			while(attr(R,"rank")<p) {
				lambda <- lambda*10
				R <- chol(XVX + lambda*I, pivot=TRUE)
			}
			j <- attr(R,"pivot")
			dbeta[j] <- backsolve(R,backsolve(R,dl[j],transpose=TRUE))
			betanew <- beta + dbeta
			munew <- drop(exp(X %*% betanew + offset))
			devnew <- total.deviance(y,munew,phi,weights)
			if(devnew <= dev || devnew/ymax < 1e-13) {
				beta <- betanew
				mu <- munew
				dev <- devnew
				break
			}

#			step not successful so increase damping
			lambda <- 2*lambda
			if(trace) cat("Damping increased to",lambda,"\n")

#			exit if too much damping
			lambdabig <- (lambda > lambdaceiling)
			if(lambdabig) {
				warning("Too much damping - convergence tolerance not achievable")
				break
			}
		}

#		iteration output
		if(trace) cat("Iter =",iter,", Dev =",dev," Beta",beta,"\n")

#		keep exiting if too much damping
		if(lambdabig) break

#		test for convergence
		scoresquare <- crossprod(dl,dbeta)
		if(trace) cat("Convergence criterion",scoresquare,dl,dbeta,"\n")
		if( scoresquare < tol || dev/ymax < 1e-12) break

#		decrease damping if successful at first try
		if(lev==1) lambda <- lambda/10
	}

	beta <- drop(beta)
	names(beta) <- colnames(X)
	convergence <- "converged"
	if(lambdabig) convergence <- "lambdabig"
	if(iter>maxit) convergence <- "maxit"
	list(coefficients=beta,fitted.values=as.vector(mu),deviance=dev,iter=iter,convergence=convergence)
}

