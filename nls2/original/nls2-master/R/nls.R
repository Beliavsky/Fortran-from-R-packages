#  File src/library/stats/R/nls.R
#  Part of the R package, https://www.R-project.org
#
#  Copyright (C) 2000-2020 The R Core Team
#  Copyright (C) 1999-1999 Saikat DebRoy, Douglas M. Bates, Jose C. Pinheiro
#
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
#  https://www.R-project.org/Licenses/

###
###            Nonlinear least squares for R
###

nlsModel.plinear <- function(form, data, start, wts)
{
    thisEnv <- environment()
    env <- new.env(hash = TRUE, parent=environment(form))
    for(i in names(data)) assign(i, data[[i]], envir = env)
    ind <- as.list(start)
    p2 <- 0
    for(i in names(ind)) {
        temp <- start[[i]]
        storage.mode(temp) <- "double"
        assign(i, temp, envir = env)
        ind[[i]] <- p2 + seq_along(start[[i]])
        p2 <- p2 + length(start[[i]])
    }
    lhs <- eval(form[[2L]], envir = env)
    storage.mode(lhs) <- "double"
    rhs <- eval(form[[3L]], envir = env)
    storage.mode(rhs) <- "double"
    .swts <- if(!missing(wts) && length(wts))
        sqrt(wts) else 1 # more efficient than  rep_len(1, NROW(rhs))
    assign(".swts", .swts, envir = env)
    p1 <- NCOL(rhs)
    p <- p1 + p2
    n <- length(lhs)
    fac <- (n -  p)/p
    cc <- QR.B <- NA
    useParams <- rep_len(TRUE, p2)
    if(is.null(attr(rhs, "gradient"))) {
        getRHS.noVarying <- function()
            numericDeriv(form[[3L]], names(ind), env)
        getRHS <- getRHS.noVarying
        rhs <- getRHS()
    } else {
        getRHS.noVarying <- function() eval(form[[3L]], envir = env)
        getRHS <- getRHS.noVarying
    }
    dimGrad <- dim(attr(rhs, "gradient"))
    marg <- length(dimGrad)
    if(marg > 0) {
        gradSetArgs <- vector("list", marg + 1L)
        for(i in 2:marg)
            gradSetArgs[[i]] <- rep_len(TRUE, dimGrad[i-1])
        useParams <- rep_len(TRUE, dimGrad[marg])
    } else {
        gradSetArgs <- vector("list", 2L)
        useParams <- rep_len(TRUE, length(attr(rhs, "gradient")))
    }
    gradSetArgs[[1L]] <- (~attr(ans, "gradient"))[[2L]]
    gradCall <-
        switch(length(gradSetArgs) - 1L,
               call("[", gradSetArgs[[1L]], gradSetArgs[[2L]]),
               call("[", gradSetArgs[[1L]], gradSetArgs[[2L]], gradSetArgs[[2L]]),
               call("[", gradSetArgs[[1L]], gradSetArgs[[2L]], gradSetArgs[[2L]],
                    gradSetArgs[[3L]]),
               call("[", gradSetArgs[[1L]], gradSetArgs[[2L]], gradSetArgs[[2L]],
                    gradSetArgs[[3L]], gradSetArgs[[4L]]))
    getRHS.varying <- function()
    {
        ans <- getRHS.noVarying()
        attr(ans, "gradient") <- eval(gradCall)
        ans
    }
    QR.rhs <- qr(.swts * rhs)
    lin <- qr.coef(QR.rhs, .swts * lhs)
    resid <- qr.resid(QR.rhs, .swts * lhs)
    topzero <- double(p1)
    dev <- sum(resid^2)
    if(marg <= 1) {
        ddot <- function(A, b) A %*% b
        dtdot <- function(A, b) t(A) %*% b
    } else if(marg == 2) {
        if(p1 == 1) {
            ddot <- function(A, b) as.matrix(A*b)
            dtdot <- function(A, b) t(b) %*% A
        } else if(p2 == 1) {
            ddot <- function(A, b) A %*% b
            dtdot <- function(A, b) t(A) %*% b
        }
    } else {
        ddot <- function(A, b) apply(A, MARGIN = 3L, FUN="%*%", b)
        dtdot <- function(A, b) apply(A, MARGIN = c(2L,3L), FUN = "%*%", b)
    }

    getPars.noVarying <- function() unlist(mget(names(ind), env))
    getPars.varying   <- function() unlist(mget(names(ind), env))[useParams]
    getPars <- getPars.noVarying

    internalPars <- getPars()
    setPars.noVarying <- function(newPars)
    {
        assign("internalPars", newPars, envir = thisEnv)
        for(i in names(ind)) {
            assign(i, unname(newPars[ ind[[i]] ]), envir = env )
        }
    }
    setPars.varying <- function(newPars)
    {
        internalPars[useParams] <- newPars
        for(i in names(ind))
            assign(i, unname(internalPars[ ind[[i]] ]), envir = env)
    }
    setPars <- setPars.noVarying
    getPred <-
        if(is.matrix(rhs)) function(X) as.numeric(X %*% lin)
        else function(X) X * lin

    m <-
        list(resid = function() resid,
             fitted = function() getPred(rhs),
             formula = function() form,
             deviance = function() dev,
             lhs = function() lhs,
             gradient = function() attr(rhs, "gradient"),
             conv = function() {
                 assign("cc", c(topzero, qr.qty(QR.rhs, .swts * lhs)[ -(1L:p1)]),
                        envir = thisEnv)
                 rr <- qr.qy(QR.rhs, cc)
                 B <- qr.qty(QR.rhs, .swts * ddot(attr(rhs, "gradient"), lin))
                 B[1L:p1, ] <- dtdot(.swts * attr(rhs, "gradient"), rr)
                 R <- t( qr.R(QR.rhs)[1L:p1, ] )
                 if(p1 == 1) B[1, ] <- B[1, ]/ c(R)
                 else B[1L:p1, ] <- forwardsolve(R, B[1L:p1, ])
                 assign("QR.B", qr(B), envir = thisEnv)
                 rr <- qr.qty(QR.B, cc)
                 sqrt( fac*sum(rr[1L:p1]^2) / sum(rr[-(1L:p1)]^2) )
             },
             incr = function() qr.solve(QR.B, cc),
             setVarying = function(vary = rep_len(TRUE, length(useParams))) {
                 assign("useParams", if(is.character(vary)) {
                     temp <- logical(length(useParams))
                     temp[unlist(ind[vary])] <- TRUE
                     temp
                 } else if(is.logical(vary) && length(vary) != length(useParams))
                        stop("setVarying : 'vary' length must match length of parameters")
                 else {
                     vary
                 }, envir = thisEnv)
                 gradCall[[length(gradCall)]] <<- useParams
                 if(all(useParams)) {
                     assign("setPars", setPars.noVarying, envir = thisEnv)
                     assign("getPars", getPars.noVarying, envir = thisEnv)
                     assign("getRHS", getRHS.noVarying, envir = thisEnv)
                 } else {
                     assign("setPars", setPars.varying, envir = thisEnv)
                     assign("getPars", getPars.varying, envir = thisEnv)
                     assign("getRHS", getRHS.varying, envir = thisEnv)
                 }
             },
             setPars = function(newPars) {
                 setPars(newPars)
                 assign("QR.rhs",
                        qr(.swts * assign("rhs", getRHS(), envir = thisEnv)),
                        envir = thisEnv)
                 assign("resid", qr.resid(QR.rhs, .swts * lhs),
                        envir = thisEnv)
                 assign("dev", sum(resid^2), envir = thisEnv )
                 if(QR.rhs$rank < p1) {
                     return(1)
                 } else {
                     assign("lin", qr.coef(QR.rhs, .swts * lhs),
                            envir = thisEnv)
                     return(0)
                 }
             },
             getPars = function() getPars(),
             getAllPars = function() c( getPars(), c( .lin = lin ) ),
             getEnv = function() env,
             trace = function() {
                 cat(format(dev),": ", format(c(getPars(), lin)))
                 cat("\n")
             },
             Rmat = function()
             qr.R(qr(.swts * cbind(ddot(attr(rhs, "gradient"), lin), rhs))),
             predict = function(newdata = list(), qr = FALSE)
             getPred(eval(form[[3L]], as.list(newdata), env))
             )
    class(m) <- c("nlsModel.plinear", "nlsModel")
    m$conv()
    on.exit( remove( data, i, m, marg, n, p, start, temp, gradSetArgs) )
    m
}

nlsModel <- function(form, data, start, wts, upper=NULL)
{
    thisEnv <- environment()
    env <- new.env(hash = TRUE, parent = environment(form))
    for(i in names(data)) assign(i, data[[i]], envir = env)
    ind <- as.list(start)
    parLength <- 0
    for(i in names(ind) ) {
        temp <- start[[i]]
        storage.mode(temp) <- "double"
        assign(i, temp, envir = env)
        ind[[i]] <- parLength + seq_along(start[[i]])
        parLength <- parLength + length(start[[i]])
    }
    getPars.noVarying <- function() unlist(mget(names(ind), env))
    getPars <- getPars.noVarying
    internalPars <- getPars()

    if(!is.null(upper)) upper <- rep_len(upper, parLength)
    useParams <- rep_len(TRUE, parLength)
    lhs <- eval(form[[2L]], envir = env)
    rhs <- eval(form[[3L]], envir = env)
    .swts <- if(!missing(wts) && length(wts))
        sqrt(wts) else rep_len(1, length(rhs))
    assign(".swts", .swts, envir = env)
    resid <- .swts * (lhs - rhs)
    dev <- sum(resid^2)
    if(is.null(attr(rhs, "gradient"))) {
        getRHS.noVarying <- function() {
            if(is.null(upper))
                numericDeriv(form[[3L]], names(ind), env)
            else
                numericDeriv(form[[3L]], names(ind), env,
                             ifelse(internalPars < upper, 1, -1))
        }
        getRHS <- getRHS.noVarying
        rhs <- getRHS()
    } else {
        getRHS.noVarying <- function() eval(form[[3L]], envir = env)
        getRHS <- getRHS.noVarying
    }
    dimGrad <- dim(attr(rhs, "gradient"))
    marg <- length(dimGrad)
    if(marg > 0L) {
        gradSetArgs <- vector("list", marg + 1L)
        for(i in 2L:marg)
            gradSetArgs[[i]] <- rep_len(TRUE, dimGrad[i-1])
        useParams <- rep_len(TRUE, dimGrad[marg])
    } else {
        gradSetArgs <- vector("list", 2L)
        useParams <- rep_len(TRUE, length(attr(rhs, "gradient")))
    }
    npar <- length(useParams)
    gradSetArgs[[1L]] <- (~attr(ans, "gradient"))[[2L]]
    gradCall <-
        switch(length(gradSetArgs) - 1L,
               call("[", gradSetArgs[[1L]], gradSetArgs[[2L]], drop = FALSE),
               call("[", gradSetArgs[[1L]], gradSetArgs[[2L]], gradSetArgs[[2L]],
                    drop = FALSE),
               call("[", gradSetArgs[[1L]], gradSetArgs[[2L]], gradSetArgs[[2L]],
                    gradSetArgs[[3L]], drop = FALSE),
               call("[", gradSetArgs[[1L]], gradSetArgs[[2L]], gradSetArgs[[2L]],
                    gradSetArgs[[3L]], gradSetArgs[[4L]], drop = FALSE))
    getRHS.varying <- function()
    {
        ans <- getRHS.noVarying()
        attr(ans, "gradient") <- eval(gradCall)
        ans
    }
    if(length(gr <- attr(rhs, "gradient")) == 1L)
		    attr(rhs, "gradient") <- gr <- as.vector(gr)
    QR <- qr(.swts * gr)
    qrDim <- min(dim(QR$qr))
    if(QR$rank < qrDim)
        stop("singular gradient matrix at initial parameter estimates")

    getPars.varying <- function() unlist(mget(names(ind), env))[useParams]
    setPars.noVarying <- function(newPars)
    {
        assign("internalPars", newPars, envir = thisEnv)
        for(i in names(ind))
            assign(i, unname(newPars[ ind[[i]] ]), envir = env)
    }
    setPars.varying <- function(newPars)
    {
        internalPars[useParams] <- newPars
        for(i in names(ind))
            assign(i, unname(internalPars[ ind[[i]] ]), envir = env)
    }
    setPars <- setPars.noVarying

    on.exit(remove(i, data, parLength, start, temp, m))
    ## must use weighted resid for use with "port" algorithm.
    m <-
	list(resid = function() resid,
	     fitted = function() rhs,
	     formula = function() form,
	     deviance = function() dev,
	     lhs = function() lhs,
	     gradient = function() .swts * attr(rhs, "gradient"),
	     conv = function() {
		 if(npar == 0) return(0)
		 rr <- qr.qty(QR, resid) # rotated residual vector
		 sqrt( sum(rr[1L:npar]^2) / sum(rr[-(1L:npar)]^2))
	     },
	     incr = function() qr.coef(QR, resid),
	     setVarying = function(vary = rep_len(TRUE, length(useParams))) {
		 assign("useParams",
			if(is.character(vary)) {
			    temp <- logical(length(useParams))
			    temp[unlist(ind[vary])] <- TRUE
			    temp
			} else if(is.logical(vary) &&
				  length(vary) != length(useParams))
			stop("setVarying : 'vary' length must match length of parameters")
			else {
			    vary
			}, envir = thisEnv)
		 gradCall[[length(gradCall) - 1L]] <<- useParams
		 if(all(useParams)) {
		     assign("setPars", setPars.noVarying, envir = thisEnv)
		     assign("getPars", getPars.noVarying, envir = thisEnv)
		     assign("getRHS", getRHS.noVarying, envir = thisEnv)
		     assign("npar", length(useParams), envir = thisEnv)
		 } else {
		     assign("setPars", setPars.varying, envir = thisEnv)
		     assign("getPars", getPars.varying, envir = thisEnv)
		     assign("getRHS", getRHS.varying, envir = thisEnv)
                     ## FIXME this is which(useParams)
		     assign("npar", length(seq_along(useParams)[useParams]),
			    envir = thisEnv)
		 }
	     },
	     setPars = function(newPars) {
		 setPars(newPars)
		 assign("resid", .swts *
			(lhs - assign("rhs", getRHS(), envir = thisEnv)),
			envir = thisEnv)
		 assign("dev", sum(resid^2), envir = thisEnv)
		 if(length(gr <- attr(rhs, "gradient")) == 1L) gr <- c(gr)
		 assign("QR", qr(.swts * gr), envir = thisEnv )
		 (QR$rank < min(dim(QR$qr))) # to catch the singular gradient matrix
	     },
	     getPars = function() getPars(),
	     getAllPars = function() getPars(),
	     getEnv = function() env,
	     trace = function() {
                 cat(format(dev),": ", format(getPars()))
                 cat("\n")
             },
	     Rmat = function() qr.R(QR),
	     predict = function(newdata = list(), qr = FALSE)
	     eval(form[[3L]], as.list(newdata), env)
	     )
    class(m) <- "nlsModel"
    m
}


coef.nls <- function(object, ...) object$m$getAllPars()

summary.nls <-
    function (object, correlation = FALSE, symbolic.cor = FALSE, ...)
{
    r <- as.vector(object$m$resid()) # These are weighted residuals.
    w <- object$weights
    n <- if (!is.null(w)) sum(w > 0) else length(r)
    param <- coef(object)
    pnames <- names(param)
    p <- length(param)
    rdf <- n - p
    resvar <- if(rdf <= 0) NaN else deviance(object)/rdf
    XtXinv <- chol2inv(object$m$Rmat())
    dimnames(XtXinv) <- list(pnames, pnames)
    se <- sqrt(diag(XtXinv) * resvar)
    tval <- param/se
    param <- cbind(param, se, tval, 2 * pt(abs(tval), rdf, lower.tail = FALSE))
    dimnames(param) <-
        list(pnames, c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
    ans <- list(formula = formula(object), residuals = r, sigma = sqrt(resvar),
                df = c(p, rdf), cov.unscaled = XtXinv,
                call = object$call,
                convInfo = object$convInfo,
                control = object$control,
                na.action = object$na.action,
                coefficients = param,
                parameters = param)# never documented, for back-compatibility
    if(correlation && rdf > 0) {
        ans$correlation <- (XtXinv * resvar)/outer(se, se)
        ans$symbolic.cor <- symbolic.cor
    }
    ## if(identical(object$call$algorithm, "port"))
    ##     ans$message <- object$message
    class(ans) <- "summary.nls"
    ans
}

.p.nls.convInfo <- function(x, digits,
			    show. = getOption("show.nls.convergence", TRUE))
{
    if(!is.null(x$convInfo)) # older fits will not have this
        with(x$convInfo,
         {
             if(identical(x$call$algorithm, "port"))
                 cat("\nAlgorithm \"port\", convergence message: ",
                     stopMessage, "\n", sep = "")
             else {
                 if(!isConv || show.) {
                     cat("\nNumber of iterations",
                         if(isConv) "to convergence:" else "till stop:", finIter,
                         "\nAchieved convergence tolerance:",
                         format(finTol, digits = digits))
                     cat("\n")
                 }
                 if(!isConv) {
                     cat("Reason stopped:", stopMessage)
                     cat("\n")
                 }
             }
         })

    invisible()
}

print.nls <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
    cat("Nonlinear regression model\n")
    cat("  model: ", deparse(formula(x)), "\n", sep = "")
    cat("   data: ", deparse(x$data), "\n", sep = "")
    print(x$m$getAllPars(), digits = digits, ...)
    cat(" ", if(!is.null(x$weights) && diff(range(x$weights))) "weighted ",
	"residual sum-of-squares: ", format(x$m$deviance(), digits = digits),
	"\n", sep = "")
    .p.nls.convInfo(x, digits = digits)
    invisible(x)
}

print.summary.nls <-
  function (x, digits = max(3L, getOption("digits") - 3L),
            symbolic.cor = x$symbolic.cor,
            signif.stars = getOption("show.signif.stars"), ...)
{
    cat("\nFormula: ",
	paste(deparse(x$formula), sep = "\n", collapse = "\n"),
        "\n", sep = "")
    df <- x$df
    rdf <- df[2L]
    cat("\nParameters:\n")
    printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars,
                 ...)
    cat("\nResidual standard error:",
        format(signif(x$sigma, digits)), "on", rdf, "degrees of freedom")
    cat("\n")
    correl <- x$correlation
    if (!is.null(correl)) {
        p <- NCOL(correl)
        if (p > 1) {
            cat("\nCorrelation of Parameter Estimates:\n")
	    if(is.logical(symbolic.cor) && symbolic.cor) {
		print(symnum(correl, abbr.colnames = NULL))
            } else {
                correl <- format(round(correl, 2), nsmall = 2L, digits = digits)
                correl[!lower.tri(correl)] <- ""
                print(correl[-1, -p, drop=FALSE], quote = FALSE)
            }
        }
    }

    .p.nls.convInfo(x, digits = digits)

    if(nzchar(mess <- naprint(x$na.action))) cat("  (", mess, ")\n", sep = "")
    cat("\n")
    invisible(x)
}

weights.nls <- function(object, ...) object$weights

predict.nls <-
  function(object, newdata, se.fit = FALSE, scale = NULL, df = Inf,
           interval = c("none", "confidence", "prediction"), level = 0.95,
           ...)
{
    if (missing(newdata)) return(as.vector(fitted(object)))
    if(!is.null(cl <- object$dataClasses)) .checkMFClasses(cl, newdata)
    object$m$predict(newdata)
}

fitted.nls <- function(object, ...)
{
    val <- as.vector(object$m$fitted())
    if(!is.null(object$na.action)) val <- napredict(object$na.action, val)
    lab <- "Fitted values"
    if (!is.null(aux <- attr(object, "units")$y)) lab <- paste(lab, aux)
    attr(val, "label") <- lab
    val
}

formula.nls <- function(x, ...) x$m$formula()

residuals.nls <- function(object, type = c("response", "pearson"), ...)
{
    type <- match.arg(type)
    if (type == "pearson") {
        val <- as.vector(object$m$resid())
        std <- sqrt(sum(val^2)/(length(val) - length(coef(object))))
        val <- val/std
        if(!is.null(object$na.action)) val <- naresid(object$na.action, val)
        attr(val, "label") <- "Standardized residuals"
    } else {
        val <- as.vector(object$m$lhs() - object$m$fitted())
        if(!is.null(object$na.action))
            val <- naresid(object$na.action, val)
        lab <- "Residuals"
        if (!is.null(aux <- attr(object, "units")$y)) lab <- paste(lab, aux)
        attr(val, "label") <- lab
    }
    val
}

logLik.nls <- function(object, REML = FALSE, ...)
{
    if (REML)
        stop("cannot calculate REML log-likelihood for \"nls\" objects")
    res <- object$m$resid() # These are weighted residuals.
    N <- length(res)
    w <- object$weights %||% rep_len(1, N)
    ## Note the trick for zero weights
    zw <- w == 0
    N <- sum(!zw)
    val <-  -N * (log(2 * pi) + 1 - log(N) - sum(log(w + zw))/N + log(sum(res^2)))/2
    ## the formula here corresponds to estimating sigma^2.
    attr(val, "df") <- 1L + length(coef(object))
    attr(val, "nobs") <- attr(val, "nall") <- N
    class(val) <- "logLik"
    val
}

df.residual.nls <- function(object, ...)
{
    w <- object$weights
    n <- if(!is.null(w)) sum(w != 0) else length(object$m$resid())
    n - length(coef(object))
}

deviance.nls <- function(object, ...) object$m$deviance()

vcov.nls <- function(object, ...)
{
    sm <- summary(object)
    sm$cov.unscaled * sm$sigma^2
}


anova.nls <- function(object, ...)
{
    if(length(list(object, ...)) > 1L) return(anovalist.nls(object, ...))
    stop("anova is only defined for sequences of \"nls\" objects")
}

anovalist.nls <- function (object, ..., test = NULL)
{
    objects <- list(object, ...)
    responses <- as.character(lapply(objects,
				     function(x) formula(x)[[2L]]))
    sameresp <- responses == responses[1L]
    if (!all(sameresp)) {
	objects <- objects[sameresp]
        warning(gettextf("models with response %s removed because response differs from model 1",
                         sQuote(deparse(responses[!sameresp]))),
                domain = NA)
    }
    ## calculate the number of models
    nmodels <- length(objects)
    if (nmodels == 1L)
        stop("'anova' is only defined for sequences of \"nls\" objects")

    models <- as.character(lapply(objects, function(x) formula(x)))

    ## extract statistics
    df.r <- unlist(lapply(objects, df.residual))
    ss.r <- unlist(lapply(objects, deviance))
    df <- c(NA, -diff(df.r))
    ss <- c(NA, -diff(ss.r))
    ms <- ss/df
    f <- p <- rep_len(NA_real_, nmodels)
    for(i in 2:nmodels) {
	if(df[i] > 0) {
	    f[i] <- ms[i]/(ss.r[i]/df.r[i])
	    p[i] <- pf(f[i], df[i], df.r[i], lower.tail = FALSE)
	}
	else if(df[i] < 0) {
	    f[i] <- ms[i]/(ss.r[i-1]/df.r[i-1])
	    p[i] <- pf(f[i], -df[i], df.r[i-1], lower.tail = FALSE)
	}
	else {                          # df[i] == 0
            ss[i] <- 0
	}
    }
    table <- data.frame(df.r,ss.r,df,ss,f,p)
    dimnames(table) <- list(1L:nmodels, c("Res.Df", "Res.Sum Sq", "Df",
					 "Sum Sq", "F value", "Pr(>F)"))
    ## construct table and title
    title <- "Analysis of Variance Table\n"
    topnote <- paste0("Model ", format(1L:nmodels), ": ", models,
                      collapse = "\n")

    ## calculate test statistic if needed
    structure(table, heading = c(title, topnote),
	      class = c("anova", "data.frame")) # was "tabular"
}

"%||%" <- function (L, R) if (is.null(L)) R else L

