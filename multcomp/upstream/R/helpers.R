
# $Id$

### model.matrix.coxph doesn't return contrasts etc.
#model.matrix.coxph <- function(object, ...) {
#    mm <- model.matrix(delete.response(terms(object)),
#                       data = model.frame(object))
#    at <- attributes(mm)
#    mm <- mm[,-1]
#    at$dim[2] <- at$dim[2] - 1
#    at$dimnames[[2]] <- at$dimnames[[2]][-1]
#    at$assign <- at$assign[-1]
#    attributes(mm) <- at
#    mm
#}

model.matrix.coxph.penal <- function(object, ...) {

    class(object) <- "coxph"
    mm <- model.matrix(object)
    at <- attributes(mm)
    indx <- grep("frailty", colnames(mm))
    ret <- mm[ , -indx, drop = FALSE]
    attr(ret, "assign") <- at$assign[-indx]
    attr(ret, "contrasts") <- at$contrasts
    ret
}

model.frame.coxph.penal <- function(formula, ...) {

    object <- formula
    tm <- terms(object)
    class(object) <- "coxph"
    mf <- model.frame(object)
    ret <- cbind(mf[[1]], model.frame(delete.response(tm), data = mf))
    colnames(ret)[1] <- colnames(mf)[1]
    ret
}

terms.coxph.penal <- function(x, ...) {

    class(x) <- "coxph"           
    tm <- terms(x)
    ctm <- as.character(tm)
    x <- strsplit(ctm[3], "+", fixed = TRUE)[[1]]
    x <- x[-grep("frailty", x)]
    fm <- paste(ctm[2], "~", paste(x, collapse = "+"))
    terms(as.formula(fm))
}

coxph.penalcoef <- function(object, ...) {

    mm <- model.matrix(object)
    class(object) <- "coxph"
    cf <- coef(object)
    cf[1:ncol(mm)]
}

coxph.penalvcov <- function(object, ...) {
    
    mm <- model.matrix(object)
    class(object) <- "coxph"
    vc <- vcov(object)        
    vc[1:ncol(mm), 1:ncol(mm), drop = FALSE]
}


#model.matrix.survreg <- function(object, ...) {
#   model.matrix(delete.response(terms(object)),
#                       data = model.frame(object))
#}

### coxme objects
model.matrix.coxme <- function(object, ...) {
    class(object) <- "coxph"
    model.matrix(object)
}

### coxme objects
model.frame.coxme <- function(formula, ...) {
    object <- formula
    class(object) <- "coxph"
    model.frame(object)
}

model.matrix.aovlist <- function(object, ...)
    stop(sQuote("glht"), " does not support objects of class ", 
         sQuote("aovlist"))

model.matrix.lme <- function(object, ...)
    model.matrix(terms(object), data = model.frame(object), ...)

model.frame.lme <- function(formula, ...) {
    object <- formula
    ret <- object$data
    if (is.null(ret)) stop("object does not contain any data")
    ret
}

### extract coefficients, covariance matrix and 
### degrees of freedom (if available) from `model'
modelparm <- function(model, coef., vcov., df, ...) 
    UseMethod("modelparm")

modelparm.default <- function(model, coef. = coef, 
                              vcov. = function(x) vcov(x, complete = FALSE), 
                              df = NULL, ...) 
{

    ### allow specification of coef and vcov directly
    if (!is.function(coef.)) {
        beta <- coef.
        coef. <- function(model) return(beta)
    }
    if (!is.function(vcov.)) {
        sigma <- vcov.
        vcov. <- function(model) return(sigma)
    }

    ### extract coefficients and their covariance matrix
    beta <- try(coef.(model, ...))
    if (inherits(beta, "try-error"))
        stop("no ", sQuote("coef"), " method for ",
             sQuote("model"), " found!")

    sigma <- try(vcov.(model, ...))
    if (inherits(sigma, "try-error"))
        stop("no ", sQuote("vcov"), " method for ",
             sQuote("model"), " found!")       
    sigma <- as.matrix(sigma)

    if (any(length(beta) != dim(sigma))) 
        stop("dimensions of coefficients and covariance matrix don't match")

    ### determine degrees of freedom
    if (is.null(df)) {
        df <- 0
        ### check if a linear model was supplied
        if (class(model)[1] %in% c("aov", "lm")) {
            class(model) <- "lm"
            df <- summary(model)$df[2]
        }
        if (inherits(model, "parm"))
            df <- model$df
    } else {
        if (df < 0) stop(sQuote("df"), " is not positive")
    }

    ### try to identify non-estimable coefficients
    ### coef.aov removes NAs, thus touch coefficients 
    ### directly
    ocoef <- coef.(model, ...)
    if (inherits(model, "aov")) ocoef <- model$coefficients
    estimable <- rep(TRUE, length(ocoef))
    if (any(is.na(ocoef))) {
        estimable[is.na(ocoef)] <- FALSE
        beta <- ocoef[estimable]
    }

    ### just in case...
    if (length(beta) != ncol(sigma) || nrow(sigma) != sum(estimable))
        stop("could not extract coefficients and covariance matrix from ", 
             sQuote("model"))

    RET <- list(coef = beta, vcov = sigma, df = df, estimable = estimable)
    class(RET) <- "modelparm"
    RET
}

### mixed effects models (package `lme4')
modelparm.mer <- function(model, coef. = lme4::fixef, vcov. = vcov, df = NULL, ...)
    modelparm.default(model, coef. = coef., vcov. = vcov., df = df, ...)

### mixed effects models (package `lme4Eigen')
modelparm.merMod <- function(model, coef. = lme4::fixef, vcov. = vcov, df = NULL, ...)
    modelparm.default(model, coef. = coef., vcov. = vcov., df = df, ...)

### package `nlme'
modelparm.lme <- function(model, coef. = nlme::fixef, vcov. = vcov, df = NULL, ...)
    modelparm.default(model, coef. = coef., vcov. = vcov., df = df, ...)

### package glmmTMB
modelparm.glmmTMB <- function(model, 
    coef. = function(object) glmmTMB::fixef(object)$cond, 
    vcov. = function(object) vcov(object)$cond, 
    df = NULL, ...)
    modelparm.default(model, coef. = coef., vcov. = vcov., df = df, ...)

### survreg models (package `survival')
vcovsurvreg <- function(object, ...) {
    sigma <- vcov(object)
    p <- length(coef(object))
    return(sigma[1:p, 1:p])
}

### nlme:::gls
model.matrix.gls <- function(object, ...)
    model.matrix(terms(object), data = nlme::getData(object), ...)

model.frame.gls <- function(formula, ...)
    model.frame(formula(formula), data = nlme::getData(formula), ...)

terms.gls <- function(x, ...)
    terms(model.frame(x), ...)

modelparm.survreg <- function(model, coef. = coef, vcov. = vcovsurvreg, df = NULL, ...)
    modelparm.default(model, coef. = coef., vcov. = vcov., df = df, ...)

modelparm.aovlist <- function(model, coef. = coef, vcov. = vcov, df = NULL, ...)
    stop(sQuote("glht"), " does not support objects of class ", sQuote("aovlist"))

modelparm.coxme <- function(model, coef. = coef, vcov. = vcov, df = NULL, ...)
    modelparm.default(model, coef. = coef., vcov. = vcov., df = df, ...)

modelparm.coxph.penal <- function(model, coef. = coxph.penalcoef, 
                                  vcov. = coxph.penalvcov, df = NULL, ...)
    modelparm.default(model, coef. = coef., vcov. = vcov., df = df, ...)

model.matrix.polr <- function(object, ...) {
    mm <- model.matrix(delete.response(terms(object)),
                      data = model.frame(object))
    at <- attributes(mm)
    mm <- mm[,-1]
    at$dim[2] <- at$dim[2] - 1
    at$dimnames[[2]] <- at$dimnames[[2]][-1]
    at$assign <- at$assign[-1]
    attributes(mm) <- at
    mm
}


polrvcov <- function(object) {
   cf <- coef(object)
   vcov <- vcov(object)
   vcov[names(cf), names(cf)]
}

modelparm.polr <- function(model, coef. = coef, vcov. = polrvcov, df = NULL, ...)
    modelparm.default(model, coef. = coef., vcov. = vcov., df = df, ...)

### fixed effects models (package fixest). Contributed by Grant McDermott 2021-12-17
modelparm.fixest <- function(model, coef. = coef, vcov. = vcov, df = NULL, ...) {
    model <- summary(model, vcov = vcov.)
    vcov. <- vcov(model)
    if (is.null(df))
        df <- fixest::degrees_freedom(model, type = "resid")
    modelparm.default(model, coef. = coef., vcov. = vcov., df = df, ...)
}

### gamlss (donated by Marcio A Diniz)
model.matrix.gamlss <- function(object, ...) {
    cf <- na.exclude(coef(object))
    
    ### extract model matrix, frame and terms
    mm <- model.matrix(terms(object),
                       data = model.frame(object))
    aux <- list(assign = attributes(mm)[["assign"]],
                contrasts = attributes(mm)[["contrasts"]])
    
    mm <- mm[, 1:length(cf)]
    attr(mm, "assign") <- aux$assign[1:length(cf)]
    attr(mm, "contrasts") <- aux$contrasts
    mm
}


gamlss.coef <- function(object, ...) {
    dots <- list(...)
    #class(object) <- class(object)[1]
    cf <- na.exclude(coef(object, what = dots$what))
    cf
}

gamlss.vcov <- function(object, ...) {
    dots <- list(...)
    #class(object) <- class(object)[1]
    
    p <- match(dots$what, object$parameters)
    
    vc <- vcov(object, what = dots$what) 
    index <- which(cumsum(rownames(vc) == "(Intercept)") == p)
    vc[index, index, drop = FALSE]
}

### modified from package MASS  
MPinv <- function (X, tol = sqrt(.Machine$double.eps))
{
    if (length(dim(X)) > 2 || !(is.numeric(X) || is.complex(X)))
        stop("X must be a numeric or complex matrix")
    if (!is.matrix(X))
        X <- as.matrix(X)
    Xsvd <- svd(X)
    if (is.complex(X))
        Xsvd$u <- Conj(Xsvd$u)
    Positive <- Xsvd$d > max(tol * Xsvd$d[1], 0)
    if (all(Positive)) 
        RET <- Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))   
    else if (!any(Positive))
        RET <- array(0, dim(X)[2:1])
    else RET <- Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) *
        t(Xsvd$u[, Positive, drop = FALSE]))
    return(list(MPinv = RET, rank = sum(Positive)))
}

### meaningless ...
chkdots <- function(...) {

    lst <- list(...)
    if (length(lst) > 0) {
        warning("Argument(s) ", sQuote(names(lst)), " passed to ", sQuote("..."), 
                " are ignored", call. = TRUE)
    }
}
