#
#  Copyright (C) 2004-2016 Friedrich Leisch and Bettina Gruen
#  $Id: models.R 5079 2016-01-31 12:21:12Z gruen $
#

FLXMRglm <- function(formula=.~.,
                     family=c("gaussian", "binomial", "poisson", "Gamma"),
                     offset=NULL)
{
    family <- match.arg(family)
    glmrefit <- function(x, y, w) {
      fit <- c(glm.fit(x, y, weights=w, offset=offset,
                       family=get(family, mode="function")()),
               list(call = sys.call(), offset = offset,
                    control = eval(formals(glm.fit)$control),            
                    method = "weighted.glm.fit"))
      fit$df.null <- sum(w) + fit$df.null - fit$df.residual - fit$rank
      fit$df.residual <- sum(w) - fit$rank
      fit$x <- x
      fit
    }
                
    z <- new("FLXMRglm", weighted=TRUE, formula=formula,
             name=paste("FLXMRglm", family, sep=":"), offset = offset,
             family=family, refit=glmrefit)
    z@preproc.y <- function(x){
      if (ncol(x) > 1)
        stop(paste("for the", family, "family y must be univariate"))
      x
    }

    if(family=="gaussian"){
      z@defineComponent <- function(para) {
        predict <- function(x, ...) {
          dotarg = list(...)
          if("offset" %in% names(dotarg)) offset <- dotarg$offset
          p <- x %*% para$coef
          if (!is.null(offset)) p <-  p + offset
          p
        }

        logLik <- function(x, y, ...)
          dnorm(y, mean=predict(x, ...), sd=para$sigma, log=TRUE)

        new("FLXcomponent",
            parameters=list(coef=para$coef, sigma=para$sigma),
            logLik=logLik, predict=predict,
            df=para$df)
      }

      z@fit <- function(x, y, w, component){
          fit <- lm.wfit(x, y, w=w, offset=offset)
          z@defineComponent(para = list(coef = coef(fit), df = ncol(x)+1,
                                sigma =  sqrt(sum(fit$weights * fit$residuals^2 /
                                                      mean(fit$weights))/ (nrow(x)-fit$rank))))
      }
    }
    else if(family=="binomial"){
      z@preproc.y <- function(x){
        if (ncol(x) != 2)
          stop("for the binomial family, y must be a 2 column matrix\n",
               "where col 1 is no. successes and col 2 is no. failures")
        if (any(x < 0))
          stop("negative values are not allowed for the binomial family")
        x
      }     
      z@defineComponent <- function(para) {
        predict <- function(x, ...) {
          dotarg = list(...)
          if("offset" %in% names(dotarg)) offset <- dotarg$offset
          p <- x %*% para$coef
          if (!is.null(offset)) p <- p + offset
          get(family, mode = "function")()$linkinv(p)
        }
        logLik <- function(x, y, ...)
          dbinom(y[,1], size=rowSums(y), prob=predict(x, ...), log=TRUE)

        new("FLXcomponent",
            parameters=list(coef=para$coef),
            logLik=logLik, predict=predict,
            df=para$df)
      }

      z@fit <- function(x, y, w, component){
        fit <- glm.fit(x, y, weights=w, family=binomial(), offset=offset, start=component$coef)
        z@defineComponent(para = list(coef = coef(fit), df = ncol(x)))
      }
    }
    else if(family=="poisson"){
      z@defineComponent <- function(para) {
        predict <- function(x, ...) {
          dotarg = list(...)
          if("offset" %in% names(dotarg)) offset <- dotarg$offset
          p <- x %*% para$coef
          if (!is.null(offset)) p <- p + offset
          get(family, mode = "function")()$linkinv(p)
        }
        logLik <- function(x, y, ...)
          dpois(y, lambda=predict(x, ...), log=TRUE)
        
        new("FLXcomponent",
            parameters=list(coef=para$coef),
            logLik=logLik, predict=predict,
            df=para$df)
      }
          
      z@fit <- function(x, y, w, component){
        fit <- glm.fit(x, y, weights=w, family=poisson(), offset=offset, start=component$coef)
        z@defineComponent(para = list(coef = coef(fit), df = ncol(x)))
      }
    }
    else if(family=="Gamma"){
      z@defineComponent <- function(para) {
        predict <- function(x, ...) {
          dotarg = list(...)
          if("offset" %in% names(dotarg)) offset <- dotarg$offset
          p <- x %*% para$coef
          if (!is.null(offset)) p <- p + offset
          get(family, mode = "function")()$linkinv(p)
        }
        logLik <- function(x, y, ...)
          dgamma(y, shape = para$shape, scale=predict(x, ...)/para$shape, log=TRUE)
        
        new("FLXcomponent", 
            parameters = list(coef = para$coef, shape = para$shape),
            predict = predict, logLik = logLik,
            df = para$df)
      }

      z@fit <- function(x, y, w, component){
        fit <- glm.fit(x, y, weights=w, family=Gamma(), offset=offset, start=component$coef)
        z@defineComponent(para = list(coef = coef(fit), df = ncol(x)+1,
                              shape = sum(fit$prior.weights)/fit$deviance))
        
      }
    }
    else stop(paste("Unknown family", family))
    z
}

###**********************************************************

FLXMCmvnorm <- function(formula=.~., diagonal=TRUE)
{
    z <- new("FLXMC", weighted=TRUE, formula=formula,
             dist = "mvnorm", name="model-based Gaussian clustering")

    z@defineComponent <- function(para) {
      logLik <- function(x, y)
        mvtnorm::dmvnorm(y, mean=para$center, sigma=para$cov, log=TRUE)
    
      predict <-  function(x, ...)
        matrix(para$center, nrow=nrow(x), ncol=length(para$center),
               byrow=TRUE)
      new("FLXcomponent", parameters=list(center = para$center, cov = para$cov),
          df=para$df, logLik=logLik, predict=predict)
    }
    
    z@fit <- function(x, y, w, ...){
      para <- cov.wt(y, wt=w)[c("center","cov")]
      para$df <- (3*ncol(y) + ncol(y)^2)/2
      if(diagonal){
        para$cov <- diag(diag(para$cov))
        para$df <- 2*ncol(y)
      }
      z@defineComponent(para)
    }
    z
}

FLXMCnorm1 <- function(formula=.~.)
{
    z <- new("FLXMC", weighted=TRUE, formula=formula,
             dist = "mvnorm", name="model-based univariate Gaussian clustering")

    z@defineComponent <- function(para) {
      logLik <- function(x, y)
        dnorm(y, mean=para$center, sd=sqrt(para$cov), log=TRUE)
    
      predict <-  function(x, ...)
        matrix(para$center, nrow=nrow(x), ncol=1,
               byrow=TRUE)
      new("FLXcomponent",
          parameters=list(mean = as.vector(para$center), sd = as.vector(sqrt(para$cov))),
          df=para$df, logLik=logLik, predict=predict)
    }
    
    z@fit <- function(x, y, w, ...){
      para <- cov.wt(as.matrix(y), wt=w)[c("center","cov")]
      z@defineComponent(c(para, list(df = 2)))
    }
    z
}


###**********************************************************

FLXMCmvbinary <- function(formula=.~., truncated = FALSE) {
  if (truncated) return(MCmvbinary_truncated(formula))
  else return(MCmvbinary(formula))
}

MCmvbinary <- function(formula=.~.)
{
    z <- new("FLXMC", weighted=TRUE, formula=formula,
             dist = "mvbinary", name="model-based binary clustering")

    ## make sure that y is binary
    z@preproc.y <- function(x){
        storage.mode(x) <- "logical"
        storage.mode(x) <- "integer"
        x
    }
    z@defineComponent <- function(para) {
      predict <- function(x, ...){
        matrix(para$center, nrow=nrow(x), ncol=length(para$center),
               byrow=TRUE)
      }
        
      logLik <- function(x, y){
        p <- matrix(para$center, nrow=nrow(x), ncol=length(para$center),
                    byrow=TRUE)
        rowSums(log(y*p+(1-y)*(1-p)))
      }
            
      new("FLXcomponent", parameters=list(center=para$center), df=para$df,
          logLik=logLik, predict=predict)
    }

    z@fit <- function(x, y, w, ...)
       z@defineComponent(list(center = colSums(w*y)/sum(w), df = ncol(y)))
    z
}




###**********************************************************

binary_truncated <- function(y, w, maxit = 200, epsilon = .Machine$double.eps) {
  r_k <- colSums(y*w)/sum(w)
  r_0 <- 0
  llh.old <- -Inf
  for (i in seq_len(maxit)) {
    p <- r_k/(1+r_0)
    llh <- sum((r_k*log(p))[r_k > 0])+ sum(((1 - r_k + r_0) * log(1-p))[(1-r_k+r_0) > 0])
    if (abs(llh - llh.old)/(abs(llh) + 0.1) < epsilon) break    
    llh.old <- llh
    prod_p <- prod(1-p)
    r_0 <- prod_p/(1-prod_p)
  }
  p
}

MCmvbinary_truncated <- function(formula=.~.)
{
    z <- MCmvbinary(formula=formula)
    z@defineComponent <- function(para) {
      predict <- function(x, ...) {
        matrix(para$center, nrow = nrow(x), ncol = length(para$center), 
               byrow = TRUE)
      }
      logLik <- function(x, y) {
        p <- matrix(para$center, nrow = nrow(x), ncol = length(para$center), 
                    byrow = TRUE)
        rowSums(log(y * p + (1 - y) * (1 - p))) - log(1 - prod(1-para$center))
      }
      new("FLXcomponent", parameters = list(center = para$center), df = para$df, 
          logLik = logLik, predict = predict)
    }
    z@fit <- function(x, y, w, ...){
      z@defineComponent(list(center = binary_truncated(y, w), df = ncol(y)))
    }   
    z
}


###**********************************************************

setClass("FLXMCmvcombi",
         representation(binary = "vector"),
         contains = "FLXMC")


FLXMCmvcombi <- function(formula=.~.)
{
    z <- new("FLXMCmvcombi", weighted=TRUE, formula=formula,
             dist = "mvcombi",
             name="model-based binary-Gaussian clustering")

    z@defineComponent <- function(para) {
      predict <- function(x, ...){
        matrix(para$center, nrow=nrow(x), ncol=length(para$center),
               byrow=TRUE)
      }
      
      logLik <- function(x, y){
        if(any(para$binary)){
          p <- matrix(para$center[para$binary], nrow=nrow(x),
                      ncol=sum(para$binary), byrow=TRUE)
          z <- rowSums(log(y[,para$binary,drop=FALSE]*p +
                               (1-y[,para$binary,drop=FALSE])*(1-p)))
        } else z <- rep(0, nrow(x))
        if(!all(para$binary)){
          if(sum(!para$binary)==1)
            z <- z + dnorm(y[,!para$binary],
                           mean=para$center[!para$binary], sd=sqrt(para$var),
                           log=TRUE)
          else
            z <- z + mvtnorm::dmvnorm(y[,!para$binary,drop=FALSE],
                                      mean=para$center[!para$binary], sigma=diag(para$var),
                                      log=TRUE)
        }
        z
      }
            
      new("FLXcomponent", parameters=list(center=para$center, var=para$var), df=para$df,
          logLik=logLik, predict=predict)
    }

    z@fit <- function(x, y, w, binary, ...){
      para <- cov.wt(y, wt=w)[c("center","cov")]
      para <- list(center = para$center, var = diag(para$cov)[!binary],
                   df = ncol(y) + sum(!binary),
                   binary = binary)
      z@defineComponent(para)
    }
    z
}

setMethod("FLXgetModelmatrix", signature(model="FLXMCmvcombi"),
          function(model, data, formula, lhs=TRUE, ...)
{

  model <- callNextMethod(model, data, formula, lhs)
  model@binary <- apply(model@y, 2, function(z) all(unique(z) %in% c(0,1)))
  model
})

setMethod("FLXmstep", signature(model = "FLXMCmvcombi"),
          function(model, weights, components)
{
   return(sapply(seq_len(ncol(weights)),
                 function(k) model@fit(model@x, model@y, weights[,k], model@binary)))
})

