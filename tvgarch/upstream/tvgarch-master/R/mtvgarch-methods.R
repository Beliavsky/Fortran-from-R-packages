#####################################################
## This file contains S3 methods for objects
## of class "mtvgarch":
##
## coef.mtvgarch()
## fitted.mtvgarch()
## logLik.mtvgarch()
## nobs.mtvgarch()
## plot.mtvgarch()
## predict.mtvgarch()
## print.mtvgarch()
## quantile.mtvgarch()
## residuals.mtvgarch()
## summary.mtvgarch()
## toLatex.mtvgarch()
## vcov.mtvgarch() 
##
#####################################################


#####################################################
coef.mtvgarch <- function (object, spec = c("tvgarch", "garch", "tv", "cc"), 
                           ...)
{
  spec <- match.arg(spec)
  if (spec == "tvgarch") {
    if (!is.null(object$order.g) && any(object$order.g != 0)) {
      results <- cbind(object$par.g, object$par.h)
    }
    else {
      results <- object$par.h
    }
  }
  else if (spec == "tv") {
    if (!is.null(object$order.g) && any(object$order.g != 0)) {
      results <- object$par.g
    }
    else {
      stop ("Only GARCH models have been estimated.")
    }
  }
  else if (spec == "garch") {
    if (!is.null(object$order.g) && any(object$order.g != 0)) {
      results <- object$par.h
    }
    else {
      results <- object$par.h
      warning ("Only GARCH models have been estimated.")
    }
  }  
  else if (spec == "cc") {
    if (is.null(object$par.dcc)) {
      results <- object$ccc[1,]
    }
    if (!is.null(object$par.dcc)) {
      results <- NULL
      results$Qbar <- object$ccc[1,]
      results$par.dcc = object$par.dcc
    }
  } 
  return(results)
}

#####################################################
fitted.mtvgarch <- function (object, spec = c("tvgarch", "garch", "tv", "cc"), 
                             as.zoo = TRUE, ...)
{
  spec <- match.arg(spec)
  if (spec == "tvgarch") {
    if (as.zoo == TRUE) {
      object$sigma2 <- zoo(object$sigma2, order.by = index(object$y))
    }
    results <- object$sigma2
  }
  else if (spec == "tv") {
    if (is.null(object$order.g) || any(object$order.g == 0)) {
      stop ("Only GARCH models have been estimated.")
    }
    if (as.zoo == TRUE) {
      object$g <- zoo(object$g, order.by = index(object$y))
    }
    results <- object$g
  }
  else if (spec == "garch") {
    if (is.null(object$order.g) || any(object$order.g == 0)) {
      warning ("Only GARCH models have been estimated.")
    }
    if (as.zoo == TRUE) {
      object$h <- zoo(object$h, order.by = index(object$y))
    }
    results <- object$h
  }  
  else {
    if (!is.null(object$par.dcc)) {
      if (as.zoo == TRUE) {
        object$dcc <- zoo(object$dcc, order.by = index(object$y))
      }
      results = object$dcc
    }
    if (is.null(object$par.dcc)) {
      if (as.zoo == TRUE) {
        object$ccc <- zoo(object$ccc, order.by = index(object$y))
      }
      results = object$ccc
    }
  }
  return(results)
}

#####################################################
logLik.mtvgarch <- function (object, ...)
{
  results <- list()
  m <- ncol(object$y)
  n <- nrow(object$y)
  npar <- 0
  for(i in 1:m){
    object.i <- object$Objs[[paste("obj", i, sep = "")]]
    name <- object$names.y[i]
    if (!is.null(object$order.g) && object$order.g[i,1] != 0) {
      spec <- "tvgarch"
    }
    if (is.null(object$order.g[i,1]) || object$order.g[i,1] == 0){ 
      spec <- "garch"
    }
    results[[paste("logLik", paste(spec), paste(name), sep = "_")]] <- 
      logLik.tvgarch(object = object.i)
    npar <- npar + attr(results[[paste("logLik", paste(spec), paste(name), 
                                       sep = "_")]], "df")
  }
  if (is.null(object$par.dcc)){
    z <- object$residuals
    sigma2 <- object$sigma2
    R <- matrix(1, m, m)
    R[lower.tri(R, FALSE)] <- R[upper.tri(R, FALSE)] <- object$ccc[1,]
    lndetR <- log(det(R))
    invR <- solve(R)
    lf <- -0.5*n*m*log(2*pi) - 0.5*sum(log(sigma2)) - 0.5*n*lndetR - 
      0.5*sum((z %*% invR)*z)
    class(lf) <- "logLik"
    attr(lf, "df") <- npar
    attr(lf, "nobs") <- n
    results[["ccc"]] <- lf
  } 
  if (!is.null(object$par.dcc)) {  
    lf <- object$logLik.dcc
    npar <- npar + 2
    class(lf) <- "logLik"
    attr(lf, "df") <- npar
    attr(lf, "nobs") <- n
    results[["dcc"]] <- lf
  }
  return(results)
}

#####################################################
nobs.mtvgarch <- function (object, ...)
{
  return(nrow(object$y))
}

#####################################################
plot.mtvgarch <- function (x, spec = c("tvgarch", "garch", "tv"), ...) 
{
  for (i in 1:ncol(x$y)) {
    xObj <- x$Objs[[paste("obj", i, sep = "")]]
    xObj$y <- matrix(xObj$y, length(xObj$y), 1)
    colnames(xObj$y) <- paste(colnames(x$y)[i])
    plot.tvgarch(x = xObj, spec = spec)
  }
}

#####################################################
predict.mtvgarch <- function (object, n.ahead = 10, newxtv = NULL, 
                              newxreg = NULL, newindex = NULL, n.sim = 5000, 
                              as.zoo = TRUE, verbose = FALSE, ...)
{
  m <- ncol(object$y)
  n <- nrow(object$y)
  ncol.i <- 0
  if (verbose == FALSE) {
    v <- m
  }
  if (verbose == TRUE) {
    v <- m*(n.sim + 1)
    if (!is.null(object$order.g)) v <- v + colSums(object$order.g != 0)[1]
  }
  if (!is.null(object$xreg)) {
    newxreg <- as.matrix(newxreg)
    if ((NROW(newxreg) != n.ahead)) {
      stop("NROW(newxreg) is unequal to n.ahead")
    }
  }
  results <- matrix(NA, n.ahead, v)
  colnames(results) <- paste("p", 1:v, sep = "_")
  if (object$spillovers == FALSE) {
    for (i in 1:m) {
      object.i <- object$Objs[[paste("obj", i, sep = "")]]
      name <- object$names.y[i]
      if (!is.null(object$xreg)) {
        if (any(object$order.x[i,] != 0)) {
          newxreg.i <- newxreg[,which(object$order.x[i,] == 1)]
        }
        else newxreg.i <- NULL
      }
      result <- predict.tvgarch(object = object.i, n.ahead = n.ahead,
                                newxtv = newxtv, newxreg = newxreg.i,
                                newindex = newindex, n.sim = n.sim,
                                as.zoo = as.zoo, verbose = verbose)
      results[,(ncol.i+1):(ncol.i+NCOL(result))] <- result
      colnames(results)[(ncol.i+1)] <- paste("p_sigma2", name, sep = "_")
      if (verbose == TRUE) {
          if (is.null(object$order.g) || object$order.g[i,1] == 0) {
            colnames(results)[(ncol.i+2):(ncol.i+NCOL(result))] <- 
              c(paste("p_h", 1:n.sim, name, sep = "_"))
          }
        if (!is.null(object$order.g) && object$order.g[i,1] != 0) {
          colnames(results)[(ncol.i+2):(ncol.i+NCOL(result))] <- 
            c(paste("p_g", name, sep = "_"), paste("p_h", 1:n.sim, name, 
                                                   sep = "_"))
        }
      }
      if (i != m) ncol.i <- ncol.i + NCOL(result)
    }
  }
  if (object$spillovers == TRUE) {
    max.p <- max(object$order.h[,1])
    max.q <- max(object$order.h[,2])
    max.r <- max(object$order.h[,3])
    if (any(c(max.q, max.q, max.q) > 1)) {
      stop("Maximum order.h = (1,1,1) with volatility spillovers when predicting
           variances.")
    }
    max.x <- ncol(object$order.x)
    object$par.h[is.na(object$par.h)] <- 0
    intercept.h <- as.vector(object$par.h[,1])
    ARCH <- matrix(0, m, m)
    diag(ARCH) <- object$par.h[,2]
    par.xreg <- object$par.h[,(2+max.q+max.p+max.r):(1+max.q+max.p+max.r+max.x)]
    ARCH[which(object$order.x == 1)] <- par.xreg[which(object$order.x == 1)]
    GARCH <- matrix(0, m, m)
    ASYM <- matrix(0, m, m)
    if (max.p == 1) {
      diag(GARCH) <- object$par.h[,(2+max.q):(1+max.q+max.p)]
    }
    if (max.r == 1) {
      diag(ASYM) <- object$par.h[,(2+max.q+max.p):(1+max.q+max.p+max.r)]
    }
    backcast.values <- list()
    if (max(object$order.h) > 0) {
      backcast.values$innovations <- tail(object$residuals,
                                          n = max(object$order.h))
      backcast.values$z2 <- backcast.values$innovations^2
      backcast.values$Ineg <- backcast.values$innovations < 0
      backcast.values$h <- tail(object$h, n = max(object$order.h))
    }
    etahat <- object$residuals
    if (verbose == TRUE) predictions.h.sim <- NULL
    predictions.h <- matrix(0, n.ahead, m)
    for (i in 1:n.sim) {
      draws <- runif(n.ahead, min = 0.5+1e-09, max = nrow(etahat)+0.5-1e-09)
      draws <- round(draws, digits = 0)
      innovations <- matrix(etahat[draws,], n.ahead, m)
      Ineg <- innovations < 0   
      y2 <- as.vector(backcast.values$z2[1,])
      ht <- as.vector(backcast.values$h[1:max.p,])
      Ineg.y2 <- as.vector((backcast.values$Ineg*y2)[1:max.r,])    
      h <- matrix(NA_real_, n.ahead, m)
      for (t in 1:n.ahead) {
        ht <- intercept.h + ARCH%*%y2 + GARCH%*%ht + ASYM%*%(Ineg.y2)
        h[t,] <- ht
        y2 <- h[t,]*innovations[t,]^2
        Ineg.y2 <- Ineg[t,]*y2
      }
      predictions.h <- predictions.h + h
      if (verbose == TRUE) {
        colnames(h) <- paste(paste("p_h", paste(object$names.y), sep = "_"), i, 
                             sep = "_")
        predictions.h.sim <- cbind(predictions.h.sim, h)
      }
    }
    predictions.h <- predictions.h/n.sim
    if (!is.null(object$order.g) && any(object$order.g != 0)) {
      coefs.g <- object$par.g
      if (!is.null(newxtv)) {
        if ((NROW(newxtv) != n.ahead)) {
          stop("NROW(newxtv) is unequal to n.ahead")
        }
        predictions.g <- matrix(NA_real_, n.ahead, m)
        newxtv <- as.matrix(newxtv)
        for (i in 1:m) {
          predictions.g[,i] <- tvObj(par.g = coefs.g[i,], fixed.par.g = NULL, 
                                 xtv = newxtv, opt = object$opt, 
                                 order.g = object$order.g[i,], fixed.h = NULL, 
                                 y = NULL, iter0 = TRUE, flag = 2)
        }
      }
      if (is.null(newxtv)) {
        predictions.g <- matrix(coefs.g[,1], n.ahead, m, byrow = TRUE)
      }
      colnames(predictions.g) <- paste("p_g", object$names.y, sep = "_")
      results <- predictions.h*predictions.g
      colnames(results) <- paste("p_sigma2", object$names.y, sep = "_")
      if (verbose == TRUE) {
        results <- cbind(results, predictions.g, predictions.h.sim)
      }
    }
    if (is.null(object$order.g) || any(object$order.g == 0)) {
      results <- predictions.h
      colnames(results) <- c(paste("p_sigma2", object$names.y, sep = "_"))
      if (verbose == TRUE) {
        results <- cbind(results, predictions.h.sim)
      }
    }
  }
  if (is.null(newindex)) {
    newindex <- 1:n.ahead
  }
  if (as.zoo == TRUE) {
    results <- zoo(results, order.by = newindex)
  }
  return(results)
}

##################################################### 
print.mtvgarch <- function (x, ...)
{
  cat("\n")
  cat("**************************** \n")
  cat("* The Volatility Component * \n")
  cat("**************************** \n")
  for (i in 1:ncol(x$y)) {
    cat("\n")
    cat("*** Equation ", i, " (",colnames(x$y)[i], ") ***\n", sep = "")
    print.tvgarch(x = x$Objs[[paste("obj", i, sep = "")]])
  }
  cat("\n")
  cat("***************************** \n")
  cat("* The Correlation Component * \n")
  cat("***************************** \n")
  if (is.null(x$par.dcc)) {
    cat("\n")
    cat("Constant Conditional Correlations: \n")
    cat(" \n")
    z <- x$residuals
    nobs <- nrow(x$y)
    m <- ncol(x$y)
    sigma2 <- x$sigma2
    lndetR <- log(det(x$R))
    invR <- solve(x$R)
    lf <- -0.5*nobs*m*log(2*pi) - 0.5*sum(log(sigma2)) - 0.5*nobs*lndetR - 
      0.5*sum((z %*% invR)*z)
    print(round(x$R, digits = 4))
    cat(" \n")
    loglik <- as.matrix(lf)
    rownames(loglik) <- "Log-likelihood:"
    colnames(loglik) <- ""
    print(round(loglik, digits = 4))
  }
  if (!is.null(x$par.dcc)) {
    cat("\n")
    cat("Dynamic Conditional Correlations \n")
    cat(" \n")
    cat(" \n")
    cat("Q matrix \n")
    cat(" \n")
    print(round(x$R, digits = 4))
    cat("\n")
    cat("Coefficients \n")
    cat(" \n")
    if (x$turbo == TRUE) {
      x$se.dcc <- 
        matrix(sqrt(diag(vcov.mtvgarch(object = x, spec = "cc"))), 1, 2)
    }
    estimates.dcc <- as.matrix(rbind(x$par.dcc, x$se.dcc))
    rownames(estimates.dcc) <- c("Estimate:", "Std. Error:")
    colnames(estimates.dcc) <- c("alpha", "beta")
    print(round(estimates.dcc,4))
    loglik <- as.matrix(x$logLik.dcc)
    rownames(loglik) <- "Log-likelihood:"
    colnames(loglik) <- ""
    print(round(loglik, digits = 4))
    cat("\n")
  }
}
  
#####################################################
quantile.mtvgarch <- function (x, probs = 0.025, type = 7, as.zoo = TRUE, ...)
{
  m <- ncol(x$y)
  n <- nrow(x$y)
  ncol.i <- 0
  iCols <- length(probs)
  if (iCols == 1) {
    results <- matrix(NA, n, m)
  }
  if (iCols > 1) {
    results <- matrix(NA, n, m*iCols)
    colnames(results) <- paste("q", 1:ncol(results), sep = "_")
  }  
  for(i in 1:m){
    result.i <- quantile.tvgarch(x = x$Objs[[paste("obj", i, sep = "")]], 
                                 probs = probs, names = TRUE, type = type, 
                                 as.zoo = as.zoo)
    name <- x$names.y[i]
    results[,(ncol.i+1):(ncol.i+iCols)] <- result.i
    if (iCols > 1) {
      colnames(results)[(ncol.i+1):(ncol.i+iCols)] <-  
        paste(colnames(result.i), name, sep = "_")
    }
    if (i != m) ncol.i <- ncol.i + iCols
  }
  if (iCols == 1) {
    colnames(results) <- paste(paste(probs*100, "%", sep = ""), 
                               x$names.y, sep = "_")
  }
  return(results)
}

#####################################################
residuals.mtvgarch <- function (object, as.zoo = TRUE, ...)
{
  if (as.zoo == TRUE) {
    object$residuals <- zoo(object$residuals, order.by = index(object$y))
  }
  return(object$residuals)
}

#####################################################
summary.mtvgarch <- function (object, ...)
{
  for (i in 1:ncol(object$y)) {
    cat("\n")
    cat(colnames(object$y)[i], ": \n", sep = "")
    summary.tvgarch(object = object$Objs[[paste("obj", i, sep = "")]], ...)
  }
  if (is.null(object$par.dcc)) {
    cat("\n")
    cat("Constant Conditional Correlations: \n")
    cat(" \n")
    z <- object$residuals
    nobs <- nrow(object$y)
    m <- ncol(object$y)
    sigma2 <- object$sigma2
    lndetR <- log(det(object$R))
    invR <- solve(object$R)
    lf <- -0.5*nobs*m*log(2*pi) -0.5*sum(log(sigma2)) -0.5*nobs*lndetR -
      0.5*sum((z %*% invR)*z)
    print(round(object$R, digits = 4))
    cat(" \n")
    loglik <- as.matrix(lf)
    rownames(loglik) <- "Log-likelihood:"
    colnames(loglik) <- ""
    print(round(loglik, digits = 4))
  }
  if (!is.null(object$par.dcc)) {
    cat("\n")
    cat("Dynamic Conditional Correlations \n")
    cat(" \n")
    cat(" \n")
    cat("Q matrix \n")
    cat(" \n")
    print(round(object$R, digits = 4))
    cat("\n")
    cat("Coefficients \n")
    cat(" \n")
    if (object$turbo == TRUE) {
      object$se.dcc <- 
        matrix(sqrt(diag(vcov.mtvgarch(object = object, 
                                       spec = "cc"))), 1, 2)
    }
    estimates.dcc <- as.matrix(rbind(object$par.dcc, object$se.dcc))
    rownames(estimates.dcc) <- c("Estimate:", "Std. Error:")
    colnames(estimates.dcc) <- c("alpha", "beta")
    print(round(estimates.dcc,4))
    loglik <- as.matrix(object$logLik.dcc)
    rownames(loglik) <- "Log-likelihood:"
    colnames(loglik) <- ""
    print(round(loglik, digits = 4))
    cat("\n")
  }
}

#####################################################
toLatex.mtvgarch <- function (object, digits = 4, ...)
{
  m <- ncol(object$y)
  names.y <- colnames(object$y)
  names.sigma2 <- colnames(object$sigma2)
  if(!is.null(object$xreg) && object$spillovers == FALSE) {
    names.x <- colnames(object$xreg)
  }
  if (!is.null(object$xreg) && object$spillovers == TRUE) names.x <- names.y
  cat("Variance equations \n\n") 
  cat("\\begin{eqnarray} \n")
  for (i in 1:m) {
    object.i <- object$Objs[[paste("obj", i, sep = "")]]
    object.i$date <- object$date
    object.i$turbo <- object$turbo
    object.i$trace <- object$trace
    coefs.h <- as.numeric(object.i$par.h)
    if (object.i$order.h[1] != 0) {
      if (is.null(object$order.g) || object$order.g[i,1] == 0) {
        names.h <- c("", paste("y^{2}_{",paste(i, sep = ""),",t-",
                               paste(seq(1:object.i$order.h[2]), sep = ""),"}", sep = ""))
      }
    }
    if (object.i$order.h[1] != 0 && !is.null(object$order.g) && 
        object$order.g[i,1] != 0) {
      names.h <- c("", paste("\\dfrac{y^{2}_{",paste(i, sep = ""),",t-",
                             paste(seq(1:object.i$order.h[2]), sep = ""),"}}{\\widehat{g}_{",
                             paste(i, sep = ""),",t-", paste(seq(1:object.i$order.h[2]), 
                                                             sep = ""),"}}", sep = ""))
    }
    if (object.i$order.h[1] != 0 && !is.null(object$order.g) && 
        object$order.g[i,1] != 0) {
      names.h <- c(names.h, paste("\\widehat{h}_{",
                                  paste(i, sep = ""),",t-",paste(seq(1:object.i$order.h[1]), 
                                                                 sep = ""),"}", sep = "")) 
    }
    if (object.i$order.h[1] != 0) {
      if (is.null(object$order.g) || object$order.g[i,1] == 0) {
        names.h <- c(names.h, paste("\\widehat{\\sigma}^2_{",
                                    paste(i, sep = ""),",t-",paste(seq(1:object.i$order.h[1]), 
                                                                   sep = ""),"}", sep = "")) 
      }
    }
    if (object.i$order.h[3] != 0) {
      names.h <- c(names.h, paste("\\dfrac{y^{2}_{",paste(i, sep = ""),",t-",
                                  paste(seq(1:object.i$order.h[3]), sep = ""),"}}{\\widehat{g}_{", 
                                  paste(i, sep = ""),",t-", 
                                  paste(seq(1:object.i$order.h[3]), sep = ""),"}}\\text{I}(",
                                  paste("y_{",paste(i, sep = ""),",t-",
                                        paste(seq(1:object.i$order.h[3]), sep = ""),"}", 
                                        sep = "")," < 0)", sep = "")) 
    }
    if (!is.null(object$xreg) && object$spillovers == FALSE) {
      names.h <- c(names.h, paste("x_{",paste(1:ncol(object$xreg),
                                              sep = ""),"t}", sep = ""))
    }
    if (!is.null(object$xreg) && object$spillovers == TRUE){
      for (x in which(object$order.x[i,] == 1)) {
        if (object$order.g[x,1] == 0) {
          names.h <- c(names.h, paste("y^{2}_{",paste(x,sep = ""),",t-1}", 
                                      sep = ""))
        }
        if (object$order.g[x,1] == 1) {
          names.h <- c(names.h, 
                       paste("\\dfrac{y^{2}_{", paste(x, sep = ""),
                             ",t-1}}{\\widehat{g}_{", paste(x, sep = ""),
                             ",t-1}}", sep = ""))
        }
      }
    }  
    coefsNames.h <- names.h
    coefs.h <- as.numeric(coefs.h)
    if (object.i$turbo == TRUE) {
      if (is.null(object$order.g) || object$order.g[i,1] == 0) {
        object.i$se.h <- sqrt(diag(vcov.tvgarch(object = object.i)))
      }
      if (!is.null(object$order.g) && object$order.g[i,1] != 0) {
        object.i$se.h <- sqrt(diag(vcov.tvgarch(object = object.i, 
                                                spec = "garch")))
        object.i$se.g <- sqrt(diag(vcov.tvgarch(object = object.i, 
                                                spec = "tv")))
        s <- length(object.i$order.g)
      }
    }
    stderrs <- as.numeric(object.i$se.h)
    eqtxt.h <- NULL
    for (j in 1:length(coefs.h)) {
      ifpluss <- ifelse(j==1, "", " + ")
      eqtxt.h <- paste(eqtxt.h, ifelse(coefs.h[j]<0, " - ",ifpluss), 
                       "\\underset{(", format(round(stderrs[j], digits=digits), 
                                              nsmall=digits),")}{",
                       format(round(abs(coefs.h[j]), digits=digits), 
                              nsmall=digits),"}",coefsNames.h[j], sep="")
    }
    txtAddEq1 <- " \\\\[1mm]"
    txtAddEq2 <- " \\nonumber \\\\[1mm]"
    if (is.null(object$order.g) || object$order.g[i,1] == 0) {
      eqtxt.h <- paste0("  \\sigma^2_{",paste(i,sep = ""),"t} &=& ", 
                        eqtxt.h, "", txtAddEq2, " \n")
    }
    if (!is.null(object$order.g) && object$order.g[i,1] != 0) {
      eqtxt.h <- paste0("  \\widehat{h}_{",paste(i,sep = ""),"t} &=& ", eqtxt.h, 
                        "", txtAddEq2, " \n")
      if (!is.null(object.i$order.g)) {
        s <- length(object.i$order.g)
        coefs.g <- object.i$par.g[1:(s+1)]
        coefs.tv <- object.i$par.g[-(1:(s+1))]
        coefsNames.g <- ""
        for(j in 1:s){
          coefsNames.g <- c(coefsNames.g, paste("\\widehat{G}_{",
                                                paste(i, sep = ""),
                                                paste(j, sep = ""),"}", 
                                                sep = ""))
        }
        coefs.g <- as.numeric(coefs.g)
        stderrs.g <- as.numeric(object.i$se.g[1:s])
        stderrs.tv <- as.numeric(object.i$se.g[-(1:s)])
        eqtxt.g <- NULL
        for (j in 1:length(coefs.g)){
          ifpluss <- ifelse(j == 1, "", " + ")
          eqtxt.g <- paste(eqtxt.g, ifelse(coefs.g[j] < 0, " - ", ifpluss))
          if (j == 1) eqtxt.g <- paste(eqtxt.g, "\\underset{(-)}{", 
                                       format(round(abs(coefs.g[j]), 
                                                    digits = digits), 
                                              nsmall = digits),"}", sep = "")
          if (j > 1) {
            eqtxt.g <- paste(eqtxt.g, "\\underset{(", 
                             format(round(stderrs.g[j-1], digits = digits), 
                                    nsmall = digits),")}{", 
                             format(round(abs(coefs.g[j]), digits = digits), 
                                    nsmall = digits),"}","\\widehat{G}_{", 
                             paste(j-1, sep = ""),"}(", "\\underset{(", 
                             format(round(stderrs.tv[j-1], digits = digits), 
                                    nsmall = digits),")}{", 
                             format(round(abs(coefs.tv[j-1]), digits = digits), 
                                    nsmall = digits),"};", sep = "")
            for(k in 1:object.i$order.g[j-1]){
              eqtxt.g <- paste(eqtxt.g, "\\underset{(", 
                               format(round(stderrs.tv[s+k], digits = digits), 
                                      nsmall = digits),")}{", 
                               format(round(abs(coefs.tv[s+k]), 
                                            digits = digits), nsmall = digits), 
                               "}", sep = "")
              if(k < object.i$order.g[j-1]) eqtxt.g <- paste(eqtxt.g, ",", 
                                                             sep = "")
            }
            if (colnames(object.i$xtv) == "time") {
              eqtxt.g <- paste(eqtxt.g, "; t/n)", sep = "")
            }
            else eqtxt.g <- paste(eqtxt.g, "; s_{t})", sep = "")
          }
        }
        txtAddEq <- " \\\\[1mm]"
        eqtxt.g <- paste0("  \\widehat{g}_{", paste(i, sep = ""),"t} &=& ", 
                          eqtxt.g, "", txtAddEq, " \n")
      }
      else eqtxt.g <- paste0("  \\widehat{g}_{", paste(i, sep = ""),
                             "t} &=& 1 \\\\[1mm] \n")
    }
    goftxt <- NULL
    goftxt <- "&&"
    iT <- length(object.i$sigma2)
    goftxt <- paste(goftxt, " \\text{Log-likelihood: }", 
                    format(round(as.numeric(logLik.tvgarch(object = object.i)), 
                                 digits=digits), nsmall=digits), "\\qquad n = ", 
                    iT, " \\nonumber \n", sep = "")
    if( i != m) goftxt <- paste(goftxt, " \\\\ [1mm] \n", sep = "")
    cat(eqtxt.h)
    if (!is.null(object$order.g) && object$order.g[i,1] != 0) cat(eqtxt.g)
    cat(goftxt)
  }
  cat("\\end{eqnarray} \n")
  cat(" \n")
  if (is.null(object$par.dcc)) {
    cat("Constant conditional correlations: \n")
    cat(" \n")
    cat("\\begin{table}[!htbp] \\centering \n")
    cat("\\caption{Constant conditional correlations.} \n") 
  }
  if (!is.null(object$par.dcc)) {
    cat("Dynamic conditional correlations: \n")
    cat(" \n")
    cat("\\begin{table}[!htbp] \\centering \n")
    cat("\\caption{Dynamic conditional correlations (unconditional matrix)} \n") 
  }
  cat("\\label{} \n") 
  cat("\\begin{tabular}{@{\\extracolsep{5pt}}", paste(rep("c", m + 1), 
                                                      sep = ""), "} \n")
  cat("\\\\[-1.8ex] \\hline \n")
  cat("\\hline \\\\[-1.8ex] \n") 
  namesR <- paste(" & ") 
  mtxR <- NULL 
  for (i in 1:m) {
    if (i < m) namesR <- paste(namesR, colnames(object$R)[i], " & ")
    if (i == m) namesR <- paste(namesR, colnames(object$R)[i], " \\\\ \n ")
    for (j in 1:m) {
      if (j == 1) {
        mtxR <- paste(mtxR, colnames(round(object$R, digits))[i], " & ")
      }
      if (j < m) {
        mtxR <- paste(mtxR, round(object$R, digits)[i,j], " & ")
      }
      if (j == m) {
        mtxR <- paste(mtxR, round(object$R, digits)[i,j], " \\\\ \n ")
      }
    }
  }
  cat(namesR)
  cat(mtxR)
  cat("\\hline \\\\[-1.8ex] \n") 
  cat("\\end{tabular} \n") 
  cat("\\end{table} \n")   
  if (!is.null(object$par.dcc)) {
    cat("Dynamic conditional correlations: \n")
    coefs <- as.numeric(object$par.dcc)
    stderrs <- as.numeric(object$se.dcc)
    cat("\\begin{table}[!htbp] \\centering \n")
    cat("\\caption{} \n") 
    cat("\\label{} \n") 
    cat("\\begin{tabular}{@{\\extracolsep{5pt}} cc} \n")
    cat("\\\\[-1.8ex] \\hline \n")
    cat("\\hline \\\\[-1.8ex] \n") 
    cat(" $ \\alpha $ & $ \\beta $ \\\\ \n")
    cat("\\hline \\\\[-1.8ex] \n")
    cat(paste(round(coefs[1], digits), sep = "")," & ", 
        paste(round(coefs[2], digits), sep = ""),"\\\\ \n")
    cat("(",paste(round(stderrs[1], digits), sep = ""),") & (", 
        paste(round(stderrs[2], digits), sep = ""),")\\\\ \n")
    cat("\\hline \\\\[-1.8ex] \n") 
    cat("\\end{tabular} \n") 
    cat("\\end{table} \n") 
  }
}

#####################################################
vcov.mtvgarch <- function (object, spec = c("tvgarch", "garch", "tv", "cc"), 
                           ...)
{
  spec <- match.arg(spec)
  if (spec != "cc") {
    results <- list()
    m <- ncol(object$y)
    for (i in 1:m) {
      name <- object$names.y[i]
      object.i <- object$Objs[[paste("obj", i, sep = "")]]
      if (is.null(object$order.g) || any(object$order.g == 0)) {
        if (spec == "tv" || spec == "garch") {
          stop ("Only GARCH models have been estimated.")
        }
      }
      if (is.null(object$order.g) || object$order.g[i,1] == 0) {
          spec.i <- "sigma2"
      }
      if (!is.null(object$order.g) && object$order.g[i,1] != 0) {
        spec.i <- spec
      }
      results[[paste("vcov", paste(spec.i), name, sep = "_")]] <- 
        vcov.tvgarch(object = object.i, spec = spec.i)
    }    
  }
  if (spec == "cc") {  
    if (object$turbo == TRUE) {
      jac.dcc <- jacobian(func = dccObj, x = object$par.dcc, 
                          z = object$residuals, sigma2 = object$sigma2, 
                          flag = 0)
      J.dcc <- crossprod(jac.dcc)  
      H.dcc <- optimHess(par = object$par.dcc, fn = dccObj, 
                         z = object$residuals, sigma2 = object$sigma2, flag = 1)
      solHdcc <- solve(-H.dcc)
      object$vcov.dcc <- solHdcc %*% J.dcc %*% solHdcc
    }
    results <- object$vcov.dcc
  }
  return(results)
}
