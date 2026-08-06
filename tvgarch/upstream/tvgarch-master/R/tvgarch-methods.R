#####################################################
## This file contains S3 methods for objects
## of class "tvgarch":
##
## coef.tvgarch()
## fitted.tvgarch()
## logLik.tvgarch()
## nobs.tvgarch()
## plot.tvgarch()
## predict.tvgarch()
## print.tvgarch()
## quantile.tvgarch()
## residuals.tvgarch()
## summary.tvgarch()
## toLatex.tvgarch()
## vcov.tvgarch() 
##
#####################################################


#####################################################
coef.tvgarch <- function (object, spec = c("tvgarch", "garch", "tv"), ...)
{
  if (!is.null(object$order.g) && object$order.g[1] != 0) {
    spec <- match.arg(spec)
    if (spec == "tvgarch") {
      results <- c(object$par.g,object$par.h)
    }
    else if (spec == "tv") results <- object$par.g
    else if (spec == "garch") results <- object$par.h
  } 
  if (is.null(object$order.g) || object$order.g[1] == 0) results <- object$par.h
  return(results)
}

#####################################################
fitted.tvgarch <- function (object, spec = c("tvgarch", "garch", "tv"), 
                            as.zoo = TRUE, ...)
{
  if (!is.null(object$order.g) && object$order.g[1] != 0) {
    spec <- match.arg(spec)
    if (spec == "tv") results <- object$g 
    else if (spec == "garch") results <- object$h 
    else if (spec == "tvgarch") results <- object$sigma2 
  }
  if (is.null(object$order.g) || object$order.g[1] == 0) {
    results <- object$sigma2
  }
  if (as.zoo == TRUE) {
    results <- zoo(results, order.by = index(object$y))
  }
  return(results)
}

#####################################################
logLik.tvgarch <- function (object, ...)
{
    result <- object$logLik
    class(result) <- "logLik"
    attr(result, "df") <- length(coef.tvgarch(object = object))
    attr(result, "nobs") <- length(object$sigma2)
    return(result)
}

#####################################################
nobs.tvgarch <- function (object, ...)
{
  return(length(object$sigma2))
}

#####################################################
plot.tvgarch <- function (x, spec = c("tvgarch", "garch", "tv"), ...) 
{
  if (!is.null(x$order.g) && x$order.g[1] != 0) {
    spec <- match.arg(spec)
    if (spec == "tvgarch") {
      par(mfrow = c(3,1))
    }
    if (spec == "tvgarch" || spec == "tv") {
      garchEst <- tvgarch(x$y, order.g = NULL, order.h = x$order.h, 
                          xreg = x$xreg)
      plot(x$y.index, sqrt(fitted(garchEst)), 
           type = "l", ylab = "", xlab = "", ...)
      lines(x$y.index, sqrt(x$g), col = "blue", ...)
      if (spec == "tvgarch") {
        title(paste("Series: ", paste(colnames(x$y)), "\n\n GARCH"))
      }
      else title("GARCH (black) and TV in TV-GARCH (blue)")
    }
    if (spec == "tvgarch") {
      plot(x$y.index, sqrt(x$sigma2), type = "l", 
           ylab = "", xlab = "", ...)
      title("TV-GARCH")
    }
    if (spec == "tvgarch" || spec == "garch") {
      plot(x$y.index, sqrt(x$h), type = "l", ylab = "",
           xlab = "", ...)
      title("GARCH in TV-GARCH")
    }
    par(mfrow=c(1,1))
  }
  if (is.null(x$order.g) || x$order.g[1] == 0) {
    plot(x$y.index, sqrt(x$h), type = "l", ylab = "", 
         xlab = "", ...)
    title(paste("Series: ", paste(colnames(x$y)), "\n\n GARCH"))
  }
}

#####################################################
predict.tvgarch <- function (object, n.ahead = 10, newxtv = NULL, 
                             newxreg = NULL, newindex = NULL, n.sim = 5000, 
                             as.zoo = TRUE, verbose = FALSE, ...)
{
  n <- length(object$y)
  if (!is.null(object$order.g) && object$order.g[1] != 0) {
      coefs.h <- as.numeric(coef.tvgarch(object = object, spec = "garch"))
    }
  if (is.null(object$order.g) || object$order.g[1] == 0) {
      coefs.h <- as.numeric(coef.tvgarch(object = object, spec = "tvgarch"))
    }
  interceptCoef <- coefs.h[1]
  archCoef <- coefs.h[(1+1):(1+object$order.h[2])]
  if (object$order.h[1] != 0) {
      garchCoef <- coefs.h[(1+object$order.h[2]+1):(1+sum(object$order.h[1:2]))]
    }
  if (object$order.h[1] == 0) garchCoef <- NULL
  if (object$order.h[3] != 0) {
      asymCoef <- 
        coefs.h[(1+sum(object$order.h[1:2])+1):(1+sum(object$order.h))]
    }
  if (object$order.h[3] == 0) asymCoef <- NULL
  if (!is.null(object$xreg)) xregCoef <- coefs.h[-(1:(1+sum(object$order.h)))]
  if (is.null(object$xreg)) xregCoef <- NULL
  backcast.values <- list()
  if (max(object$order.h) > 0) {
    backcast.values$innovations <- tail(as.numeric(object$residuals), 
                                        n = max(object$order.h))
    backcast.values$z2 <- backcast.values$innovations^2
    backcast.values$Ineg <- as.numeric(backcast.values$innovations < 0)
    backcast.values$sigma2 <- tail(as.numeric(object$h), 
                                   n = max(object$order.h))
    if (!is.null(object$xreg)) {
      backcast.values$xreg <- tail(as.numeric(object$xreg%*%xregCoef), 
                                   n = max(object$order.h))
    }
  }
  xreg <- NULL
  if (!is.null(object$xreg)) {
      if ((NROW(newxreg) != n.ahead)) {
        stop("NROW(newxreg) is unequal to n.ahead")
      }
      xreg <- cbind(newxreg) %*% xregCoef
    }
  etahat <- coredata(object$residuals)
  draws <- runif(n.ahead*n.sim, min = 0.5+1e-09, 
                 max = length(etahat)+0.5-1e-09)
  draws <- round(draws, digits = 0)
  innovations <- matrix(etahat[draws], n.ahead, n.sim)
  predictions.h <- matrix(NA, n.ahead, n.sim)
  for (i in 1:n.sim) {
      predictions.h[,i] <- garchxSim(n = n.ahead, intercept = interceptCoef, 
                                   arch = archCoef, garch = garchCoef, 
                                   asym = asymCoef, xreg = xreg, 
                                   innovations = innovations[,i], 
                                   backcast.values = backcast.values, 
                                   verbose = TRUE, as.zoo = FALSE)[,"sigma2"]
    }
  if (!is.null(object$order.g) && object$order.g[1] != 0) {
      coefs.g <- as.numeric(coef.tvgarch(object = object, spec = "tv"))
      if (!is.null(newxtv)) {
        if ((NROW(newxtv) != n.ahead)) {
          stop("NROW(newxtv) is unequal to n.ahead")
        }
        newxtv <- as.matrix(newxtv)
        predictions.g <- tvObj(par.g = coefs.g, fixed.par.g = NULL, 
                               xtv = newxtv, opt = object$opt, 
                               order.g = object$order.g, fixed.h = NULL, 
                               y = NULL, iter0 = TRUE, flag = 2)
      }
      if (is.null(newxtv)) predictions.g <- rep(coefs.g[1], n.ahead)
      result <- as.matrix(rowMeans(predictions.h)*predictions.g)
      if (verbose == TRUE) {
        result <- cbind(result, as.vector(predictions.g), predictions.h)
        colnames(result) <- c("p_sigma2", "p_g", paste("p_h", 
                                                       1:NCOL(predictions.h),
                                                       sep = "_"))
      }
    }
  if (is.null(object$order.g) || object$order.g[1] == 0) {
      result <- as.matrix(rowMeans(predictions.h))
      if (verbose == TRUE) {
        result <- cbind(result, predictions.h)
        colnames(result) <- c("p_sigma2", paste("p_h", 1:NCOL(predictions.h), 
                                              sep = "_"))
      }
    }
  if (is.null(newindex)) {
      newindex <- 1:n.ahead
    }
  if (as.zoo == TRUE) {
      result <- zoo(result, order.by = newindex)
    }
  return(result)
}

#####################################################
print.tvgarch <- function (x, ...)
{
  if (!is.null(x$order.g) && x$order.g[1] != 0) {
    loglik <- as.matrix(x$logLik)
    rownames(loglik) <- "Log-likelihood:"
    colnames(loglik) <- ""
    cat("\n")
    cat("Date:", x$date, "\n")
    model <- NULL
    if (length(x$order.g) == 1) model <- paste(model, "TV(", x$order.g, ")-", 
                                               sep = "") 
    else {
      for (i in 1:length(x$order.g)){
        if (i == 1) model <- paste(model, "TV(", x$order.g[i], sep = "") 
        else if (i == length(x$order.g)) {
          model <- paste(model, ",", x$order.g[i], ")-", sep = "") 
        }
        else model <- paste(model, ",", x$order.g[i], sep = "")
      }
    }
    if (x$order.h[3] != 0) {
      model <- paste(model, "GJR(", x$order.h[3], ")-", sep = "") 
    }
    if (x$order.h[1] != 0) {
      model <- paste(model, "G", sep = "") 
    }
    model <- paste(model, "ARCH(", sep = "") 
    if (x$order.h[1] != 0) {
      model <- paste(model, x$order.h[1], ",", sep = "") 
    }
    model <- paste(model, x$order.h[2], ")", sep = "") 
    if (!is.null(x$xreg)) model <- paste(model, "-X", sep = "") 
    cat("Model:", model, "\n")
    cat("Method: maximization by parts \n")
    cat("No. of iterations:", x$iter, "\n")
    cat("No. of observations:", length(x$sigma2), "\n")
    cat("\n")  
    cat("* TV specification (long-term component) * \n")
    cat("\n")
    cat("Optimization: linearly-constrained \n")
    if (!is.null(x$message.g)) cat("Message (constrOptim):", x$message.g, "\n")
    cat("Transition variable:", colnames(x$xtv), "\n")
    cat("\n")
    cat("intercept:", x$par.g[1], "(fixed from first iteration) \n")
    cat("\n")
    if (x$turbo == TRUE) {
      x$se.g <- sqrt(diag(vcov.tvgarch(object = x, spec = "tv")))
    }
    else x$se.g <- x$se.g[-1]
    estimates.g <- as.matrix(rbind(x$par.g[-1], x$se.g))
    rownames(estimates.g) <- c("Estimate:", "Std. Error:")
    colnames(estimates.g) <- x$names.g[-1] 
    print(round(estimates.g, 7))
    cat("\n")
    cat("* GARCH specification (short-term component) * \n")
    cat("\n")
    cat("Optimization: box-constrained \n")
    if (!is.null(x$message.h)) cat("Message (nlminb):", x$message.h, "\n")
    cat("\n")
    if (x$turbo == TRUE) {
      x$se.h <- sqrt(diag(vcov.tvgarch(object = x, spec = "garch")))
    }
    estimates.h <- as.matrix(rbind(x$par.h, x$se.h))
    rownames(estimates.h) <- c("Estimate:", "Std. Error:")
    colnames(estimates.h) <- x$names.h
    colnames(estimates.h)[1] <- "intercept"
    print(round(estimates.h, 7))
    print(round(loglik, digits = 4))
    cat("\n")
  }
  if (is.null(x$order.g) || x$order.g[1] == 0) {
    pars <- coef.tvgarch(object = x)
    vcovmat <- vcov.tvgarch(object = x)
    out1 <- rbind(pars, sqrt(diag(vcovmat)))
    rownames(out1) <- c("Estimate:", "Std. Error:")
    out2 <- as.data.frame(matrix(NA_real_, 1, 1))
    out2[1, 1] <- as.character(round(logLik.tvgarch(object = x), digits = 4))
    rownames(out2) <- "Log-likelihood:"
    colnames(out2) <- ""
    cat("\n")
    cat("\n")
    cat("Date:", x$date, "\n")
    model <- NULL
    if (x$order.h[3] != 0) {
      model <- paste(model, "GJR(", x$order.h[3], ")-", sep = "") 
    }
    if (x$order.h[1] != 0) {
      model <- paste(model, "G", sep = "") 
    }
    model <- paste(model, "ARCH(", sep = "") 
    if (x$order.h[1] != 0) {
      model <- paste(model, x$order.h[1], ",", sep = "") 
    }
    model <- paste(model, x$order.h[2], ")", sep = "") 
    if (!is.null(x$xreg)) model <- paste(model, "-X", sep = "") 
    cat("Model:", model, "\n")
    cat("Method: normal ML\n")
    cat("Message (nlminb):", x$message.h, "\n")
    cat("No. of observations:", length(x$sigma2),
        "\n")
    cat("\n")  
    print(out1)
    print(out2)
    cat("\n")  
  }
}


#####################################################
quantile.tvgarch <- function (x, probs = 0.025, names = TRUE, type = 7, 
                              as.zoo = TRUE, ...)
{
  etahat <- residuals.tvgarch(object = x)
  sigma <- sqrt(fitted.tvgarch(object = x))
  qvals <- quantile(etahat, probs = probs, names = names, type = type)
  iN <- NROW(etahat)
  iCols <- length(qvals)
  result <- matrix(NA, iN, iCols)
  colnames(result) <- names(qvals)
  for (i in 1:iCols) {
    result[,i] <- sigma*qvals[i]
  }
  if (iCols == 1) {
    result <- as.vector(result)
  }
  if (as.zoo == TRUE) {
    result <- zoo(result, order.by = index(etahat))
  }
  return(result)
}

#####################################################
residuals.tvgarch <- function (object, as.zoo = TRUE, ...)
{
  if (as.zoo == TRUE) {
    object$residuals <- zoo(object$residuals, order.by = index(object$y))
  }
  return(object$residuals)
}

#####################################################
summary.tvgarch <- function (object, ...)
{
  loglik <- as.matrix(object$logLik)
  rownames(loglik) <- "Log-likelihood:"
  colnames(loglik) <- ""
  if (!is.null(object$order.g) && object$order.g[1] != 0) {
    cat("\n")
    model <- NULL
    if (length(object$order.g) == 1) {
      model <- paste(model, "TV(", object$order.g, ")-", sep = "") 
    }
    else {
      for (i in 1:length(object$order.g)){
        if (i == 1) model <- paste(model, "TV(", object$order.g[i], sep = "") 
        else if (i == length(object$order.g)) {
          model <- paste(model, ",", object$order.g[i], ")-", sep = "") 
        }
        else model <- paste(model, ",", object$order.g[i], sep = "")
      }
    }
    if (object$order.h[3] != 0) {
      model <- paste(model, "GJR(", object$order.h[3], ")-", sep = "") 
    }
    if (object$order.h[1] != 0) model <- paste(model, "G", sep = "") 
    model <- paste(model, "ARCH(", sep = "") 
    if (object$order.h[1] != 0) model <- paste(model, object$order.h[1], ",", 
                                               sep = "") 
    model <- paste(model, object$order.h[2], ")", sep = "") 
    if (!is.null(object$xreg)) model <- paste(model, "-X", sep = "") 
    cat("Model:", model, "\n")
    cat("Transition variable:", colnames(object$xtv), "\n")
    cat("intercept.g:", object$par.g[1], "(fixed) \n")
    cat("\n")
    if (object$turbo == TRUE) {
      object$se.g <- sqrt(diag(vcov.tvgarch(object = object, spec = "tv")))
    }
    else object$se.g <- object$se.g[-1]
    estimates.g <- as.matrix(rbind(object$par.g[-1], object$se.g))
    rownames(estimates.g) <- c("Estimate:", "Std. Error:")
    colnames(estimates.g) <- object$names.g[-1] 
    if (object$turbo == TRUE) {
      object$se.h <- sqrt(diag(vcov.tvgarch(object = object, spec = "garch")))
    }
    estimates.h <- as.matrix(rbind(object$par.h, object$se.h))
    rownames(estimates.h) <- c("Estimate:", "Std. Error:")
    colnames(estimates.h) <- object$names.h
    estimates <- cbind(estimates.g, estimates.h)
    print(round(estimates, 7))
    print(round(loglik, digits = 4))
    cat("\n")
  }
  if (is.null(object$order.g) || object$order.g[1] == 0) {
    cat("\n")
    model <- NULL
    if (object$order.h[3] != 0) {
      model <- paste(model, "GJR(", object$order.h[3], ")-", sep = "") 
    }
    if (object$order.h[1] != 0) model <- paste(model, "G", sep = "") 
    model <- paste(model, "ARCH(", sep = "") 
    if (object$order.h[1] != 0) model <- paste(model, object$order.h[1], ",", 
                                               sep = "") 
    model <- paste(model, object$order.h[2], ")", sep = "") 
    if (!is.null(object$xreg)) model <- paste(model, "-X", sep = "") 
    cat("Model:", model, "\n")
    cat("\n")
    if (object$turbo == TRUE) object$se.h <- sqrt(diag(vcov.tvgarch(object = 
                                                                      object)))
    estimates.h <- as.matrix(rbind(object$par.h, object$se.h))
    rownames(estimates.h) <- c("Estimate:", "Std. Error:")
    colnames(estimates.h) <- object$names.h
    print(round(estimates.h, 7))
    print(round(loglik, digits = 4))
    cat("\n")
    
  }
}

#####################################################
toLatex.tvgarch <- function (object, digits = 4, ...)
{
  if (!is.null(object$order.g) && object$order.g[1] != 0) {
    if (object$turbo == TRUE) {
      object$se.g <- sqrt(diag(vcov.tvgarch(object = object, spec = "tv")))
      object$se.h <- sqrt(diag(vcov.tvgarch(object = object, spec = "garch")))
    }
    else object$se.g <- object$se.g[-1]
    coef.h <- object$par.h
    if (object$order.h[1] != 0 && is.null(object$order.g)) names.h <- 
      c("", paste("y^{2}_{t-",paste(seq(1:object$order.h[2]), sep = ""),"}", 
                  sep = ""))
    if (object$order.h[1] != 0 && !is.null(object$order.g)) names.h <- 
      c("", paste("\\dfrac{y^{2}_{t-",paste(seq(1:object$order.h[2]), sep = ""),
                  "}}{\\widehat{g}_{t-",
                  paste(seq(1:object$order.h[2]), sep = ""),"}}", sep = ""))
    if (object$order.h[1] != 0 && !is.null(object$order.g)) names.h <- 
      c(names.h, paste("\\widehat{h}_{t-",paste(seq(1:object$order.h[1]), 
                                                sep = ""),"}", sep = "")) 
    if (object$order.h[1] != 0 && is.null(object$order.g)) names.h <- 
      c(names.h, paste("\\widehat{\\sigma}^2_{t-",
                       paste(seq(1:object$order.h[1]), sep = ""),"}", sep = "")) 
    if (object$order.h[3] != 0) names.h <- 
      c(names.h, paste("\\dfrac{y^{2}_{t-",paste(seq(1:object$order.h[3]), 
                                                 sep = ""),
                       "}}{\\widehat{g}_{t-", paste(seq(1:object$order.h[3]), 
                                                    sep = ""),
                       "}}\\text{I}(",paste("y_{t-",
                                            paste(seq(1:object$order.h[3]), 
                                                  sep = ""),"}", sep = ""),
                       " < 0)", sep = "")) 
    if (!is.null(object$xreg)) names.h <- 
      c(names.h, paste("x_{",paste(1:ncol(object$xreg),sep = ""),",t}", 
                       sep = ""))
    coefsNames.h <- names.h
    coef.h <- as.numeric(coef.h)
    stderrs <- as.numeric(object$se.h)
    eqtxt.h <- NULL
    for (i in 1:length(coef.h)) {
      ifpluss <- ifelse(i==1, "", " + ")
      eqtxt.h <- paste(eqtxt.h, ifelse(coef.h[i]<0, " - ", ifpluss), 
                       "\\underset{(", format(round(stderrs[i], digits = digits), 
                                              nsmall = digits), ")}{", 
                       format(round(abs(coef.h[i]), 
                                    digits = digits), nsmall = digits),"}", coefsNames.h[i], 
                       sep = "")
    }
    txtAddEq1 <- " \\\\[1mm]"
    txtAddEq2 <- " \\\\[1mm]"
    eqtxt.h <- paste0("  \\widehat{h}_t & = ", eqtxt.h, "", txtAddEq2, " \n")
    s <- length(object$order.g)
    coef.g <- object$par.g[1:(s+1)]
    coef.G <- object$par.g[-(1:(s+1))]
    coefsNames.g <- ""
    for (i in 1:s) {
      coefsNames.g <- c(coefsNames.g, paste("\\widehat{G}_{", paste(i, 
                                                                    sep = ""), 
                                            "}", sep = ""))
    }
    coef.g <- as.numeric(coef.g)
    stderrs.g <- as.numeric(object$se.g[1:s])
    stderrs.tv <- as.numeric(object$se.g[-(1:s)])
    eqtxt.g <- NULL
    for (i in 1:length(coef.g)) {
      ifpluss <- ifelse(i == 1, "", " + ")
      eqtxt.g <- paste(eqtxt.g, ifelse(coef.g[i] < 0, " - ",ifpluss))
      if (i == 1) eqtxt.g <- paste(eqtxt.g, "\\underset{(-)}{", 
                                   format(round(abs(coef.g[i]), 
                                                digits = digits), 
                                          nsmall = digits), "}", sep = "")
      else {
        eqtxt.g <- paste(eqtxt.g, "\\underset{(", format(round(stderrs.g[i-1], 
                                                               digits = digits), 
                                                         nsmall = digits), 
                         ")}{", format(round(abs(coef.g[i]), digits = digits), 
                                       nsmall = digits), "}",
                         "\\widehat{G}_{", paste(i-1, sep = ""), "}(", 
                         "\\underset{(", 
                         format(round(stderrs.tv[i-1], digits = digits), 
                                nsmall = digits), ")}{", 
                         format(round(abs(coef.G[i-1]), digits = digits),
                                nsmall = digits), "},", sep = "")
        for (j in 1:object$order.g[i-1]) {
          eqtxt.g <- paste(eqtxt.g, "\\underset{(", 
                           format(round(stderrs.tv[s+j], digits = digits), 
                                  nsmall = digits), ")}{", 
                           format(round(abs(coef.G[s+j]), digits = digits), 
                                  nsmall = digits), "}", sep = "")
          if (j < object$order.g[i-1]) eqtxt.g <- paste(eqtxt.g, ",", sep = "")
        }
        if (colnames(object$xtv) == "time") eqtxt.g <- paste(eqtxt.g, "; t/n)", 
                                                             sep = "")
        else eqtxt.g <- paste(eqtxt.g, "; s_{t})", sep = "")
      }
    }
    txtAddEq <- " \\\\[1mm]"
    eqtxt.g <- paste0("  \\widehat{g}_t & = ", eqtxt.g, "", txtAddEq, " \n")
    goftxt <- NULL
    goftxt <- "&"
    iT <- length(object$sigma2)
    goftxt <- paste(goftxt, " \\text{Log-likelihood: }", 
                    format(round(as.numeric(object$logLik), digits = digits), 
                           nsmall = digits), "\\qquad n = ", iT, 
                    " \n", sep = "")
    cat("%%the model was estimated:", object$date, "\n")
    cat("\\begin{align*}\n")
    cat(eqtxt.h)
    cat(eqtxt.g)
    cat(goftxt)
    cat("\\end{align*}\n")
  }
  if (is.null(object$order.g) || object$order.g[1] == 0) {
    coefs <- coef.tvgarch(object = object)
    coefsNames <- names(coefs)
    coefsNames[1] <- ""
    coefs <- as.numeric(coefs)
    if (object$turbo == TRUE) {
      object$se.h <- sqrt(diag(vcov.tvgarch(object = object)))
    }
    stderrs <- as.numeric(object$se.h)
    eqtxt <- NULL
    for (i in 1:length(coefs)) {
      ifpluss <- ifelse(i == 1, "", " + ")
      eqtxt <- paste(eqtxt, ifelse(coefs[i] < 0, " - ", ifpluss), 
                     "\\underset{(", format(round(stderrs[i], digits = digits), 
                                            nsmall = digits), ")}{", 
                     format(round(abs(coefs[i]), digits = digits), 
                            nsmall = digits), "}", coefsNames[i], sep = "")
    }
    txtAddEq <- " \\\\[1mm]"
    eqtxt <- paste0("  \\widehat{\\sigma}_t^2 &=& ", eqtxt, "", txtAddEq, " \n")
    goftxt <- NULL
    goftxt <- "   &&"
    iT <- nobs.tvgarch(object = object)
    goftxt <- paste(goftxt, " \\text{Log-likelihood: }", 
                    format(round(as.numeric(logLik.tvgarch(object = object)),
                                 digits = digits), nsmall = digits), 
                    "\\qquad n = ", iT, " \\nonumber \n", sep = "")
    cat("%%note: the 'eqnarray' environment requires the 'amsmath' package\n")
    cat("%%the model was estimated:", object$date, "\n")
    cat("\\begin{eqnarray}\n")
    cat(eqtxt)
    cat(goftxt)
    cat("\\end{eqnarray}\n")  }
}

#####################################################
vcov.tvgarch <- function (object, spec = c("tvgarch", "garch", "tv"), ...)
{
  if (!is.null(object$order.g) && object$order.g[1] != 0) { 
    spec <- match.arg(spec)
    npar.g <- length(object$par.g)-1
    if (object$turbo == TRUE) {
      if (spec == "tv") {
        s <- length(object$order.g)
        object$vcov.g <- matrix(NA_real_, npar.g, npar.g)
        jac.g <- jacobian(func = tvObj, x = object$par.g[-1], 
                          fixed.par.g = object$par.g[1], 
                          xtv = object$xtv, opt = object$opt, 
                          order.g = object$order.g, 
                          fixed.h = object$h, y = object$y, iter0 = FALSE, 
                          flag = 0)
        J.g <- crossprod(jac.g)  
        H.g <-  optimHess(par = object$par.g[-1], fn = tvObj, 
                          fixed.par.g = object$par.g[1], 
                          xtv = object$xtv, opt = object$opt, 
                          order.g = object$order.g, 
                          fixed.h = object$h, y = object$y, iter0 = FALSE, 
                          flag = 1)
        solHG <- solve(-H.g)
        vcov.g <- solHG %*% J.g %*% solHG
        object$vcov.g <- vcov.g
        rownames(object$vcov.g) <- object$names.g[-1]
        colnames(object$vcov.g) <- object$names.g[-1]
      } 
      if (spec == "garch") {
        object$vcov.h <- vcov.garchx(object = object$aux.h, 
                                     vcov.type = "robust")
      }
      if (spec == "tvgarch") {
        jac <- jacobian(func = tvgarchObj, x = c(object$par.g[-1],object$par.h),
                        fixed.par.g = object$par.g[1], y = object$y, 
                        order.g = object$order.g, xtv = object$xtv, 
                        opt = object$opt, iter.fit.h = object$iter.fit.h, 
                        flag = 0)
        JJ <- crossprod(jac)  
        HH <-  optimHess(par = c(object$par.g[-1],object$par.h), 
                         fn = tvgarchObj, fixed.par.g = object$par.g[1], 
                         y = object$y, order.g = object$order.g, 
                         xtv = object$xtv, opt = object$opt, 
                         iter.fit.h = object$iter.fit.h, flag = 1)
        solHG <- solve(-HH)
        object$vcov <- solHG %*% JJ %*% solHG
        colnames(object$vcov) <- c(object$names.g[-1], object$names.h)
        rownames(object$vcov) <- c(object$names.g[-1], object$names.h)
      }
    }
    if (spec == "tv") {
      return(object$vcov.g)
    }
    else if (spec == "garch") {
      return(object$vcov.h)
    }
    else if (spec == "tvgarch") {
      return(object$vcov)
    }
  }
  if (is.null(object$order.g) || object$order.g[1] == 0) {
    if (object$turbo == TRUE) {
      object$vcov.h <- vcov.garchx(object = object$aux.h, 
                                   vcov.type = "robust")
    }
    return(object$vcov.h)  
  }
}
