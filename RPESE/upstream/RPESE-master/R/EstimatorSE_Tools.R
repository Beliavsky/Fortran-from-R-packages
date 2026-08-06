# ===========================
# EStimatorSE Tool Functions
# SE Types
# ===========================
#'
#' @importFrom utils packageDescription
#' @importFrom stats sd
#' @import xts
#' @import zoo
#'
# Compute the standard error for the xts object
SE.xts <-function(x, se.fun, myfun, myfun.IF,
                  prewhiten = FALSE, cleanOutliers = FALSE, fitting.method = c("Exponential", "Gamma")[1],
                  freq.include = c("All", "Decimate", "Truncate")[1], freq.par = 0.5, d.GLM.EN = 5,
                  ...){

  if (is.vector(x) || is.null(ncol(x)) || ncol(x) == 1) {
    x <- as.numeric(x)
    return(se.fun(x = x, myfun = myfun, myfun.IF = myfun.IF,
                  prewhiten = prewhiten, cleanOutliers = cleanOutliers, fitting.method = fitting.method,
                  freq.include = freq.include, freq.par = freq.par, d.GLM.EN = d.GLM.EN,
                  ...))
  }
  else{
    x <- coredata(x)
    return(apply(x, 2, se.fun, myfun = myfun, myfun.IF = myfun.IF,
                 prewhiten = prewhiten, cleanOutliers = cleanOutliers, fitting.method = fitting.method,
                 freq.include = freq.include, freq.par = freq.par, d.GLM.EN = d.GLM.EN,
                 ...))
  }
}

# Compute the standard error using influence function approach for a vector
SE.IF.iid <- function(x, myfun.IF, ...){
  N <- length(x)
  x.IF <- myfun.IF(x, ...)
  x.IF.2 <- x.IF^2
  fit.out <- mean(x.IF.2)
  return(sqrt(fit.out/N))
}

# Compute the standard error of the measure by iid bootstrapping
SE.BOOT.iid <- function(x, myfun, prewhiten = FALSE, nsim = 100, ...){

  res <- boot(data = x, statistic = function(x,i, ...) myfun(x[i], ...), R = nsim, ... = ...)
  return(sd(res$t))
}

# Compute the standard error of the measure by tsboot()
SE.BOOT.cor <- function(x, myfun, myfun.IF, prewhiten = FALSE, nsim = 1000,
                        sim = "fixed", l = round(length(x)/5), ...){

  res = tsboot(tseries = x, statistic = function(x, ...) myfun(x, ...), R = nsim,
               sim = sim, l = l, ...)
  return(sd(res$t))
}

# Compute the standard error using GLM-EN approach for serially correlated data using RPEGLMEN
SE.IF.cor <- function(x, myfun.IF,
                      d.GLM.EN = 5, alpha.EN = 0.5, keep = 1,
                      standardize = FALSE,
                      prewhiten = FALSE, cleanOutliers = FALSE, fitting.method = c("Exponential", "Gamma")[1],
                      freq.include = c("All", "Decimate", "Truncate")[1], freq.par = 0.5,
                      return.coef = FALSE,
                      ...){

  data.IF <- myfun.IF(x, prewhiten = FALSE, cleanOutliers = cleanOutliers, ...)
  if(prewhiten){
    ar.coeffs <- as.numeric(arima(x = data.IF, order = c(1,0,0), include.mean = TRUE)$coef[1])
    data.IF <- as.numeric(arima(x = data.IF, order = c(1,0,0), include.mean = TRUE)$res)
  } else{
    ar.coeffs <- NULL
  }
  fit.out <- SE.GLMEN(data.IF, standardize = standardize,
                      d = d.GLM.EN, alpha.EN = alpha.EN, keep = keep,
                      fitting.method = fitting.method,
                      freq.include = freq.include, freq.par = freq.par,
                      prewhiten = prewhiten, ar.coeffs = ar.coeffs,
                      return.coef = return.coef,
                      ...)
  if(return.coef){
    coeffs <- fit.out[[2]]
    fit.out <- fit.out[[1]]
    return(list(se = sqrt(fit.out), coef = coeffs))
  }
  return(sqrt(fit.out))
}

# Tool function for xts objects
sd.xts <- xts:::sd.xts
