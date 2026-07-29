coef.acdFit <- function(object, returnCoef = "all", ...){
  returnCoef <- match.arg(returnCoef, c("all", "distribution", "model"))
  
  switch(returnCoef,
                        all = c(object$mPara, object$dPara),
                        distribution = object$dPara,
                        model = object$mPara)
}

residuals.acdFit <- function(object, ...){
  object$residuals
}

predict.acdFit <- function(object, N = 10, ...){
  
  k <- max(object$order)
  endMu = utils::tail(object$muHats, k) #the end of the estimated expected durations of the fitted model
  if(length(object$durations$adjDur) != 0) 
    endDurations <- utils::tail(object$durations$adjDur, k) #the end of the durations of the fitted model
  else #no 'adjDur' column
    endDurations <- utils::tail(object$durations$durations, k) #the end of the durations of the fitted model
  
  errorExpectation <- 1
  #if the fitted model didn't have a forced error expectation = 1, the mean of the residuals is instead used:
  if(object$forceErrExpec == FALSE) errorExpectation <- mean(object$residuals) 
  
  #"simulates" with error terms equal to their expectation, starting from the endpoints of the original data set
  sim_ACD(N = N, model = object, Nburn = length(endDurations), startX = endDurations,
          startMu = endMu, errors = errorExpectation)
  
}

print.acdFit <- function(x, ...){
  
  if(x$distribution == "exponential") {
    cat("\nACD model estimation by (Quasi) Maximum Likelihood \n")
  } else {
    cat("\nACD model estimation by Maximum Likelihood \n")
  }
  
  cat("\nCall:\n") 
  cat(" ", deparse(x$call), "\n")
  cat("\nModel:\n") 
  cat(" ", x$model)
  cat("(")
  cat(x$order[1])
  for(i in 2:length(x$order)) cat("", x$order[i], sep = ", ")
  cat(")")
  if(length(x$SNIACDbp) != 0) cat("\n  Break points:", x$SNIACDbp)
  cat("\n")
  cat("\nDistribution:\n") 
  cat(" ", x$distribution)
  cat("\n\nN:", x$N) 
  cat("\n\nParameter estimate:\n") 
  print(format(x$parameterInference, digits = 3, scientific = F))
  if(length(x$comments) > 0){
    cat("\nNote:", x$comments) 
  }
  if(length(x$forcedDistPara) > 0){
    cat("\n\nFixed mean distribution parameter: \n") 
    cat(" ", names(x$forcedDistPara), ": ", x$forcedDistPara, sep = "")
  }
  if(length(x$bootErr) != 0){
    cat("\n\nBootstrap correlations:\n")
    print(format(data.frame(x$bootCorr), digits = 3, scientific = F))
  }
  if(length(x$robustCorr) != 0){
    cat("\n\nQML robust correlations:\n")
    print(format(data.frame(x$robustCorr), digits = 3, scientific = F))
  }
  cat("\n\nGoodness of fit:\n")
  print.data.frame(x$goodnessOfFit)
  cat("\nConvergence:", x$convergence, "\n")  
  cat("\nNumber of log-likelihood function evaluations:", x$evals, "\n")  
  if(length(x$bootErr) == 0) cat("\nEstimation time:", round(x$estimationTime, digits = 4), attributes(x$estimationTime)$units, "\n") 
  else cat("\nTotal estimation time (including bootstrap simulations):", round(x$estimationTime, digits = 4), attributes(x$estimationTime)$units, "\n") 
  cat("\nDescription:", x$description)
  cat("\n\n")
}

tidy.acdFit <- function(x, ...) {
  terms <- c(x$mPara, x$dPara)
  data.frame(
    term = names(terms),
    estimate = as.numeric(terms),         
    std.error = x$parameterInference$SE,
    p.value = x$parameterInference$PV
  )
}

glance.acdFit <- function(x, ...) {
  data.frame(
    logLik = x$goodnessOfFit$value[[1]],
    AIC = x$goodnessOfFit$value[[2]],
    BIC = x$goodnessOfFit$value[[3]],
    MSE = x$goodnessOfFit$value[[4]],
    nobs = x$N
  )
}
