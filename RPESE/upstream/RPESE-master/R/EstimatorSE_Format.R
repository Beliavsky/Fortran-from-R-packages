# =================================
# EStimatorSE Formatting Function
# =================================

# Format output of EstimatorSE function
EstimatorSE.format <- function(SE.out, data, a, b, return.coef, se.method){

  # Adjusting the output for coefficients
  if(return.coef && (se.method == "IFcor" || se.method == "IFcorPW"))
    SE.out <- list(se = sapply(SE.out, function(l) l[[1]]), coef = sapply(SE.out, function(l) l[[2]])) else{
      if(se.method == "IFcorAdapt"){
        SE.out$cor <- list(se = sapply(SE.out$cor, function(l) l[[1]]))
        SE.out$corPW <- list(se = sapply(SE.out$corPW, function(l) l[[1]]))
      } else
        SE.out <- list(se = SE.out, coef = NULL)
    }

  # Adaptive method computation of weighted estimates
  if(se.method == "IFcorAdapt"){
    if(is.vector(data) || sum(ncol(data)>1) == 0){
      ar1.param <- arima(x = data, order = c(1,0,0), include.mean = TRUE)[[1]][1]
      if(0 <= ar1.param & ar1.param<a)
        w <- 0 else if(a <= ar1.param & ar1.param <= b)
          w <- (ar1.param - a)/(b - a) else
            w <- 1
          SE.out$corAdapt$se <- (1-w)*SE.out$cor$se + w*SE.out$corPW$se
          names(SE.out$corAdapt$se) <- colnames(data)
    } else{
      for(my.col in 1:ncol(data)){
        temp.data <- data[, my.col]
        ar1.param <- arima(x = temp.data, order = c(1,0,0), include.mean = TRUE)[[1]][1]
        if(0 <= ar1.param & ar1.param<a)
          w <- 0 else if(a <= ar1.param & ar1.param <= b)
            w <- (ar1.param - a)/(b - a) else
              w <- 1
        SE.out$corAdapt$se[my.col] <- (1-w)*SE.out$cor$se[my.col] + w*SE.out$corPW$se[my.col]
      }
      names(SE.out$corAdapt$se) <- colnames(data)
    }
    return(list(se = SE.out$corAdapt$se, coef = NULL))
  } else
    return(SE.out)
}
