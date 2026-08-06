# ================================
# RPE IF Plots
# IF Shape Evaluation
# ================================

# Shape evaluation 
EvaluateShape <- function(estimator, family, eff,
                          retVals, returns, k, nuisPars,
                          IFplot, IFprint, ...){
  
  # Function evaluation
  if(is.null(retVals))
    if(!is.null(returns))
      retVals <- seq(mean(returns)-k*sd(returns), mean(returns)+k*sd(returns), by=0.001) else
        retVals <- seq(0.005-k*0.07, 0.005+k*0.07, by=0.001)
  IFvals <- cbind(retVals, IF.fn(retVals, estimator = estimator, returns = returns, 
                                 nuisance.par = nuisPars, family = family, eff = eff, ...))
  colnames(IFvals) <- c("r", "IFvals")
  if(isTRUE(IFplot)){
    plot(IFvals[,1], IFvals[,2], type="l", 
         xlab="r", ylab="IF", col="blue", lwd=1, 
         main=estimator,
         panel.first=grid(), cex.lab=1.25)
    abline(h=0, v=0)
  }
  
  # Returning IF values
  return(IFvals)
}