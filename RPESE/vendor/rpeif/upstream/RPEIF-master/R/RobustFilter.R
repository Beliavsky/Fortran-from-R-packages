# ================================
# Robust Filtering
# ================================

# Implementation of the robust filter for the IF TS
robust.cleaning <- function(returns, robust.method=c("locScaleRob")[1], 
                            alpha.robust=0.05, normal.efficiency=0.95){
  
  # Data cleaning using method from RobStatTM
  if(robust.method=="locScaleRob"){ 
    
    if(!requireNamespace("RobStatTM", quietly = TRUE)) {
      stop("Package \"RobStatTM\" needed for this function to work. Please install it.",
           call. = FALSE)
    }
  
    # Length of the IF TS
    n <- length(returns)
    # Casting as a numerci
    returns <- as.numeric(returns)
    # Location and scale parameters
    mu <- RobStatTM::locScaleM(returns, psi="mopt", eff=normal.efficiency)$mu
    s <- RobStatTM::locScaleM(returns, psi="mopt", eff=normal.efficiency)$disper
    
    # Efficiency parameter
    if(normal.efficiency==0.95)
      efficiency.param <- 3 else if(normal.efficiency==0.99)
        efficiency.param <- 3.568 else if(normal.efficiency==0.999)
          efficiency.param <- 4.21 else
            warning("Invalid value for the normal distribution efficiency.")
    
    # Computing the limits
    uplim <- rep(mu + s*efficiency.param, n) 
    dnlim <- rep(mu - s*efficiency.param, n)
    
    # Computing the filtered time series
    xcl <- ifelse(returns >= uplim, uplim, returns)
    xcl <- ifelse(returns <= dnlim, dnlim, xcl)
    
    return(xcl)
    
  } else
    warning("Invalid name for robust method.") # Otherwise invalid method name for robust data cleaning
  
}



