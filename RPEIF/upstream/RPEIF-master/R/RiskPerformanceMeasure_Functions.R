# ==================
# RPE Estimators
# Functions for IF
# ==================

# Wrapper function for IF RPE
IF.fn <- function(x, estimator, returns, nuisance.par, family, eff, ...){
  
  # Available estimators
  estimator.available <- c("DSR", "ES", "ESratio", "LPM",
                           "Mean", "OmegaRatio", "RachevRatio", "robMean",
                           "SemiSD", "SD", "SR", "SoR", 
                           "VaR", "VaRratio")
  
  # Checking if the specified estimator is available
  if(!(estimator %in% estimator.available))
    stop("The specified estimator is not available.")
  
  # Plot for the specified estimator
  switch(estimator,
         DSR = IF.DSR.fn(x, returns, nuisance.par, ...),
         ES = IF.ES.fn(x, returns, nuisance.par, ...),
         ESratio = IF.ESratio.fn(x, returns, nuisance.par, ...),
         LPM = IF.LPM.fn(x, returns, nuisance.par, ...),
         Mean = IF.Mean.fn(x, returns, nuisance.par, ...),
         OmegaRatio = IF.OmegaRatio.fn(x, returns, nuisance.par, ...),
         RachevRatio = IF.RachevRatio.fn(x, returns, nuisance.par, ...),
         robMean = IF.robMean.fn(x, returns, nuisance.par, family, eff, ...),
         SemiSD = IF.SemiSD.fn(x, returns, nuisance.par, ...),
         SD = IF.SD.fn(x, returns, nuisance.par, ...),
         SR = IF.SR.fn(x, returns, nuisance.par, ...),
         SoR = IF.SoR.fn(x, returns, nuisance.par, ...),
         VaR = IF.VaR.fn(x, returns, nuisance.par, ...),
         VaRratio = IF.VaRratio.fn(x, returns, nuisance.par, ...))
}

# Downside SR IF Function
IF.DSR.fn <- function(x, returns, parsDSR.IF, rf = 0){
  
  # IF for null returns
  if(is.null(returns)){
    
    # Parameters for null returns
    mu.hat <- parsDSR.IF$mu
    semisd <- parsDSR.IF$semisd
    semimean <- parsDSR.IF$semimean
    dsr <- parsDSR.IF$dsr
    
  } else{
    
    # Computing the mean of the returns
    mu.hat <- mean(returns)
    # Computing the SemiSD
    semisd <- sqrt((1/length(returns))*sum((returns-mu.hat)^2*(returns <= mu.hat)))
    # Computing the SemiMean
    semimean <- (1/length(returns))*sum((returns-mu.hat)*(returns <= mu.hat))
    # Computing DSR of the returns
    dsr <- (mu.hat-rf)/(semisd*sqrt(2))
  }
  
  # IF computation
  IF.DSR <- (-dsr*(x - mu.hat)^2*(x <= mu.hat)/(2*semisd^2) + (x - mu.hat)*(1/semisd + dsr*semimean/semisd^2) + dsr/2)/sqrt(2)
  return(IF.DSR)
}

# ES IF Function
IF.ES.fn <- function(x, returns, parsES.IF, alpha.ES = 0.05){
  
  # IF for null returns
  if(is.null(returns)){
    
    # Finding the quantile of the density fit based on the desired tail probability
    quantile.alpha <- parsES.IF$q.alpha
    # Computing the ES parameter based on the desired tail probability
    ES.tail <- parsES.IF$es.alpha
    
  } else{
    
    # Finding the quantile of the density fit based on the desired tail probability
    quantile.alpha <- quantile(returns, alpha.ES)
    # Computing the ES parameter based on the desired tail probability
    ES.tail <- -mean(returns[returns <= quantile.alpha])
  }
  
  # IF Computation
  IF.ES <- 1/alpha.ES*((-x+quantile.alpha)*(x <= quantile.alpha))-quantile.alpha-ES.tail
  return(IF.ES)
}

# ES Ratio IF Function
IF.ESratio.fn <- function(x, returns, parsESratio.IF, alpha = 0.1, rf = 0){
  
  # IF for null returns
  if(is.null(returns)){
    
    # Parameters for null returns
    mu.hat <- parsESratio.IF$mu
    q.alpha <- parsESratio.IF$q.alpha
    ES.hat <- parsESratio.IF$es.alpha
    ESratio.hat <- parsESratio.IF$ES.ratio
    
  } else{
    
    # Computing the mean of the returns
    mu.hat <- mean(returns)
    # Computing the SD of the returns
    sigma.hat <- mean((returns-mu.hat)^2)
    # Storing the negative value of the VaR based on the desired alpha
    q.alpha <- as.numeric(quantile(returns, alpha))
    # Computing the ES of the returns
    ES.hat <- -mean(returns[returns <= q.alpha])
    # Computing the ESratio of the returns
    ESratio.hat <- (mu.hat-rf)/ES.hat
  }
  
  # IF computation 
  IF.ESratio <- (x - mu.hat)/ES.hat - ESratio.hat/ES.hat*((-x + q.alpha)*(x <= q.alpha)/alpha - q.alpha - ES.hat)
  return(IF.ESratio)
}

# LPM IF Function
IF.LPM.fn <- function(x, returns, parsLPM.IF, const = 0, order = 1){
  
  # Nuisance paramters if returns null
  if(is.null(returns)){

    lpm1 <- parsLPM.IF$lpm1
    lpm2 <- parsLPM.IF$lpm2
        
  } else{
    
    lpm1 <- LPM(returns, const = const, order = 1)
    lpm2 <- LPM(returns, const = const, order = 2)

  }
    
    # IF computation
    if(order == 1)
      IF.LPM <- (const - x)*(x <=  const) - lpm1 else if(order == 2)
        IF.LPM <- (const - x)^2*(x <=  const) - lpm2
    return(IF.LPM)
}

# Mean IF Function
IF.Mean.fn <- function(x, returns, parsMean.IF){
  
  # Computing the nuisance parameters
  if(is.null(returns))
    mu.hat <- parsMean.IF$mu else
      mu.hat <- mean(returns)
  
  # Computing the IF for the mean of the returns
  IF.Mean <- x - mu.hat
  
  # Return value for the estimator
  return(IF.Mean)
}

# Omega Ratio IF Function
IF.OmegaRatio.fn <- function(x, returns, parsOmega.IF, const = 0){
  
  # IF for null returns
  if(is.null(returns)){
    
    # Parameters for null returns
    lpm1 <- parsOmega.IF$lpm1
    upm1 <- parsOmega.IF$upm1
    omega <- parsOmega.IF$omega
    
  } else{
    
    # Returning length of returns vector
    N <- length(returns)
    # Computing Partial moments 
    lpm1 <- LPM(returns, const = const, order = 1) 
    upm1 <- UPM(returns, const = const, order = 1)
    omega <- 1 + (mean(returns)-const)/lpm1
  }
  
  # IF computation 
  IF.Omega <- 1/lpm1*((x-const)*(x >= const)-upm1) - omega/lpm1*((const-x)*(x <= const)-lpm1)
  return(IF.Omega)
}

# Rachev Ratio IF Function
IF.RachevRatio.fn <- function(x, returns, parsRachev.IF, alpha = 0.1, beta = 0.1){
  
  # IF for null returns
  if(is.null(returns)){
    
    # Parameters for null returns
    q.alpha <- parsRachev.IF$q.alpha
    es.alpha <- parsRachev.IF$es.alpha
    q.beta <- parsRachev.IF$q.beta
    eg.beta <- parsRachev.IF$eg.beta
    rachev.ratio <- parsRachev.IF$rach.r
    
  } else{
    
  
    # Finding the quantile of the density fit based on the desired lower tail probability
    q.alpha <- as.numeric(quantile(returns, alpha))
    # Computing the ES of the returns (lower tail)
    es.alpha <- -mean(returns[returns <= q.alpha])
    
    # Finding the quantile of the density fit based on the desired upper tail probability
    q.beta <- as.numeric(quantile(returns, 1-beta))
    # Computing the ES of the returns (upper tail)
    eg.beta <- mean(returns[returns >= q.beta])
    # Computing Rachev ratio
    rachev.ratio <- eg.beta/es.alpha
  }
  
  # IF computation 
  IF.RachevRatio <- (1/es.alpha)*((x >= q.beta)*(x-q.beta)/beta + q.beta - eg.beta) - 
    (rachev.ratio/es.alpha)*((-1)*(x <= q.alpha)*(x-q.alpha)/alpha - q.alpha - es.alpha)
  return(IF.RachevRatio)
}

# Rachev Ratio IF Function
IF.robMean.fn <- function(x, returns, parsrobMean.IF, family, eff){
  
  # Extracting tuning parameters
  tuning.parameters <- 
    switch(family,
         mopt = RobStatTM::mopt(eff),
         opt = RobStatTM::opt(eff),
         bisquare = RobStatTM::bisquare(eff))
  
  if(is.null(returns)){
    
    rob.mu <- parsrobMean.IF$mu
    rob.sd <- parsrobMean.IF$sd
    psi.prime <- parsrobMean.IF$psi.prime
    
    # IF computation 
    IF.robMean <- rob.sd * (RobStatTM::rhoprime((x - rob.mu) / rob.sd, family = family, cc = tuning.parameters)) / 
      psi.prime

  } else{
    
    # Computing the robust estimates for location and scale
    rob.mu <- RobStatTM::locScaleM(returns, psi = family, eff = eff)$mu
    rob.sd <- mad(returns)
    
    # IF computation 
    IF.robMean <- rob.sd * (RobStatTM::rhoprime((x - rob.mu) / rob.sd, family = family, cc = tuning.parameters)) / 
      mean(RobStatTM::rhoprime2((returns - rob.mu) / rob.sd, family = family, cc = tuning.parameters))
  }

  return(IF.robMean)
}

# Semi-SD IF Function
IF.SemiSD.fn <- function(x, returns, parsSemiSD.IF){
  
  # IF for null returns
  if(is.null(returns)){
    
    mu.hat <- parsSemiSD.IF$mu
    semisd <- parsSemiSD.IF$semisd
    semimean <- parsSemiSD.IF$semimean
    
  } else{
    
    # Computing the mean of the returns
    mu.hat <- mean(returns)
    # Computing the SemiSD 
    semisd <- sqrt((1/length(returns))*sum((returns-mu.hat)^2*(returns<= mu.hat)))
    # Computing the SemiMean
    semimean <- (1/length(returns))*sum((returns-mu.hat)*(returns<= mu.hat))
  }
  
  # IF computation
  IF.SemiSD <- ((x - mu.hat)^2 * (x <=  mu.hat) - 2*semimean * (x-mu.hat) - semisd^2)/(2*semisd)
  return(IF.SemiSD)
}

# SD IF Function
IF.SD.fn <- function(x, returns, parSemiSD.IF){
  
  # Computing the nuisance parameters
  if(is.null(returns)){
    
    mu.hat <- parSemiSD.IF$mu
    sd.hat <- parSemiSD.IF$sd
    
  } else{
    
    # Computing the mean of the returns
    mu.hat <- mean(returns)
    # Computing the standard deviation of the returns
    sd.hat <- sd(returns)
  }
  
  # IF computation
  IF.SD <- ((x-mu.hat)^2-sd.hat^2)/2/sd.hat
  return(IF.SD)
}

# SoR IF Function
IF.SoR.fn <- function(x, returns, parsSoR.IF, rf = 0, const = 0, threshold = c("mean", "const")[1]){
  
  # Case where we want mean threshold
  if(threshold == "mean")
    return(IF.SoR_M.fn(x = x, returns = returns, parsSoR_M.IF = parsSoR.IF, rf = rf)) else if(threshold == "const")
      return(IF.SoR_C.fn(x = x, returns = returns, parsSoR_C.IF = parsSoR.IF, const = const))
}

# SoR (Constant Threshold) IF Function
IF.SoR_C.fn <- function(x, returns, parsSoR_C.IF, const = 0){
  
  # IF for null returns
  if(is.null(returns)){
    
    # Parameters for null returns
    mu.hat <- parsSoR_C.IF$mu
    sor.c <- parsSoR_C.IF$sor.c
    lpm2 <- parsSoR_C.IF$lpm2
    
  } else{
    
    # Computation the mean of the returns
    mu.hat <- mean(returns)
    # Computating the LPM of the returns
    lpm2 <- LPM(returns, const = const, order = 2)
    # Computing the Sortino ratio of the returns
    sor.c <- (mu.hat-const)/sqrt(lpm2)
    
  }
  
  # IF computation
  IF.SoR <- -sor.c*(x-const)^2*(x<= const)/(2*lpm2) + (x-mu.hat)/sqrt(lpm2) + sor.c/2
  return(IF.SoR)
}

# SoR (Mean Threshold) IF Function
IF.SoR_M.fn <- function(x, returns, parsSoR_M.IF, rf = 0){
  
  # IF for null returns
  if(is.null(returns)){
    
    # Parameters for null returns
    mu.hat <- parsSoR_M.IF$mu
    semisd <- parsSoR_M.IF$semisd
    semimean <- parsSoR_M.IF$semimean
    sor <- parsSoR_M.IF$sor.mu
    
  } else{
    
    # Computing the mean of the returns
    mu.hat <- mean(returns)
    # Computing the SemiSD
    semisd <- sqrt((1/length(returns))*sum((returns-mu.hat)^2*(returns<= mu.hat)))
    # Computing the SemiMean
    semimean <- (1/length(returns))*sum((returns-mu.hat)*(returns<= mu.hat))
    # Computing Sortino ratio of the returns
    sor <- (mu.hat-rf)/(semisd)
  }
  
  # IF computation
  IF.SoR <- (-sor*(x - mu.hat)^2*(x<= mu.hat)/(2*semisd^2) + (x - mu.hat)*(1/semisd + sor*semimean/semisd^2) + sor/2)
  return(IF.SoR)
}

# SR IF Function
IF.SR.fn <- function(x, returns, parsSR.IF, rf = 0){
  
  # Computing the nuisance parameters
  if(is.null(returns)){
    
    mu.hat <- parsSR.IF$mu
    sd.hat <- parsSR.IF$sd
    sr <- parsSR.IF$sr
    
  } else{
    
    # Computing the mean of the returns
    mu.hat <- mean(returns)
    # Computing the SD of the returns
    sd.hat <- sd(returns)
    # Computing the SR of the returns
    sr <- (mu.hat-rf)/sd.hat
  }
  
  # IF Computation
  IF.SR <- 1/sd.hat*(x-mu.hat)-1/2*mu.hat/sd.hat^3*((x-mu.hat)^2-sd.hat^2)
  return(IF.SR)
}

# VaR IF Function
IF.VaR.fn <- function(x, returns, parsVaR.IF, alpha = 0.05){
  
  # IF for null returns
  if(is.null(returns)){
    
    # Fitting a density function to the returns
    fq.alpha <- parsVaR.IF$fq.alpha
    # Finding the quantile of the density fit based on the desired tail probability
    quantile.alpha <- parsVaR.IF$q.alpha
    
  } else{
    
    # Fitting a density function to the returns
    density.fit <- approxfun(density(returns))
    # Finding the quantile of the density fit based on the desired tail probability
    quantile.alpha <- quantile(returns, alpha)
    # Probability for the obtained quantile
    fq.alpha <- density.fit(quantile.alpha)
  }
  
  # IF Computation
  IF.VaR <- ((x<= quantile.alpha)-alpha)/fq.alpha
  return(IF.VaR)
}

# VaR Ratio IF Function
IF.VaRratio.fn <- function(x, returns, parsVaRratio.IF, alpha = 0.1, rf = 0){
  
  # IF for null returns
  if(is.null(returns)){
    
    # Parameters for null returns
    mu.hat <- parsVaRratio.IF$mu
    VaR.hat <- -parsVaRratio.IF$q.alpha
    VaRratio.hat <- parsVaRratio.IF$VaR.ratio
    fq.alpha <- parsVaRratio.IF$fq.alpha

  } else{
    
    # Mean of returns
    mu.hat <- mean(returns)
    # Fitting a density function to the returns
    density.fit <- approxfun(density(returns))
    # Probability for the obtained quantile
    fq.alpha <- quantile(returns, alpha)
    # Finding the quantile of the density fit based on the desired tail probability
    VaR.hat <- -quantile(returns, alpha)
    # Computing the VaR ratio
    VaRratio.hat <- (mu.hat - rf)/VaR.hat
  }
  
  # IF Computation
  IF.VaRratio <- (x - mu.hat)/VaR.hat - (VaRratio.hat/VaR.hat) * ((x <=  -VaR.hat)-alpha)/fq.alpha
  return(IF.VaRratio)
}







