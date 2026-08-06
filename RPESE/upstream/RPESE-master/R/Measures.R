# ================
# RPE Estimators
# Point Estimates
# ================

# Downside SR
DSR <- function(returns, rf = 0, ...){

  # Computing the mean of the returns
  mu.hat <- mean(returns)
  # Computing the SemiSD
  semisd <- sqrt((1/length(returns))*sum((returns-mu.hat)^2*(returns <= mu.hat)))
  # Computing the SemiMean
  semimean <- (1/length(returns))*sum((returns-mu.hat)*(returns <= mu.hat))
  # Computing DSR of the returns
  DSR <- (mu.hat)/(semisd*sqrt(2))

  # Returning estimate
  return(DSR)
}

# ES
ES <- function(returns, alpha.ES = 0.05, ...){

  # Finding the quantile of the density fit based on the desired tail probability
  quantile.alpha <- as.numeric(quantile(returns, alpha.ES))
  # Computing the ES parameter based on the desired tail probability
  ES <- -mean(returns[returns <= quantile.alpha])

  # Returning estimate
  return(ES)
}

# ES Ratio
ESratio <- function(returns, alpha = 0.1, rf = 0, ...){

  # Computing the mean of the returns
  mu.hat <- mean(returns)
  # Computing the SD of the returns
  sigma.hat <- mean((returns-mu.hat)^2)
  # Storing the negative value of the VaR based on the desired alpha
  q.alpha <- as.numeric(quantile(returns, alpha))
  # Computing the ES of the returns
  ES <- -mean(returns[returns <= q.alpha])
  # Computing the ESratio of the returns
  ESratio <- (mu.hat-rf)/ES

  # Returning estimate
  return(ESratio)
}

# LPM
LPM <- function(returns, const = 0, order = 1, ...){

  # Returning estimate
  return(LPM(returns, const = const, order = 1, ...))
}

# Mean
Mean <- function(returns, ...){

  # Returning estimate
  return(mean(returns))
}

# Omega Ratio
OmegaRatio <- function(returns, const = 0, ...){

  # Computing Omega
  Omega <- 1 + (mean(returns)-const)/LPM(returns, const = const, order = 1)

  # Returning estimate
  return(Omega)
}

# Rachev Ratio
RachevRatio <- function(returns, alpha = 0.1, beta = 0.1, ...){

  # Finding the quantile of the density fit based on the desired lower tail probability
  q.alpha <- as.numeric(quantile(returns, alpha))
  # Computing the ES of the returns (lower tail)
  es.alpha <- -mean(returns[returns <= q.alpha])
  # Finding the quantile of the density fit based on the desired upper tail probability
  q.beta <- as.numeric(quantile(returns, 1-beta))
  # Computing the ES of the returns (upper tail)
  eg.beta <- mean(returns[returns >= q.beta])
  # Computing Rachev ratio
  RachevRatio <- eg.beta/es.alpha

  # Returning estimate
  return(RachevRatio)
}

# Robust Location M-estimator
robMean <- function(returns, family = c("mopt", "opt", "bisquare")[1], eff = 0.95, ...){

  # Computing robust location estimator
  robMean <- RobStatTM::locScaleM(returns, psi = family, eff = eff)$mu

  # Returning estimate
  return(robMean)
}

# Semi-SD
SemiSD <- function(returns, ...){

  # Computing the mean of the returns
  mu.hat <- mean(returns)
  # Computing the SemiSD
  SemiSD <- sqrt((1/length(returns))*sum((returns-mu.hat)^2*(returns <= mu.hat)))

  # Returning estimate
  return(SemiSD)
}

# SD
SD <- function(returns, ...){

  # Returning the estimate
  return(sd(returns))
}

# SoR
SoR <- function(returns, rf = 0, const = 0, threshold = c("mean", "const")[1], ...){

  # Case where we want mean threshold
  if(threshold == "mean")
    return(SoR_M(returns = returns, rf = rf)) else if(threshold == "const")
      return(SoR_C(returns = returns, const = const))
}

# SoR (Constant Threshold)
SoR_C <- function(returns, const = 0, ...){

  # Computation the mean of the returns
  mu.hat <- mean(returns)
  # Computating the LPM of the returns
  lpm2 <- LPM(returns, const = const, order = 2)
  # Computing the Sortino ratio of the returns
  SoR <- (mu.hat-const)/sqrt(lpm2)

  # Returning estimate
  return(SoR)
}

# SoR (Mean Threshold)
SoR_M <- function(returns, rf = 0, ...){

  # Computing the mean of the returns
  mu.hat <- mean(returns)
  # Computing the SemiSD
  semisd <- sqrt((1/length(returns))*sum((returns-mu.hat)^2*(returns <= mu.hat)))
  # Computing the SemiMean
  semimean <- (1/length(returns))*sum((returns-mu.hat)*(returns <= mu.hat))
  # Computing Sortino ratio of the returns
  SoR <- (mu.hat-rf)/(semisd*sqrt(2))

  # Returning estimate
  return(SoR)
}

# SR
SR <- function(returns, rf = 0, ...){

  # Computing the mean of the returns
  mu.hat <- mean(returns)
  # Computing the SD of the returns
  sd.hat <- sd(returns)
  # Computing the SR of the returns
  SR <- (mu.hat-rf)/sd.hat

  # Returning estimate
  return(SR)
}

# VaR
VaR <- function(returns, alpha = 0.05, ...){

  # Computing VaR
  VaR <- -quantile(returns, alpha)

  # Returning estimate
  return(VaR)
}

# VaR Ratio
VaRratio <- function(returns, alpha = 0.1, rf = 0, ...){

  # Mean of returns
  mu.hat <- mean(returns)
  # Fitting a density function to the returns
  density.fit <- approxfun(density(returns))
  # Probability for the obtained quantile
  fq.alpha <- quantile(returns, alpha)
  # Finding the quantile of the density fit based on the desired tail probability
  VaR <- -quantile(returns, alpha)
  # Computing the VaR ratio
  VaRratio <- (mu.hat - rf)/VaR

  # Returning estimate
  return(VaRratio)
}

# =================
# Function for LPM
# =================

LPM <- function(returns, const = 0, order = 1, ...){

  # Computing the LPM
  return(1/length(returns)*sum((const-returns[returns <= const])^order))
}

# =================
# Function for UPM
# =================

UPM <- function(returns, const = 0, order = 1, ...){

  # Computing the UPM
  return(1/length(returns)*sum((const-returns[returns >= const])^order))
}

