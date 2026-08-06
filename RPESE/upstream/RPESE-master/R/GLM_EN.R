# ====================================
# Functions - Penalized GLM
# Periodogram and RPEGLM Package Call
# ====================================

# Compute the periodogram as defined in H&W 1981
fit.periodogram <- function(data,
                            freq.include = c("All", "Decimate", "Truncate")[1], freq.par = 0.5,
                            twosided = FALSE, keep = 1){

  data.fft <- fft(data)
  N <- length(data)
  tmp <- Mod(data.fft[2:floor(N/2)])^2/N
  tmp <- sapply(tmp, function(x) max(0.00001,x))
  freq <- (1:(floor(N/2)-1))/N

  if(freq.include == "All"){

    tmp <- tmp[1:floor(length(tmp)*keep)]
    freq <- freq[1:floor(length(freq)*keep)]

  } else if(freq.include == "Decimate"){

    tmp <- tmp[seq(1, length(tmp), by = round(1/freq.par,0))]
    freq <- freq[seq(1, length(freq), by = round(1/freq.par,0))]

  } else if(freq.include == "Truncate"){
    tmp <- tmp[1:floor(length(tmp)*freq.par)]
    freq <- freq[1:floor(length(freq)*freq.par)]
  }

  if(twosided){
    tmp <- c(rev(tmp), tmp)
    freq <- c(-rev(freq), freq)
  }
  return(list(spec = tmp,freq = freq))
}

# Function to fit penalized GLM (Exponential or Gamma) to periodogram
SE.GLMEN <- function(data,
                     d = 7, alpha.EN = 0.5,
                     keep = 1,
                     standardize = FALSE,
                     prewhiten = FALSE, ar.coeffs = NULL,
                     fitting.method = c("Exponential", "Gamma")[1],
                     freq.include = c("All", "Decimate", "Truncate")[1], freq.par = 0.5,
                     return.coef = FALSE,
                     ...){

  # Computing the periodograms
  N <- length(data)
  my.periodogram <- fit.periodogram(data, keep  =  keep,
                                  freq.include = freq.include, freq.par = freq.par)
  my.freq <- my.periodogram$freq
  my.periodogram <- my.periodogram$spec
  nfreq <- length(my.freq)

  # creating the x.mat
  x.mat <- rep(1,length(my.freq))
  for(col.iter in 1:d){
    x.mat <- cbind(x.mat,my.freq^col.iter)
  }

  # Standardizing x.mat
  mean_vec <- apply(x.mat[,-1], 2, mean)
  sd_vec <- apply(x.mat[,-1], 2, sd)
  if (standardize)
    for(i in 2:ncol(x.mat)){
      tmp <- x.mat[,i]
      x.mat[,i] <- (tmp - mean(tmp))/sd(tmp)
    }

  # Testing whether the fitting.method vector is valid
  if(!(fitting.method %in% c("Exponential", "Gamma")))
    stop("The specified fitting method is not available.")

  # Fit the glmnet Gamma model
  if(fitting.method == "Gamma")
    GLM.out <- RPEGLMEN::fit.glmGammaNet(x.mat, my.periodogram, alpha.EN = alpha.EN, ...) else if(fitting.method == "Exponential")
      GLM.out <- RPEGLMEN::glmnet_exp(x.mat, my.periodogram, alpha.EN = alpha.EN, ...)


  # Return the estimated variance, and coeffs if return.coef  =  TRUE
  if(return.coef){
    if(standardize){
      variance <- exp(sum(GLM.out*c(1, -mean_vec/sd_vec)))/N
      if(prewhiten)
        variance <- variance/(1 - sum(ar.coeffs))^2
      coeffs <- GLM.out
      return(list(variance, coeffs))
    }
    variance <- exp(GLM.out[1])/N
    if(prewhiten)
      variance <- variance/(1 - sum(ar.coeffs))^2
    coeffs <- GLM.out
    return(list(variance, coeffs))
  }

  if(standardize){
    variance <- exp(sum(GLM.out*c(1,-mean_vec/sd_vec)))/N
    if(prewhiten)
      variance <- variance/(1 - sum(ar.coeffs))^2
    return(variance)
  }
  variance <- exp(GLM.out[1])/N
  if(prewhiten)
    variance <- variance/(1 - sum(ar.coeffs))^2
  return(variance)
}




