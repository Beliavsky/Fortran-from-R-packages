## Functions for the skew-t distribution and copula #############

## 0. Helper functions from 'ghyp' ##################################

#' Computing the univariate df given the pdf
#' @param q vector of evaluation points
#' @param pdf function interpreted as the univariate pdf
#' @param pdf.args list of arguments for pdf()
#' @param lower lower integration values
#' @param upper upper integration values
#' @param ...
#' @return vector of probabilities
#' @note from 'ghyp' package
p_default_ <- function(q, pdf, pdf.args, lower, upper, ...){
   int.pars <- if(missing(upper)) list(f = pdf, lower = lower, upper = q) else
      list(f = pdf, lower = q, upper = upper)
   tmp.prob <- try(do.call("integrate", c(pdf.args, int.pars, list(...))),
                   silent = TRUE)
   if(inherits(tmp.prob, "try-error")){
      ## Use CMC
      W <- rigamma(1000, df = pdf.args$df)
      return(mean(pnorm( (q - pdf.args$gamma * W) / (sqrt(W)), lower.tail = TRUE)))
      #warning("Failed to determine probability with 'q = ", q,
      #        "'!\nMessage: ", as.character(tmp.prob), "\n")
      #return(NA)
   } else return(tmp.prob$value)
}

#' Computing the univariate quantile function given the pdf
#' @param p vector of probabilities
#' @param pdf function interpreted as the univariate pdf
#' @param pdf.args list of arguments for pdf()
#' @param interval interval containing the quantile
#' @param p.lower lower integration limit
#' @param upper upper integration values
#' @param ...
#' @return vector of quantiles
#' @note from 'ghyp' package
q_default_ <- function(p, cdf, interval, tol = 1e-7, cdf.args, ...){
   if(p > 0 & p < 1){
      dist.func <- function(x) cdf(x) - p
      tmp.quantile <- try(uniroot(dist.func, interval = interval, tol = tol))
      if(inherits(tmp.quantile, "try-error")){
         warning("Failed to determine quantile with 'probs = ", p,
                 "'!\nMessage: ", as.character(tmp.quantile), "\n")
         return(NA)
      } else return(tmp.quantile$root)
   } else if(p == 0) return(-Inf) else if(p == 1) return(Inf) else return(NA)
}




## 1. Inverse-gamma wrappers #########################################

qigamma <- function(p, df)
   1 / qgamma(1 - p, shape = df/2, rate = df/2)

rigamma <- function(n, df)
   1 / rgamma(n, shape = df/2, rate = df/2)

digamma <- function(x, df, log = FALSE){
   lres <- dgamma(1/x, shape = df/2, rate = df/2, log = TRUE) - 2 * log(x)
   if(log) lres else exp(lres)
}

## 2. Univariate skew-t distribution ###########################################

#' Standardize arguments for pskewtd1d()
#' @param p vector of evaluation points
#' @param sig scale parameter; > 0
#' @param gamma skewness parameter
#' @param loc location parameter
#' @author Erik Hintz
#' @return 2-list with standardized eval points 'p' and standardized 'gamma'
standardize_args <- function(q, sig, gamma, loc){
   if(loc != 0) q <- q - loc
   if(sig != 0) {
      q <- q/sig
      gamma <- gamma/sig
   }
   list(q = q, gamma = gamma) # return
}


#' Density of univariate skew-t distribution
#' @param x evaluation points
#' @param gamma skewness
#' @param df degrees-of-freedom
#' @param log logical wether log-density shall be returned
#' @return vector of (log-)density values
#' @author Erik Hintz using code from Yoshiba (2018)
#' @note More efficient than dskewt(...) in the univariate case
dskewt1d <- function(x, gamma, df, log = FALSE){
   if (!is.matrix(x))  x <- cbind(x) # 1-column matrix if 'x' is a vector
   d <- ncol(x) # dimension
   stopifnot(d == 1)
   mahax <- x * x
   mahagam <- gamma * gamma
   xSigmainvgam <- x * gamma
   dfQxgam <- (df + mahax) * mahagam
   ## Compute log-density
   ldens <- log(besselK(sqrt(dfQxgam), nu = (df+d)/2)) + xSigmainvgam +
      (df+d)/4*log(dfQxgam) - (df+d)/2*log(1 + mahax/df) + (1-(df+d)/2)*log(2) -
      lgamma(df/2) - d/2 * log(pi * df)
   ## Return
   if(log) ldens else exp(ldens)
}


## Distribution function of the univariate skew-t distribution
#' @param q vector of evaluation points
#' @param gamma skewness vector
#' @param df degrees-of-freedom, >0
#' @param log logical if log-density shall be returned
#' @author Erik Hintz using code from Yoshiba (2018)
#' @return n-vector of probabilities
qskewt1d_yos <- function(u, gamma, df, method = c("interpol", "integration", "splines"),
                        subdivisions = 100, spline.points = 200, interpol.points = 150,
                        root.tol = .Machine$double.eps^0.5, rel.tol = root.tol^1.5,
                        abs.tol = rel.tol, lower.tail = TRUE)
{
   if(!is.vector(u)) u <- as.vector(u)
   ## Trivial cases
   if(all(gamma == 0))
      if(is.infinite(df)) return(qnorm(u)) else return(qt(u, df = df))
   ## Set-up
   method <- match.arg(method)
   p.raw <- u
   p.raw[p.raw < 0 | p.raw > 1] <- NaN
   p.raw[p.raw == 1] <- Inf
   p.raw[p.raw == 0] <- -Inf
   p <- p.raw[is.finite(p.raw)]
   if(length(p) == 0) return(p.raw)
   ## Set up cdf
   mycdf <- function(x) pskewt1d_yos(q = x, gamma = gamma, df = df, loc = 0,
                                    sig = 1, subdivisions = subdivisions,
                                    rel.tol = rel.tol, abs.tol = abs.tol,
                                    lower.tail = TRUE)
   ## Helper function: Use Newton's method to find the range of the quantiles
   internal.bisection <- function(df, gamma, p, tol, rel.tol, abs.tol, subdivisions)
   {
      iter <- 0
      range.found <- FALSE
      tvar <- df/(df-2)
      step.size <- if(df > 4) sqrt(2*gamma^2*tvar^2/(df - 4)+ tvar) else 1
      if(!is.finite(step.size)) step.size <- 1
      q.0 <- qt(p, df = df)
      q.range <- c(q.0 - step.size, q.0 + step.size)
      while(!range.found & iter < 5000){
         iter <- iter + 1
         p.range <- pskewt1d_yos(q = q.range, gamma = gamma, df = df,
                                subdivisions = 100, rel.tol = rel.tol,
                                abs.tol = abs.tol) - p
         if(any(is.na(p.range))){
            warning("NA returned by pskewt1d()")
            return(NA)
         }
         lower <- p.range[1]
         upper <- p.range[2]
         if(upper < 0 & lower < 0){
            q.range[1] <- q.range[2]
            q.range[2] <- q.range[2] + step.size
            next
         }
         if(upper > 0 & lower > 0){
            q.range[2] <- q.range[1]
            q.range[1] <- q.range[1] - step.size
            next
         }
         if(upper > 0 & lower < 0){
            range.found <- TRUE
         }
      }
      if(iter >= 1000){
         warning("Unable to determine interval where the quantiles are in-between.\n",
                 "Perhaps the skewness 'gamma' is too large!")
      }
      q.root <- q_default_(p, mycdf, interval = q.range,
                           tol = root.tol)
      return(q.root)
   } # end of 'internal.bisection()'

   ## Actual computation
   if(length(p) == 1){
      ## If a single quantile is requested use the newton method anyway
      value <- internal.bisection(df, gamma = gamma, p = p, tol = root.tol,
                                  rel.tol = rel.tol, abs.tol = abs.tol,
                                  subdivisions = subdivisions)
      p.raw[is.finite(p.raw)] <- as.numeric(value)
      return(p.raw)
   }else if(length(p) == 2){
      ## If two quantiles are requested use the newton method anyway
      value1 <- internal.bisection(df, gamma = gamma, p = p[1], tol = root.tol,
                                   rel.tol = rel.tol, abs.tol = abs.tol,
                                   subdivisions = subdivisions)
      value2 <- internal.bisection(df, gamma = gamma, p = p[2], tol = root.tol,
                                   rel.tol = rel.tol, abs.tol = abs.tol,
                                   subdivisions = subdivisions)
      p.raw[is.finite(p.raw)] <- c(value1, value2)
      return(p.raw)
   }else{
      ## If more than two quantiles are requested use the newton method
      ## to find the range where the quantiles can be found.
      q.min <- internal.bisection(df, gamma = gamma, p = min(p), tol = root.tol,
                                  rel.tol = rel.tol, abs.tol = abs.tol,
                                  subdivisions = subdivisions)
      q.max <- internal.bisection(df, gamma = gamma, p = max(p), tol = root.tol,
                                  rel.tol = rel.tol, abs.tol = abs.tol,
                                  subdivisions = subdivisions)
      interval <- c(q.min, q.max)
      if(any(is.na(interval))){ # failed to determine bounds for the quantiles
         p.raw[is.finite(p.raw)] <- NA
         return(p.raw)
      }
      ## Extend the interval by 10 percent so that 'uniroot()' does not crash
      interval <- c(interval[1] - 0.1 * diff(range(interval)),
                    interval[2] + 0.1 * diff(range(interval)))
      if(method == "integration"){ # integration method
         pdf.args <- list(gamma = gamma, df = df, log = FALSE)
         p <- matrix(p, ncol = 1)
         value <- apply(p, MARGIN = 1, FUN = q_default_, cdf = mycdf,
                        interval = interval)

      } else if(method == "splines"){ # Splines method
         interval.seq <- seq(from = min(interval), to = max(interval),
                             length = spline.points)
         ## Compute the distribution function to be interpolated by splines
         p.interval <- mycdf(interval.seq)
         ## Spline function (as fast approximation to the df)
         spline.distribution.func <- splinefun(interval.seq, p.interval)
         ## Helper function: quantile.root.func == 0
         quantile.root.func <- function(x, tmp.p) spline.distribution.func(x) - tmp.p
         value <- p
         for(i in 1:length(p))
            value[i] <- uniroot(quantile.root.func, interval = interval,
                                tol = root.tol, tmp.p = p[i])$root
      } else {
         q.seq <- seq(from = q.min, to = q.max, length.out = interpol.points)
         ## Apply df to these interpolating quantiles
         px <- sort(mycdf(q.seq))
         ## Piecewise cubic hermite interpolation
         value <- pchip(px, q.seq, as.vector(p))
      }
      p.raw[is.finite(p.raw)] <- value
      p.raw
   }

}


## Distribution function of the univariate skew-t distribution
#' @param q vector of evaluation points
#' @param gamma skewness vector
#' @param df degrees-of-freedom, >0
#' @param log logical if log-density shall be returned
#' @author Erik Hintz using code from Yoshiba (2018)
#' @return n-vector of probabilities
pskewt1d_yos <- function(q, gamma, df, loc = 0, sig = 1, subdivisions = 100,
                        rel.tol = .Machine$double.eps^0.5, abs.tol = rel.tol,
                        lower.tail = TRUE)
{
   ## Standardize
   newargs <- standardize_args(q, sig = sig, gamma = gamma, loc = loc)
   q <- newargs$q
   gamma <- newargs$gamma
   ## Special case of t and normal distributions
   if(all(gamma == 0)){ # no skewness
      if(df == Inf){
         return(pnorm(q, lower.tail = lower.tail))
      }else{
         return(pt(q, df = df, lower.tail = lower.tail))
      }
   }
   q.raw <- q
   q.finite <- q.raw[is.finite(q.raw)]
   q.mat <- matrix(q.finite, ncol = 1)
   p.raw <- rep(NA, length(q.raw))
   pdf.args <- list(gamma = gamma, df = df, log = FALSE)
   if(lower.tail){
      p.raw[q.raw == -Inf] <- 0
      p.raw[q.raw == Inf] <- 1
      value <- apply(q.mat, 1, p_default_, pdf = "dskewt1d", lower = -Inf,
                     pdf.args = pdf.args, subdivisions = subdivisions,
                     rel.tol = rel.tol, abs.tol = abs.tol)
   }else{
      p.raw[q.raw == -Inf] <- 1
      p.raw[q.raw == Inf] <- 0
      value <- apply(q.mat, 1, p_default_, pdf = "dskewt1d", upper = Inf,
                     pdf.args = pdf.args, subdivisions = subdivisions,
                     rel.tol = rel.tol, abs.tol = abs.tol)
   }
   p.raw[is.finite(q.raw)] <- value
   ## Return
   p.raw
}


## 3. Multivariate skew-t distribution #########################################


## Random number generator of multivariate GH skew-t distribution
#' @param n sample size
#' @param loc location vector
#' @param scale scale matrix
#' @param factor *upper* triangular Cholesky factor of scale
#' @param gamma skewness vector
#' @param df degrees-of-freedom, >0
#' @author Erik Hintz using code from Yoshiba (2018)
#' @return (n, d) matrix of samples from t_d(df, loc, scale, gamma)
rskewt <- function(n, loc = rep(0, d), scale = diag(2), factor = NULL,
                   gamma = rep(0, d), df = Inf,
                   method = c("PRNG", "sobol", "ghalton"), skip = 0){
   ## Dimension
   if(is.null(factor)) factor <- chol(scale) # multiplication from the right below
   d <- nrow(factor)
   inversion <- (method != "PRNG")
   ## Compute W, Z
   if(inversion){
      U <- switch(method,
                  "sobol" = {
                     qrng::sobol(n, d = d + 1, randomize = "digital.shift",
                                 skip = skip)
                  },
                  "ghalton" = {
                     qrng::ghalton(n, d = d + 1, method = "generalized")
                  })
      W <- if(df == Inf) 1 else 1 / qgamma(1 - U[, 1], shape = df/2, rate = df/2)
      Z <- qnorm(U[, -1])
   } else {
      W <- if(df == Inf) 1 else 1/rgamma(n, shape = df/2, rate = df/2)
      Z <- matrix(rnorm(n * d), ncol = d)
   }
   Y <- Z %*% factor # scale-transform
   ## Return
   sweep((t(matrix(gamma, d, n)) * W + sqrt(W) * Y), 2, loc, "+")
}

## Random number generator of the skew-t copula
#' @param n sample size
#' @param scale scale matrix
#' @param factor *upper* triangular Cholesky factor of scale
#' @param gamma skewness vector
#' @param df degrees-of-freedom, >0
#' @param pseudo logical if pseudo-copula samples shall be returned
#' @author Erik Hintz using code from Yoshiba (2018)
#' @return (n, d) matrix of samples from C^t_(df, scale, gamma)
rskewtcopula <- function(n, scale = diag(2), factor = NULL,
                         gamma = rep(0, d), df = Inf, pseudo = TRUE,
                         method = c("PRNG", "sobol", "ghalton"), skip = 0){
   method <- match.arg(method)
   ## Dimension
   if(is.null(factor)) factor <- chol(scale) # multiplication from the right below
   d <- nrow(factor)
   ## Sample from the skew-t distribution
   X <- rskewt(n, factor = factor, gamma = gamma, df = df, method = method,
               skip = skip)
   ## Compute copula observations
   U <- if(pseudo) pobs(X) else {
      sapply(1:d, function(i) pskewt1d_yos(X[, i], gamma = gamma[i],
                                       df = df))
   }
   ## Return
   U
}

## Density of the multivariate skew-t distribution
#' @param x (n, d)-matrix of evaluation points
#' @param loc location vector
#' @param scale scale matrix
#' @param gamma skewness vector
#' @param df degrees-of-freedom, >0
#' @param log logical if log-density shall be returned
#' @author Erik Hintz using code from Yoshiba (2018)
#' @return n-vector of (log-)density values
#' @note If x is a vector, it is transformed to a 1-column matrix
dskewt <- function(x, loc = rep(0, d), scale = diag(2),
                   gamma = rep(0, d), df, log = FALSE, scale.inv, ldet){
   if (!is.matrix(x))
      x <- cbind(x) # 1-column matrix if 'x' is a vector
   d <- ncol(x) # dimension
   ## If length(gamma) == 1 => equiskewed; make a d-vector for mahalanobis()
   if(length(gamma) != d){
      stopifnot(length(gamma) == 1)
      gamma <- rep(gamma, d)
   }
   if(all(gamma == 0)) return(dStudent(
      x, df = df, scale = scale, gamma = gamma, log = log))
   ## Remove 'loc'
   if(any(loc != 0)) x <- sweep(x, 2, loc)
   ## Compute Sigma^{-1} along with log det(Sigma)
   if(!hasArg(scale.inv)){
      scale.inv <- pd.solve(scale, log.det = TRUE)
      ldet <- attr(scale.inv, "log.det")
   } else stopifnot(hasArg(ldet))
   ## Compute various quantities, such as maha-distances
   mahax <- mahalanobis(x, center = FALSE, cov = scale.inv,
                        inverted = TRUE)
   mahagam <- mahalanobis(gamma, center = FALSE, cov = scale.inv,
                          inverted = TRUE)
   xSigmainvgam <- as.vector(x %*% scale.inv %*% gamma)
   dfQxgam <- (df + mahax) * mahagam
   ## Compute log-density
   ldens <- log(max(besselK(sqrt(dfQxgam), nu = (df+d)/2), .Machine$double.xmin)) + xSigmainvgam +
      (df+d)/4*log(dfQxgam) - (df+d)/2*log(1 + mahax/df) + (1-(df+d)/2)*log(2) -
      lgamma(df/2) - d/2 * log(pi * df) - ldet/2
   ## Return
   if(log) ldens else exp(ldens)
}

#' Density of the skew-t copula
#' @param u (n, d)-matrix of evaluation points
#' @param scale scale matrix
#' @param gamma d-vector of skewness parameters
#' @param df dof parameter
#' @param log logical if log-density shall be returned
#' @param scale.inv inverse of scale; computed if not provided
#' @param ldet log(det(scale)); computed if not provided
#' @return n-vector with (log-)density values
#' @author Erik Hintz using code from Yoshiba (2018)
dskewtcopula <- function(u, scale = diag(2), gamma = rep(0, d), df, log = FALSE,
                         scale.inv, ldet){
   if (!is.matrix(u))
      u <- cbind(u) # 1-column matrix if 'u' is a vector
   d <- ncol(u) # dimension
   if(all(gamma == 0))
      return(dStudentcopula(u, df = df, scale = scale, log = log))
   ## Equiskewed?
   equiskewed <- all(gamma == gamma[1])
   pseudoskewt <- if(!equiskewed){
      sapply(1:d, function(i) qskewt1d_yos(u[, i], gamma = gamma[i], df = df))
   } else { # only *one* call to qskewt1d()
      matrix(qskewt1d_yos(as.vector(u), gamma = gamma[1], df = df), ncol = d)
   }
   lnum <- dskewt(pseudoskewt, scale = scale, gamma = gamma, df = df, log = TRUE,
                  scale.inv = scale.inv, ldet = ldet)
   ldenom <- rowSums(sapply(1:d, function(i)
      dskewt(pseudoskewt[, i], scale = as.matrix(1), df = df,
             gamma = gamma[i], log = TRUE)))
   ## Return
   if(log) lnum - ldenom else exp(lnum - ldenom)
}