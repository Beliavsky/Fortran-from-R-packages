#' The (Modified) PERT Distribution
#'
#' Density, distribution function, quantile function and random generation
#' for the PERT (\emph{aka} Beta PERT) distribution with minimum equals to \samp{min}, mode equals to \samp{mode}
#' (or, alternatively, mean equals to \samp{mean}) and maximum equals to \samp{max}. 
#'
#' The PERT distribution is a \code{\link{Beta}} distribution extended to the domain \samp{[min, max]} with mean 
#' \deqn{mean=\frac{min+shape\times mode+max}{shape+2}}{mean = (min + shape * mode + max)/(shape + 2)}
#'
#' The underlying beta distribution is specified by \eqn{\alpha_{1}}{shape1} and \eqn{\alpha_{2}}{shape2} defined as
#'
#' \deqn{\alpha_{1}=\frac{(mean-min)(2\times mode-min-max)}{(mode-mean)(max-min)}}{shape1=(mean - min)*(2 mode-min-max)/((mode-mean)*(max - min))}
#'
#' \deqn{\alpha_{2}=\frac{\alpha_{1}\times (max-mean)}{mean-min}}{shape2=shape1*(max - mean)/(mean - min)}
#'
#' \samp{mode} or \samp{mean} can be specified, but not both. Note: \samp{mean} is the last parameter for back-compatibility. 
#' A warning will be provided if some combinations of \samp{min}, \samp{mean} and \samp{max} leads to impossible mode.
#' 
#' David Vose (See reference) proposed a modified PERT distribution with a shape parameter different from 4. 
#'
#' The PERT distribution is frequently used (with the \link{triangular} distribution) to translate expert estimates
#' of the min, max and mode of a random variable in a smooth parametric distribution. 
#' @param x,q Vector of quantiles.
#' @param p Vector of probabilities
#' @param n Number of observations. If length(n) > 1, the length is taken to be the number required.
#' @param min Vector of minima.
#' @param mode Vector of modes.
#' @param mean Vector of means, can be specified in place of \samp{mode} as an alternative parametrization.
#' @param max Vector of maxima.
#' @param shape Vector of scaling parameters. Default value: 4.
#' @param log,log.p Logical; if \samp{TRUE}, probabilities \samp{p} are given as \samp{log(p)}.
#' @param lower.tail Logical; if \samp{TRUE} (default), probabilities are \samp{P[X <= x]}, otherwise, \samp{P[X > x]}
#' @name pert
#' @keywords distribution
#' @references Vose D. Risk Analysis - A Quantitative Guide (2nd and 3rd editions, John Wiley and Sons, 2000, 2008).  
#' @author Regis Pouillot and Matthew Wiener
#' @examples
#' curve(dpert(x,min=3,mode=5,max=10,shape=6), from = 2, to = 11, lty=3,ylab="density")
#' curve(dpert(x,min=3,mode=5,max=10), from = 2, to = 11, add=TRUE)
#' curve(dpert(x,min=3,mode=5,max=10,shape=2), from = 2, to = 11, add=TRUE,lty=2)
#' legend(x = 8, y = .30, c("Default: 4","shape: 2","shape: 6"), lty=1:3)
#' ## Alternatie parametrization using mean
#' curve(dpert(x,min=3,mean=5,max=10), from = 2, to = 11, lty=2 ,ylab="density")
#' curve(dpert(x,min=3,mode=5,max=10), from = 2, to = 11, add=TRUE)
#' legend(x = 8, y = .30, c("mode: 5","mean: 5"), lty=1:2)
#' @seealso \code{\link{Beta}}
#' @return
#' \samp{dpert} gives the density, \samp{ppert} gives the distribution function,
#' \samp{qpert} gives the quantile function, and \samp{rpert} generates random deviates.
#' @export
dpert <- function(x,min=-1,mode=0,max=1,shape=4,log=FALSE, mean = 0){
	if(length(x) == 0) return(numeric(0))
  if (!missing(mode) && !missing(mean)) stop("specify 'mode' or 'mean' but not both")
  
  min <- as.vector(min)
	max <- as.vector(max)
	shape <- as.vector(shape)

  if (missing(mode)){
    mean <- as.vector(mean)
    mode <- ((shape+2)*mean - min - max) / shape
    if(any(mode < min | mode > max)) warning("Some values of mean lead to mode < min or mode > max.")
    
  } else {mode <- as.vector(mode)}
  
	if(any(mode < min | mode > max)) warning("Some values of mode < min or mode > max.")
	
	a1 <- 1+shape*(mode-min)/(max-min)	
	a2 <- 1+shape*(max-mode)/(max-min)
	
	oldw <- options(warn = -1)
	on.exit(options(oldw), add = TRUE)
	d <- (x-min)^(a1-1) * (max-x)^(a2-1) / beta(a=a1,b=a2) / (max-min)^(a1+a2-1)
	
  d[x < min | x > max] <- 0
	d[mode < min | max < mode] <- NaN
	d[shape <= 0] <- NaN
	if(log) d <- log(d)
	if(any(is.na(d))) warning("NaN in dpert")
  return(d)}

#' @rdname pert
#' @export
ppert <- function(q,min=-1,mode=0,max=1,shape=4,lower.tail = TRUE, log.p = FALSE, mean = 0){
	if(length(q) == 0) return(numeric(0))
  if (!missing(mode) && !missing(mean)) stop("specify 'mode' or 'mean' but not both")
  
  min <- as.vector(min)
  max <- as.vector(max)
  shape <- as.vector(shape)
  
  if (missing(mode)){
    mean <- as.vector(mean)
    mode <- ((shape+2)*mean - min - max) / shape
    if(any(mode < min | mode > max)) warning("Some values of mean lead to mode < min or mode > max.")
    
  } else {mode <- as.vector(mode)}
  
  if(any(mode < min | mode > max)) warning("Some values of mode < min or mode > max.")
  
  a1 <- 1+shape*(mode-min)/(max-min)	
	a2 <- 1+shape*(max-mode)/(max-min)
  
	oldw <- options(warn = -1)
	on.exit(options(oldw), add = TRUE)
	p <- pbeta(q=(q-min)/(max-min),shape1=a1,shape2=a2)
	
  p[q < min] <- 0
	p[q >= max] <- 1
	p[mode < min | max < mode] <- NaN
	p[shape <= 0] <- NaN
	
  if(!lower.tail) p <- 1-p
	if(log.p) p <- log(p)
	if(any(is.na(p))) warning("NaN in ppert")
  return(p)}

#' @rdname pert
#' @export
qpert <- function(p, min=-1, mode=0, max=1, shape=4, lower.tail=TRUE, log.p=FALSE, mean=0){
    if (length(p) == 0) 
      return(numeric(0))
    if (!missing(mode) && !missing(mean)) stop("specify 'mode' or 'mean' but not both")
    
    min <- as.vector(min)
    max <- as.vector(max)
    shape <- as.vector(shape)
    
    if (missing(mode)){
      mean <- as.vector(mean)
      mode <- ((shape+2)*mean - min - max) / shape
      if(any(mode < min | mode > max)) warning("Some values of mean lead to mode < min or mode > max.")
      
    } else {mode <- as.vector(mode)}
  
    if(any(mode < min | mode > max)) warning("Some values of mode < min or mode > max.")
    
    lout <- max(length(p),length(min),length(mode),length(max),length(shape))
    min <- rep(min, length.out=lout)
    mode <- rep(mode, length.out=lout)
    max <- rep(max, length.out=lout)
    shape <- rep(shape, length.out=lout)
    
    if (log.p) 
      p <- exp(p)
    if (!lower.tail) 
      p <- 1 - p
    a1 <- 1 + shape * (mode - min)/(max - min)
    a2 <- 1 + shape * (max - mode)/(max - min)
    oldw <- options(warn = -1)
    on.exit(options(oldw), add = TRUE)
    q <- qbeta(p, shape1 = a1, shape2 = a2)
    q <- q * (max - min) + min
    # Check if min == max, but avoid some of the rounding errors
    minmodemax <- min >= max
    q <- ifelse(minmodemax, min, q)
    q[p < 0 | p > 1] <- NaN
    q[mode < min | max < mode] <- NaN
    q[shape <= 0] <- NaN
    if (any(is.na(q))) 
      warning("NaN in qpert")
    return(q)
  }



#' @rdname pert
#' @export
rpert <- function(n,min=-1,mode=0,max=1,shape=4, mean=0){
  if (length(n) > 1) 
    n <- length(n)
  if (length(n) == 0 || as.integer(n) == 0) 
    return(numeric(0))
  n <- as.integer(n)
  
  if (!missing(mode) && !missing(mean)) stop("specify 'mode' or 'mean' but not both")

  min <- rep(as.vector(min),length.out=n)
  max <- rep(as.vector(max),length.out=n)
  shape <- rep(as.vector(shape),length.out=n)
  
  if (missing(mode)){
    mean <- rep(as.vector(mean),length.out=n)
    mode <- ((shape+2)*mean - min - max) / shape
    if(any(mode < min | mode > max)) warning("Some values of mean lead to mode < min or mode > max.")
    
  } else {  mode <- rep(as.vector(mode),length.out=n) }
  

  a1 <- 1 + shape * (mode - min)/(max - min)
  a2 <- 1 + shape * (max - mode)/(max - min)
  oldw <- options(warn = -1)
  on.exit(options(oldw), add = TRUE)
  r <- rbeta(n, shape1 = a1, shape2 = a2) * (max - min) + min
  # Check if min == max, but avoid some of the rounding errors
  minmodemax <- min >= max
  r <- ifelse(minmodemax, min, r)

  check <- mode < min | mode > max
  if(any(check)) {
    r[check] <- NA
    warning("Some values of mode < min or mode > max.")
  }
  
  if (any(is.na(r))) 
    warning("NaN in rpert")
  return(r)
}
