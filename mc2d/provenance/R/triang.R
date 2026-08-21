#' The Triangular Distribution
#' 
#' Density, distribution function, quantile function and random generation
#' for the triangular distribution with minimum equal to \samp{min}, mode equal \samp{mode}
#' (alternatively, mean equal \samp{mean}) and maximum equal to \samp{max}.
#' 
#' If \samp{min == mode == max}, there is no density in that case and 
#' \samp{dtriang} will return \samp{NaN} (the error condition) (Similarity with \code{\link[stats]{Uniform}}).
#'
#' \samp{mode} or \samp{mean} can be specified, but not both. Note: \samp{mean} is the last parameter for back-compatibility.
#' A warning will be provided if some combinations of \samp{min}, \samp{mean} and \samp{max} leads to impossible mode.
#' 
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If length(n) > 1, the length is taken to be the number required.
#' @param min vector of minima.
#' @param mode vector of modes.
#' @param max vector of maxima.
#' @param mean Vector of means, can be specified in place of \samp{mode} as an alternative parametrization.
#' @param log,log.p logical; if \samp{TRUE}, probabilities \samp{p} are given as \samp{log(p)}.
#' @param lower.tail logical; if \samp{TRUE} (default), probabilities are \samp{P[X <= x]}, otherwise, \samp{P[X > x]}.
#' 
#' @return
#' \samp{dtriang} gives the density, \samp{ptriang} gives the distribution function,
#' \samp{qtriang} gives the quantile function, and \samp{rtriang} generates random deviates.
#' @name triangular
#' @keywords distribution
#' @examples
#' curve(dtriang(x, min=3, mode=6, max=10), from = 2, to = 11, ylab="density")
#' ## Alternative parametrization
#' curve(dtriang(x, min=3, mean=6, max=10), from = 2, to = 11, add=TRUE, lty=2)
#' ##no density when  min == mode == max
#' dtriang(c(1,2,3),min=2,mode=2,max=2)
#' @export
dtriang <- function(x, min=-1, mode=0, max=1, log=FALSE, mean = 0){
	if(length(x) == 0) return(numeric(0))
  if (!missing(mode) && !missing(mean)) stop("specify 'mode' or 'mean' but not both")
  
  min <- as.vector(min)
  max <- as.vector(max)

  if (missing(mode)){
    mean <- as.vector(mean)
    mode <- (3*mean - min - max)
    if(any(mode < min | mode > max)) warning("Some values of mean lead to mode < min or mode > max.")
    
  } else {mode <- as.vector(mode)}
  
	# quel: x < mode or x = mode = max 
	xmaxmode <-  x >= max & mode >= max # to avoid some rounding errors
	quel <- (x < mode) | xmaxmode  
	d <- ifelse(quel,
              2*(x-min)/((mode-min)*(max-min)),
	            2 *(max-x)/((max-mode)*(max-min)))


	d[x < min | x > max] <- 0
	# For min = mode = max: provide an error like in dunif
	xminmodemax <- min >= max # to avoid some rounding errors
	d[xminmodemax] <- NaN
	
	d[mode < min | max < mode] <- NaN

	if(log) d <- log(d)
	if(any(is.na(d))) warning("NaN in dtriang")
  return(d)}

#' @rdname triangular
#' @export
ptriang <- function(q,min=-1,mode=0,max=1,lower.tail = TRUE, log.p = FALSE, mean = 0){
	if(length(q) == 0) return(numeric(0))
	# quel: q < mode or q = mode = max 
  if (!missing(mode) && !missing(mean)) stop("specify 'mode' or 'mean' but not both")
  
  min <- as.vector(min)
  max <- as.vector(max)
  
  if (missing(mode)){
    mean <- as.vector(mean)
    mode <- (3*mean - min - max)
    if(any(mode < min | mode > max)) warning("Some values of mean lead to mode < min or mode > max.")
    
  } else {mode <- as.vector(mode)}
  
  # quel: q < mode or q = mode = max 
  qmaxmode <-  q >= max & mode >= max # to avoid some rounding errors
	quel <- (q < mode) | qmaxmode  
	p <- ifelse(quel,
              (q-min)^2 / ((mode-min)*(max-min)),
	             1 - ((max-q)^2/((max-mode)*(max-min))))
	#if q = max = mode = min
	qminmodemax <- qmaxmode & q <= min
	p[qminmodemax] <- 1

	p[q < min] <- 0
	p[q > max] <- 1
	p[mode < min | max < mode] <- NaN
    if(!lower.tail) p <- 1-p
    if(log.p) p <- log(p)
	if(any(is.na(p))) warning("NaN in ptriang")
  return(p)}

#' @rdname triangular
#' @export
qtriang <- function(p, min=-1, mode=0, max=1, lower.tail=TRUE, log.p=FALSE, mean = 0){
  if (length(p) == 0) 
    return(numeric(0))
  if (!missing(mode) && !missing(mean)) stop("specify 'mode' or 'mean' but not both")
  
  min <- as.vector(min)
  max <- as.vector(max)
  
  if (missing(mode)){
    mean <- as.vector(mean)
    mode <- (3*mean - min - max)
    if(any(mode < min | mode > max)) warning("Some values of mean lead to mode < min or mode > max.")
    
  } else {mode <- as.vector(mode)}
  
  lout <- max(length(p),length(min),length(mode),length(max))
  min <- rep(min, length.out=lout)
  mode <- rep(mode, length.out=lout)
  max <- rep(max,length.out=lout)
  
  if (log.p) 
    p <- exp(p)
  if (!lower.tail) 
    p <- 1 - p
  quel <- p <= (mode - min)/(max - min)
  q <- ifelse(quel, min + sqrt(p * (mode - min) * (max - min)), 
              max - sqrt((1 - p) * (max - min) * (max - mode)))
  #if max = mode = min
  minmodemax <- min >= max # Avoid some rounding errors
  q <- ifelse(minmodemax, min, q)
  q[p < 0 | p > 1] <- NaN
  q[mode < min | max < mode] <- NaN
  if (any(is.na(q))) 
    warning("NaN in qtriang")
  return(q)
}


#' @rdname triangular
#' @export
rtriang <- function(n, min=-1, mode=0, max=1, mean = 0){
  if (length(n) > 1) 
    n <- length(n)
  if (length(n) == 0 || as.integer(n) == 0) 
    return(numeric(0))
  n <- as.integer(n)
  if (n < 0) 
    stop("integer(n) can not be negative in rtriang")
  if (!missing(mode) && !missing(mean)) stop("specify 'mode' or 'mean' but not both")
  
  min <- as.vector(min)
  max <- as.vector(max)
  
  if (missing(mode)){
    mean <- as.vector(mean)
    mode <- (3*mean - min - max)
    if(any(mode < min | mode > max)) warning("Some values of mean lead to mode < min or mode > max.")
    
  } else {mode <- as.vector(mode)}
  
  U <- runif(n)
  ow <- options(warn = -1)
  on.exit(options(ow), add = TRUE)
  q <- qtriang(U,
               min = rep(as.vector(min),length.out=n),
               mode = rep(as.vector(mode),length.out=n),
               max = rep(as.vector(max),length.out=n),
               lower.tail = TRUE,
               log.p = FALSE)
  if (any(is.na(q))) 
    warning("NaN in rtriang")
  return(q)
}
