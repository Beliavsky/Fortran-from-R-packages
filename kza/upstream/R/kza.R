.packageName <- "kza"

kz <- function(x,m,k = 3) {
    if (length(dim(x)) > 3) stop("Too many dimensions.")
    if (is.ts(x)) {
    	TS=TRUE
    	start=start(x)
    	f=frequency(x)
    } else {TS=FALSE}
    storage.mode(x) <- "double"
    x <- .Call("kz", x, as.vector(as.integer(floor(m/2))), as.integer(k))
    if (TS) {
    	x<-ts(x, start=start, frequency=f)
    }
    return (x)
}

# one 90-degree rotation; .rot90(.rot90(m, k), (4 - k) %% 4) is the identity
.rot90 <- function(m, k = 1) {
    k <- k %% 4
    for (i in seq_len(k)) m <- t(m)[ncol(m):1, , drop = FALSE]
    m
}

kza <- function(x, m, y = NULL, k = 3, min_size = round(0.05*m), tol = 1.0e-5, impute_tails = FALSE,
                symmetrize = FALSE, normalize = c("max", "quantile")) {
    if (length(dim(x)) > 3) stop("Too many dimensions.")
    if (is.null(y)) y<-kz(x,m=m,k=k)

    if (is.ts(x)) {
    	TS=TRUE
    	start=start(x)
    	f=frequency(x)
    } else {TS=FALSE}

    storage.mode(x) <- storage.mode(y) <- "double"
    # The shrink-factor normalizer: "max" is the published behavior (the
    # literal maximum of the difference metric); "quantile" uses its 99th
    # percentile so a handful of extreme pixels don't set the scale for
    # the whole image, falling back to the maximum whenever the quantile
    # collapses to zero (sparse structure on an otherwise flat field).
    normalize <- match.arg(normalize)
    nprob <- if (normalize == "quantile") 0.99 else 1.0
    # m is the (full) window width, matching kz(): the per-side radius
    # handed to the C code is floor(m/2), and the same radius drives the
    # kz() baseline above -- the published algorithm uses one q for both.
    # (Before 4.2.0 kza passed floor(m) as the radius, so the adaptive
    # window was ~double the baseline's for the same m.)
    run <- function(x., y.) .Call("kza", x., y., as.vector(as.integer(floor(m/2))), as.integer(k), as.integer(min_size), as.double(tol), as.double(nprob))
    if (symmetrize) {
        # The published head/tail window-allocation rule can pick the wrong
        # side immediately beside a sharp, high-contrast, axis-aligned edge,
        # so kza does not commute with 90-degree rotation. Averaging the
        # filter over the four rotations of the input cancels that
        # row/column bias -- exactly, for square input -- at 4x the compute.
        if (!is.matrix(x)) stop("symmetrize = TRUE is currently implemented for matrix input only")
        acc <- matrix(0, nrow(x), ncol(x))
        for (rk in 0:3) {
            acc <- acc + .rot90(run(.rot90(x, rk), .rot90(y, rk)), (4 - rk) %% 4)
        }
        kza.x <- acc / 4
    } else {
        kza.x <- run(x, y)
    }
    if (TS) {
    	kza.x<-ts(kza.x, start=start, frequency=f)
    }
    
    # Tail NA-marking is a 1D concept; applying it by linear index to a
    # matrix/array would blank the first/last m cells in column-major
    # order, at positions that mean nothing in 2D/3D.
    if (impute_tails==FALSE && is.null(dim(x))) {
    	kza.x[1:m]=NA
    	kza.x[(length(kza.x)-m):length(kza.x)]=NA
	}
    
 	structure(list(
 			time.series = x,
            kz = y,
            kza = kza.x,
            window=m, k=k, min_size=min_size, tol=tol,
            call=match.call()
            ),
        class = "kza")    
}

plot.kza <- function(x, ...)
{
	if (is.ts(x$kz) && is.ts(x$kza)) {
    	plot(cbind(
               kz = x$kz,
               kza = x$kza
               ),
         main = paste("KZA Decomposition of time series"), 
         ...)
	} else {
		par(mfrow=c(2,1))
		plot(x$kz, ylab="kz", type='l')
		plot(x$kza, ylab="kza", type='l')
		par(mfrow=c(1,1))
	}         
}

kzsv <- function(object) 
{
	if (inherits(object, "kza")) {y=object}
	else { stop("Need to use result from kza!") }
	
    # C signature: (kza_data, kz_data, window, minimum_window_length,
    # iterations, tolerance) -- iterations is accepted and unused. The old
    # call passed k as minimum_window_length and the nonexistent y$m as
    # iterations, and never passed min_size at all.
    # window is the full width; the C code takes the per-side radius
    s <- .Call("kzsv", y$kza, y$kz, as.integer(floor(y$window/2)), as.integer(y$min_size), as.integer(y$k), as.double(y$tol))

    if (is.ts(y$kz)) s<-ts(s, frequency=frequency(y$kz), start=start(y$kz))

 	structure(list(
 			kza = y,
            kzsv = s,
            call=match.call()
            ),
        class = "kzsv")    
}

plot.kzp <- function(x, ...)
{
	if (is.null(x$smooth_periodogram)) dz<-x$periodogram else dz<-x$smooth_periodogram
	omega<-seq(0:(length(x$periodogram)-1))/(x$k*(x$window-1))
	plot(omega, dz, type="l", xlab="Frequency", ylab="")
}

plot.kzsv <- function(x, ...)
{
	x<-x$kza
  plot(cbind(	kz = x$kz,	kza = x$kza, sigma= sqrt(x$kzsv/mean(x$kzsv))/1.96, 
    concavity=diff(diff(sqrt(x$kzsv/mean(x$kzsv))/1.96))), type='l',	    	
   		main = paste("KZSV Sample Variance"))
}

.cluster <- function(x, span=15)
{
	p=NULL
	m=NULL
	i=1
	z=x[i]
	for (y in x) {
		if (abs(y-z)<span) p=c(p,y) else  {m=c(m,round(mean(p))); p=c(y); } 
		z=x[i]
		i=i+1
	}
	m=c(m,round(mean(p)))
}

.peaks <- function(x, sigma=3, span=25)
{
	a<-x-mean(x)
	p2=NULL
	i=0

	sigma=sigma*sqrt(var(x))
	p<-a[a>sigma]
	p2<-rep(0,length(p))
	i=1
	for (j in p) {
		p2[i] <-	which(a==j)
		i=i+1
	}
	return (.cluster(p2, span))
}


summary.kzsv <- function(object, digits = getOption("digits"), ...)
{
    cat(" Call:\n ")
    dput(object$call, control=NULL)

	s<-sqrt(object$kzsv/mean(object$kzsv))/1.96
	d<-diff(diff(sqrt(object$kzsv/mean(object$kzsv))/1.96))
	
	p<-.peaks(s)
	
	if (is.ts(object$kzsv)) {
	    cat("\n Dates of interest:\n")
	    cat(" dates \t\t sigma\n")
	    for (m in p) {
	    	cat (" "); cat(as.integer(time(s)[m])); cat("\t");
	    	cat(" "); cat(cycle(s)[m]); cat("\t");
	    	cat(" "); cat(round(s[m],1)); cat("\n");
	    }
	} else {
	    cat("\n Periods of interest:\n")
	    cat(" period\n")
	    cat(" "); cat(p); cat("\t\t"); cat("\n");
	}
	
    invisible(object)
}

kzs <- function(y,m=NULL,k=3,t=NULL) 
{
	if (is.null(m)) {
		m=100*sqrt(mean(diff(y)^2))
		if (m>length(y)) m=2
	}		
	
	y<-kzft(y,m=m,k=k,f=0)
	a<-y

	if (is.ts(y)) ans<-ts(Re(a),start=start(y),frequency=frequency(y))
	else ans<-Re(a)
	return (ans)
}

