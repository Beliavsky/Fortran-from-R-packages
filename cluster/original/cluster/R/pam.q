#### PAM : Partitioning Around Medoids
#### --- $Id: pam.q 8260 2023-08-11 20:10:20Z maechler $
pam <- function(x, k, diss = inherits(x, "dist"),
		metric = c("euclidean", "manhattan"), ## FIXME: add "jaccard"
                medoids = if(is.numeric(nstart)) "random",
                nstart = if(variant == "faster") 1L else NA,
                stand = FALSE, cluster.only = FALSE, do.swap = TRUE,
                keep.diss = !diss && !cluster.only && n < 100,
                keep.data = !diss && !cluster.only,
                variant = c("original", "o_1", "o_2", "f_3", "f_4", "f_5", "faster"),
		pamonce = FALSE, trace.lev = 0)
{
    stopifnot(length(cluster.only) == 1, length(trace.lev) == 1)
    nMax <- 65536 # 2^16 (as 1+ n(n-1)/2 must be < max_int = 2^31-1)
    if((diss <- as.logical(diss))) {
	## check type of input vector
	if(anyNA(x)) stop("NA values in the dissimilarity matrix not allowed.")
        if(keep.data) stop("Cannot keep data when 'x' is a dissimilarity!")
	if(!inherits(x, "dissimilarity")) { # try to convert to
	    if(!is.null(dim(x))) {
		x <- as.dist(x) # or give an error
	    } else {
		## possibly convert input *vector*
		if(!is.numeric(x) || is.na(n <- sizeDiss(x)))
		    stop("'x' is not and cannot be converted to class \"dissimilarity\"")
		attr(x, "Size") <- n
	    }
	    class(x) <- dissiCl
	    if(is.null(attr(x,"Metric"))) attr(x, "Metric") <- "unspecified"
	}
	## adapt S dissimilarities to Fortran:
	## convert upper matrix, read by rows, to lower matrix, read by rows.
	n <- attr(x, "Size")
	if(n > nMax)
	    stop(gettextf("have %d observations, but not more than %d are allowed",
			  n, nMax))
	dv <- x[lower.to.upper.tri.inds(n)]
	## prepare arguments for the Fortran call
	dv <- c(0, dv) ## <- internally needed {FIXME! memory hog!}
	storage.mode(dv) <- "double"
	jp <- 1
	mdata <- FALSE
	ndyst <- 0
    }
    else {
	## check input matrix and standardize, if necessary
	x <- data.matrix(x)# dropping "automatic rownames" compatibly with daisy()
        if(!(is.numeric(x) || is.logical(x))) stop("x is not a numeric dataframe or matrix.")
	x2 <- x ; dimnames(x2) <- NULL
	n <- nrow(x2)
	if(n > nMax)
	    stop(gettextf("have %d observations, but not more than %d are allowed",
			  n, nMax))
	if(stand) x2 <- scale(x2, scale = apply(x2, 2, meanabsdev))
	## put info about metric, size and NAs in arguments for the Fortran call
	metric <- match.arg(metric)
	ndyst <- c("euclidean" = 1L, "manhattan" = 2L)[[metric]]
	jp <- ncol(x2)
	if((mdata <- any(inax <- is.na(x2)))) { # TRUE if x[] has any NAs
	    jtmd <- integer(jp)
	    jtmd[apply(inax, 2L, any)] <- -1L
	    ## VALue for MISsing DATa
	    valmisdat <- 1.1* max(abs(range(x2, na.rm=TRUE)))
	    x2[inax] <- valmisdat
	}
        storage.mode(x2) <- "double"
    }
    if((k <- as.integer(k)) < 1 || k >= n)
        stop("Number of clusters 'k' must be in {1,2, .., n-1}; hence n >= 2")
    missVari <- missing(variant)
    variant <- match.arg(variant) # incl. validity check
    if(!missVari) {
        if(!missing(pamonce))
            stop("Set either 'variant' or 'pamonce', but not both")
        pamonce <- -1L +   ##  0            1      2      3      4      5       6
            match(variant, c("original", "o_1", "o_2", "f_3", "f_4", "f_5", "faster"))
        if(missing(medoids) && variant == "faster")
            medoids <- "random"
    } ## else if(!missing(pamonce)) Deprecated("use 'variant' instead")

    if(randIni <- identical("random", medoids))
        medoids <- sample.int(n, k)
    else if(!is.null(medoids)) { # non-default: check provided medoids
        ## 'fixme': consider  sort(medoids) {and rely on it in ../src/pam.c }
        if(!is.integer(medoids))
            medoids <- as.integer(medoids)
	if(length(medoids) != k || any(medoids < 1L) || any(medoids > n) ||
           any(duplicated(medoids)))
	    stop(gettextf(
		"'medoids' must be NULL or vector of %d distinct indices in {1,2, .., n}, n=%d",
		k, n))
        ## use observation numbers  'medoids' as starting medoids for 'swap' only
    }
    nisol <- integer(if(cluster.only) 1 else k)
    if(do.swap) nisol[1] <- 1L

    pamDo <- function(medoids) {
        .Call(cl_Pam, k, n,
                 !diss, # == do_diss: compute d[i,j] them from x2[] and allocate in C
                 if(diss) dv else x2,
                 !cluster.only, ## == all_stats == "old"  obj[1+ 0] == 0
                 medoids,
                 do.swap, trace.lev, keep.diss, pamonce,
                 ## only needed if(!diss) [ <=> if(do_diss) ] :
                 if(mdata) rep(valmisdat, jp) else double(1), # valmd
                 if(mdata) jtmd else integer(jp),	      # jtmd
                 ndyst)	                                      # dist_kind
    }

    res <- pamDo(medoids)
    ## Error if have NA's in diss:
    if(!diss && is.integer(res))
        stop("No clustering performed, NAs in the computed dissimilarity matrix.")
    if(randIni && nstart >= 2) {
        it <- 0L
        for(it in 2:nstart) {
            r <- pamDo(medoids = sample.int(n, k))
            if(r$obj[2] < res$obj[2]) {
                if(trace.lev)
                    cat(sprintf("Found better objective, %g < %g (it=%d)\n",
                                r$obj[2], res$obj[2], it))
                res <- r
            }
        }
    } ## else just once

    xLab <- if(diss) attr(x, "Labels") else dimnames(x)[[1]]
    r.clu <- res$clu
    if(length(xLab) > 0)
	names(r.clu) <- xLab

    if(cluster.only)
	return(r.clu)

    ## Else, usually
    medID <- res$med
    if(any(medID <= 0))
	stop("error from .C(cl_pam, *): invalid medID's")
    sildim <- res$silinf[, 4]
    if(diss) {
	## add labels to Fortran output
	r.med <- if(length(xLab) > 0) {
	    sildim <- xLab[sildim]
	    xLab[medID]
	} else medID
    }
    else {
	if(keep.diss) {
	    ## adapt Fortran output to S:
	    ## convert lower matrix, read by rows, to upper matrix, read by rows.
	    disv <- res$dys[-1]
	    disv[disv == -1] <- NA
	    disv <- disv[upper.to.lower.tri.inds(n)]
	    class(disv) <- dissiCl
	    attr(disv, "Size") <- nrow(x)
	    attr(disv, "Metric") <- metric
	    attr(disv, "Labels") <- dimnames(x)[[1]]
	}
	## add labels to Fortran output
	r.med <- x[medID, , drop=FALSE]
	if(length(xLab) > 0)
	    sildim <- xLab[sildim]
    }
    ## add names & dimnames to Fortran output
    r.obj <- structure(res$obj, names = c("build", "swap"))
    r.isol <- factor(res$isol, levels = 0:2, labels = c("no", "L", "L*"))
    names(r.isol) <- 1:k
    r.clusinf <- res$clusinf
    dimnames(r.clusinf) <- list(NULL, c("size", "max_diss", "av_diss",
					"diameter", "separation"))
    ## construct S object
    r <-
	list(medoids = r.med, id.med = medID, clustering = r.clu,
	     objective = r.obj, isolation = r.isol,
	     clusinfo = r.clusinf,
	     silinfo = if(k != 1) {
		 silinf <- res$silinf[, -4, drop=FALSE]
		 dimnames(silinf) <-
		     list(sildim, c("cluster", "neighbor", "sil_width"))
		 list(widths = silinf,
		      clus.avg.widths = res$avsil[1:k],
		      avg.width = res$ttsil)
	     },
	     diss = if(keep.diss) { if(diss) x else disv },
	     call = match.call())
    if(keep.data) { ## have !diss
	if(mdata) x2[x2 == valmisdat] <- NA
	r$data <- structure(x2, dimnames = dimnames(x))
    }
    class(r) <- c("pam", "partition")
    r
}

### From Schubert, Dec 2020 --- but MM decides to rather implement  pam(*,  variant = "faster")
if(FALSE) ## FasterPAM : Faster Partitioning Around Medoids
fasterpam <- function(x, k, diss = inherits(x, "dist"),
		metric = c("euclidean", "manhattan"), ## FIXME: add "jaccard"
                medoids = NULL,
                stand = FALSE, cluster.only = FALSE, # do.swap = TRUE, ## (not here)
                keep.diss = !diss && !cluster.only && n < 100,
                keep.data = !diss && !cluster.only,
                ## pamonce = FALSE, ## (not here)
		trace.lev = 0)
{
	if((diss <- as.logical(diss))) {
		n <- attr(x, "Size")
	} else {
		n <- nrow(x)
	}
	if (is.null(medoids)) {
		medoids = sample.int(n, k)
	}
	pam(x = x, k = k, diss = diss, metric = metric, medoids = medoids,
            stand = stand, cluster.only = cluster.only, do.swap = TRUE,
            keep.diss = keep.diss, keep.data = keep.data, pamonce = 6, trace.lev = trace.lev)
}



## non-exported:
.print.pam <- function(x, ...) {
    cat("Medoids:\n");		print(cbind(ID = x$id.med, x$medoids), ...)
    cat("Clustering vector:\n");	print(x$clustering, ...)
    cat("Objective function:\n");	print(x$objective, ...)
}

print.pam <- function(x, ...)
{
    .print.pam(x, ...)
    cat("\nAvailable components:\n")
    print(names(x), ...)
    invisible(x)
}

summary.pam <- function(object, ...)
{
    class(object) <- "summary.pam"
    object
}

print.summary.pam <- function(x, ...)
{
    .print.pam(x, ...)
    cat("\nNumerical information per cluster:\n"); print(x$clusinfo, ...)
    cat("\nIsolated clusters:\n L-clusters: ")
    print(names(x$isolation[x$isolation == "L"]), quote = FALSE, ...)
    cat(" L*-clusters: ")
    print(names(x$isolation[x$isolation == "L*"]), quote = FALSE, ...)
    if(length(x$silinfo) != 0) {
	cat("\nSilhouette plot information:\n")
	print(x$silinfo[[1]], ...)
	cat("Average silhouette width per cluster:\n")
	print(x$silinfo[[2]], ...)
	cat("Average silhouette width of total data set:\n")
	print(x$silinfo[[3]], ...)
    }
    if(!is.null(x$diss)) { ## Dissimilarities:
	cat("\n");			print(summary(x$diss, ...))
    }
    cat("\nAvailable components:\n");	print(names(x), ...)
    invisible(x)
}

