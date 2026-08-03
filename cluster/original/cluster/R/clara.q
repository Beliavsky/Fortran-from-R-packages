#### CLARA := Clustering LARge Applications
####
#### Note that the algorithm is O(n), but O(ns^2) where ns == sampsize

clara <- function(x, k,
		  metric = c("euclidean", "manhattan", "jaccard"), # "gower"),
                  stand = FALSE, cluster.only = FALSE,
		  samples = 5, sampsize = min(n, 40 + 2 * k), trace = 0,
                  medoids.x = TRUE, keep.data = medoids.x, rngR = FALSE,
                  pamLike = FALSE, correct.d = TRUE)
{
    ## check type of input matrix and values of input numbers
    if(inherits(x, "dist"))# catch user error
	stop("'x' is a \"dist\" object, but should be a data matrix or frame")
    ## for "gower" we should allow 'type'  as in daisy() --> use that code in a *modular* way
    x <- data.matrix(x)
    if(!(is.numeric(x) || is.logical(x))) stop("x is not a numeric dataframe or matrix.")
    n <- nrow(x)
    if((k <- as.integer(k)) < 1 || k > n - 1)
	stop("The number of cluster should be at least 1 and at most n-1." ) # ==> n >= 2
    if((sampsize <- as.integer(sampsize)) < max(2,k+1))
	stop(gettextf("'sampsize' should be at least %d = max(2, 1+ number of clusters)",
                      max(2,k+1)), domain=NA)
    if(n < sampsize)
	stop(gettextf("'sampsize' = %d should not be larger than the number of objects, %d",
                      sampsize, n), domain=NA)
    if((samples <- as.integer(samples)) < 1)
	stop("'samples' should be at least 1")

    jp <- ncol(x)
    namx <- dimnames(x)[[1]]
    ## standardize, if necessary {careful not to copy unnecessarily}:
    if(medoids.x) ## need to save original 'x'
        ox <- x
    else if(keep.data)
        stop("when 'medoids.x' is FALSE, 'keep.data' must be too")
    metric <- match.arg(metric)
    if(stand)
        x <- scale(x, scale = apply(x, 2, meanabsdev))
    if(keep.data)
        data <- x
    ## put info about metric, size and NAs in arguments for the .C call

    dFlag <- -1L # not used (in C code)
    if((mdata <- any(inax <- is.na(x)))) { # TRUE if x[] has any NAs
	jtmd <- integer(jp)
	jtmd[apply(inax, 2L, any)] <- -1L
	## VALue for MISsing DATa
	valmisdat <- 1.1* max(abs(range(x, na.rm=TRUE)))
	x[inax] <- valmisdat
	if(missing(correct.d))
	    warning("Distance computations with NAs: using correct instead of pre-2016 wrong formula.
Use  'correct.d=FALSE'  to get previous results or set 'correct.d=TRUE' explicitly
to suppress this warning.")
	else if(!is.finite(dFlag <- as.integer(correct.d)))
	    stop("invalid 'correct.d'")
    } else rm(inax) # save space

    res <- .C(cl_clara,
	      n,
	      jp,
	      k, 						## 3
	      clu = as.double(x),
	      samples,			# = nran
	      sampsize, 		# = nsam		## 6
	      dis   = double(1 + (sampsize * (sampsize - 1))/2),
	      as.integer(mdata),	# = mdata
	      valmd = if(mdata) rep(valmisdat, jp) else -1.,	## 9
	      jtmd  = if(mdata) jtmd else integer(1),
	      c("euclidean" = 1L, "manhattan" = 2L, "jaccard" = 3L, "gower" = 6L)[[metric]],
					# =  diss_kind (DISS_KIND : ../src/cluster.h)
	      as.logical(rngR[1]), 	# = rng_R		## 12
	      as.logical(pamLike[1]),	# = pam_like
	      as.integer(dFlag),	# = d_flag
	      integer(sampsize),	# = nrepr		## 15
	      integer(sampsize),	# = nsel
	      sample= integer(sampsize),# = nbest
	      integer(k),		# = 			## 18
	      imed = integer(k),	# = nrx
	      double(k),		# = radus
	      double(k),		# = ttd 		## 21
	      double(k),		# = ratt
	      avdis  = double(k),	# = ttbes
	      maxdis = double(k),	# = rdbes 		## 24
	      ratdis = double(k),	# = rabes
	      size  = integer(k),	# = mtt
	      obj   = double(1), 				## 27
	      avsil = double(k),
	      ttsil = double(1),
	      silinf = matrix(0, sampsize, 4), 			## 30
	      jstop = integer(1),
	      as.integer(trace),	# = trace_lev
	      double (3 * sampsize),	# = tmp			## 33
	      integer(6 * sampsize)	# = itmp
	      )[if(cluster.only) c("clu", "jstop") else substitute()] # empty index: <all>
    ## res[] components really used below:
    ## jstop, clu, silinf, dis, sample, med, imed, obj, size, maxis, avdis, ratdis,
    ## avsil, ttsil
    ## give a warning when errors occured:
    if(res$jstop) {
	if(mdata && any(aNA <- apply(inax,1, all))) {
	    i <- which(aNA)
	    nNA <- length(i)
	    pasteC <- function(...) paste(..., collapse= ",")
	    if(nNA < 13)
		stop(sprintf(ngettext(nNA,
			      "Observation %s has *only* NAs --> omit it for clustering",
			      "Observations %s have *only* NAs --> omit them for clustering!"),
			     pasteC(i)), domain = NA)
	    else
		stop(sprintf(ngettext(nNA,
			      "%d observation (%s) has *only* NAs --> omit them for clustering!",
			      "%d observations (%s ...) have *only* NAs --> omit them for clustering!"),
			     nNA, pasteC(i[1:12])), domain = NA)
	} ## else
	if(res$jstop == 1)
	    stop("Each of the random samples contains objects between which no distance can be computed.")
	if(res$jstop == 2)
	    stop(gettextf("For each of the %d samples, at least one object was found which could not be assigned to a cluster (because of missing values).", samples))
	## else {cannot happen}
	stop("invalid 'jstop' from .C(cl_clara,.): ", res$jstop)
    }
    ## 'res$clu' is potentially too large; cut down ASAP
    clustering <- as.integer(res$clu[1:n])
    if(!is.null(namx)) names(clustering) <- namx
    if(cluster.only)
        return(clustering)
    ## adapt C output to S:
    ## convert lower matrix, read by rows, to upper matrix, read by rows.
    disv <- res$dis[-1]
    disv[disv == -1] <- NA
    disv <- disv[upper.to.lower.tri.inds(sampsize)]
    class(disv) <- dissiCl
    attr(disv, "Size") <- sampsize
    attr(disv, "Metric") <- metric
    attr(disv, "Labels") <- namx[res$sample]
    res$med <- if(medoids.x) ox[res$imed, , drop = FALSE]
    ## add labels to C output
    sildim <- res$silinf[, 4]
    if(!is.null(namx)) {
	sildim     <- namx[sildim]
	res$sample <- namx[res$sample]
    }
    r <- list(sample = res$sample, medoids = res$med, i.med = res$imed,
	      clustering = clustering, objective = res$obj,
	      clusinfo = cbind(size = res$size, "max_diss" = res$maxdis,
	      "av_diss" = res$avdis, isolation = res$ratdis),
	      diss = disv, call = match.call())
    ## add dimnames to C output
    if(k > 1) {
	dimnames(res$silinf) <- list(sildim,
				     c("cluster", "neighbor", "sil_width", ""))
	r$silinfo <- list(widths          = res$silinf[, -4],
			  clus.avg.widths = res$avsil,
			       avg.width  = res$ttsil)
    }
    if(keep.data) r$data <- data
    class(r) <- c("clara", "partition")
    r
}

print.clara <- function(x, ...)
{
    cat("Call:	", deparse1(x$call),
	"\nMedoids:\n");		print(x$medoids, ...)
    cat("Objective function:\t ", format(x$objective, ...),"\n",
	"Clustering vector: \t", sep=""); str(x$clustering, vec.len = 7)
    cat("Cluster sizes:	    \t", x$clusinfo[,1],
	"\nBest sample:\n");		print(x$sample, quote = FALSE, ...)
    cat("\nAvailable components:\n");	print(names(x), ...)
    invisible(x)
}

summary.clara <- function(object, ...)
{
    class(object) <- "summary.clara"
    object
}

print.summary.clara <- function(x, ...)
{
    cat("Object of class 'clara' from call:\n", deparse1(x$call),
	"\nMedoids:\n");		print(x$medoids, ...)
    cat("Objective function:\t ", format(x$objective, ...),
	"\nNumerical information per cluster:\n")
    print(x$clusinfo, ...)
    if(has.sil <- !is.null(x$silinfo)) {
	cat("Average silhouette width per cluster:\n")
	print(x$silinfo[[2]], ...)
	cat("Average silhouette width of best sample:",
	    format(x$silinfo[[3]], ...), "\n")
    }
    cat("\nBest sample:\n");		print(x$sample, quote = FALSE, ...)
    cat("Clustering vector:\n");	print(x$clustering, ...)
    if(has.sil) {
	cat("\nSilhouette plot information for best sample:\n")
	print(x$silinfo[[1]], ...)
    }
    if(!is.null(x$diss)) { ## Dissimilarities:
	cat("\n");			print(summary(x$diss, ...))
    }
    cat("\nAvailable components:\n");	print(names(x), ...)
    invisible(x)
}

