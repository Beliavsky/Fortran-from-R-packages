#### $Id: agnes.q 8468 2024-12-10 15:11:21Z maechler $

agnes <- function(x, diss = inherits(x, "dist"), metric = "euclidean",
		  stand = FALSE, method = "average", par.method,
                  keep.diss = n < 100, keep.data = !diss, trace.lev = 0)
{
    METHODS <- c("average", "single","complete", "ward","weighted", "flexible", "gaverage")
    ## hclust has more;  1    2         3           4       5         6         7
    meth <- pmatch(method, METHODS)
    if(is.na(meth)) stop("invalid clustering method")
    if(meth == -1) stop("ambiguous clustering method")
    cl. <- match.call()
    method <- METHODS[meth]
    if(method == "flexible") {
	## Lance-Williams formula (but *constant* coefficients):
	stopifnot((np <- length(a <- as.numeric(par.method))) >= 1)
	attr(method,"par") <- par.method <-
	    if(np == 1)## default (a1= a, a2= a, b= 1-2a, c = 0)
		c(a, a, 1-2*a, 0)
	    else if(np == 3)
		c(a, 0)
	    else if(np == 4) a
	    else stop("'par.method' must be of length 1, 3, or 4")
        ## if(any(par.method[1:2]) < 0)
        ##     warning("method \"flexible\": alpha_1 or alpha_2 < 0 can give invalid dendrograms"
    } else if (method == "gaverage") {
        attr(method,"par") <- par.method <-
            if (missing(par.method)) {
                ## Default par.method: Using beta = -0.1 as advised in Belbin et al. (1992)
                beta <- -0.1
                c(1-beta, 1-beta, beta, 0)
            } else {
                stopifnot((np <- length(b <- as.numeric(par.method))) >= 1)
                if(np == 1)## default (a1= 1-b, a2= 1-b, b= b, c= 0)
		c(1-b, 1-b, b, 0)
                else if(np == 3)
                    c(b, 0)
                else if(np == 4) b
                else stop("'par.method' must be of length 1, 3, or 4")
            }
        ## if(any(par.method[1:2]) < 0)
        ##     warning("method \"gaverage\": alpha_1 or alpha_2 < 0 can give invalid dendrograms"
    } else ## dummy (passed to C; length >= 1 : using `alpha--` there
	par.method <- double(1)

    if((diss <- as.logical(diss))) {
	## check type of input vector
	if(anyNA(x)) stop("NA-values in the dissimilarity matrix not allowed.")
	if(data.class(x) != "dissimilarity") { # try to convert to
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
	n <- attr(x, "Size")
	dv <- x[lower.to.upper.tri.inds(n)] # is *slow* [c * n^2 ; but large c]  in large cases
	## prepare arguments for the Fortran call
	dv <- c(0., dv)# "double", 1st elem. "only for Fortran" (?)
	jp <- 1L
	mdata <- FALSE
	ndyst <- 0
	x2 <- double(1)
    }
    else {
	## check input matrix and standardize, if necessary
	x <- data.matrix(x)
	if(!is.numeric(x)) stop("x is not a numeric dataframe or matrix.")
	x2 <- if(stand) scale(x, scale = apply(x, 2, meanabsdev)) else x
        storage.mode(x2) <- "double"
	ndyst <- if(metric == "manhattan") 2 else 1
	n <- nrow(x2)
	jp <- ncol(x2)
	if((mdata <- any(inax <- is.na(x2)))) { # TRUE if x[] has any NAs
	    jtmd <- integer(jp)
	    jtmd[apply(inax, 2L, any)] <- -1L
	    ## VALue for MISsing DATa
	    valmisdat <- 1.1* max(abs(range(x2, na.rm=TRUE)))
	    x2[inax] <- valmisdat
	}
	dv <- double(1 + (n * (n - 1))/2)
    }
    if(n <= 1) stop("need at least 2 objects to cluster")
    stopifnot(length(trace.lev <- as.integer(trace.lev)) == 1)
    C.keep.diss <- keep.diss && !diss
    res <- .C(twins,
		    as.integer(n),
		    as.integer(jp),
		    x2,
		    dv,
		    dis = double(if(C.keep.diss) length(dv) else 1L),
		    jdyss = if(C.keep.diss) diss + 10L else as.integer(diss),
		    if(mdata && jp) rep(valmisdat, jp) else double(1L),
		    if(mdata) jtmd else integer(jp),
		    as.integer(ndyst),
		    1L,# jalg = 1 <==> AGNES
		    meth,# integer
		    integer(n),
		    ner = integer(n),
		    ban = double(n),
		    ac = double(1), # coef
                    par.method, # = alpha (of length 1, 3, or 4)
		    merge = matrix(0L, n - 1, 2), # integer
                    trace = trace.lev)[c("dis", "jdyss", "ner", "ban", "ac", "merge")]
    if(!diss) {
	##give warning if some dissimilarities are missing.
	if(res$jdyss == -1)
	    stop("No clustering performed, NA-values in the dissimilarity matrix.\n" )
        if(keep.diss) {
            ## adapt Fortran output to S:
            ## convert lower matrix,read by rows, to upper matrix, read by rows.
            disv <- res$dis[-1]
            disv[disv == -1] <- NA
            disv <- disv[upper.to.lower.tri.inds(n)]
            class(disv) <- dissiCl
            attr(disv, "Size") <- nrow(x)
            attr(disv, "Metric") <- metric
            attr(disv, "Labels") <- dimnames(x)[[1]]
        }
	##add labels to Fortran output
	if(length(dimnames(x)[[1]]) != 0)
	    order.lab <- dimnames(x)[[1]][res$ner]
    }
    else {
        if(keep.diss) disv <- x
	##add labels to Fortran output
	if(length(attr(x, "Labels")) != 0)
	    order.lab <- attr(x, "Labels")[res$ner]
    }
    clustering <- list(order = res$ner, height = res$ban[-1], ac = res$ac,
		       merge = res$merge, diss = if(keep.diss)disv,
		       call = cl., method = METHODS[meth])
    if(exists("order.lab"))
	clustering$order.lab <- order.lab
    if(keep.data && !diss) {
	if(mdata) x2[x2 == valmisdat] <- NA
	clustering$data <- x2
    }
    class(clustering) <- c("agnes", "twins")
    clustering
}

summary.agnes <- function(object, ...)
{
    class(object) <- "summary.agnes"
    object
}

print.agnes <- function(x, ...)
{
    cat("Call:	", deparse1(x$call),
	"\nAgglomerative coefficient: ", format(x$ac, ...),
	"\nOrder of objects:\n")
    print(if(length(x$order.lab) != 0) x$order.lab else x$order,
	  quote = FALSE, ...)
    cat("Height (summary):\n");		print(summary(x$height), ...)
    cat("\nAvailable components:\n");	print(names(x), ...)
    invisible(x)
}

print.summary.agnes <- function(x, ...)
{
    ## a bit more than print.agnes() ..
    cat("Object of class 'agnes' from call:\n", deparse1(x$call),
	"\nAgglomerative coefficient: ", format(x$ac, ...),
	"\nOrder of objects:\n")
    print(if(length(x$order.lab) != 0) x$order.lab else x$order,
	  quote = FALSE, ...)
    cat("Merge:\n");			print(x$merge, ...)
    cat("Height:\n");			print(x$height, ...)
    if(!is.null(x$diss)) { ## Dissimilarities:
	cat("\n");			print(summary(x$diss, ...))
    }
    cat("\nAvailable components:\n");	print(names(x), ...)
    invisible(x)
}

as.dendrogram.twins <- function(object, ...) ## ... : really only 'hang'
    as.dendrogram(as.hclust(object), ...)
