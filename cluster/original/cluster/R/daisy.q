
daisy <- function(x, metric = c("euclidean", "manhattan", "gower"),
		  stand = FALSE, type = list(), weights = rep.int(1, p),
		  warnBin = warnType, warnAsym = warnType, warnConst = warnType,
		  warnType = TRUE)
{
    ## check type of input matrix
    if(length(dx <- dim(x)) != 2 || !(is.data.frame(x) || is.numeric(x)))
	stop("x is not a dataframe or a numeric matrix.")
    n <- dx[1]# nrow
    p <- dx[2]# ncol
    varnms <- dimnames(x)[[2]]
    pColl <- function(n) paste(n, collapse = ", ")
    if(length(type)) {
	if(!is.list(type) || is.null(ntyp <- names(type)) || any(ntyp == ""))
	    stop(gettextf("invalid %s; must be named list", sQuote("type")))
	## check each component to be valid column names or numbers:
	for(nt in ntyp) {
	    cvec <- type[[nt]]
	    ct <- paste0("type$", nt)
	    if(is.character(cvec)) {
		if(!is.null(varnms) && !all(cvec %in% varnms))
		    stop(gettextf("%s has invalid column names", ct))
	    }
	    else if(is.numeric(cvec)) {
		if(!all(1 <= cvec & cvec <= p))
		    stop(gettextf("%s must be in 1:ncol(x)", ct))
	    }
	    else stop(gettextf("%s must contain column names or numbers", ct))
	}
	tA <- type$asymm
	tS <- type$symm
	if(!is.null(tA) || !is.null(tS)) { ## check binary columns:
	    ## tA and tS might be character and integer!
	    d.bin <- cbind(as.data.frame(x[, tA, drop= FALSE]),
					 x[, tS, drop= FALSE])
	    lenB <- lengths(lapply(d.bin, function(y) levels(as.factor(y))))
	    if(any(lenB > 2))
		stop("at least one binary variable has more than 2 levels.")
	    if(any(lenB < 2))
		warning("at least one binary variable has not 2 different levels.")
	    ## Convert factors to integer, such that ("0","1") --> (0,1):
	    if(any(is.f <- vapply(d.bin, is.factor, NA)))
		d.bin[is.f] <- lapply(d.bin[is.f],
				      function(f) as.integer(as.character(f)))
	    if(!all(vapply(d.bin, FUN.VALUE=NA, function(y)
			   is.logical(y) ||
			   all(sort(unique(as.numeric(y[!is.na(y)])))%in% 0:1))))
		stop("at least one binary variable has values not in {0,1,NA}")
	}
    }
    ## transform variables and construct 'type' vector
    if(is.data.frame(x)) {
	type2 <- vapply(x, data.class, "<class>")
	x <- data.matrix(x)
    } else { ## matrix
        type2 <- rep("numeric", p)
        names(type2) <- colnames(x)
    }
    if(length(type)) {
	tT <- type$ ordratio
	tL <- type$ logratio # lg-transformed, staying numeric, -> treated as "I"
	x[, names(type2[tT])] <- unclass(as.ordered(x[, names(type2[tT])]))
	x[, names(type2[tL])] <- log10(		    x[, names(type2[tL])])
	type2[tA] <- "A"
	type2[tS] <- "S"
	type2[tT] <- "T" # was "O" (till 2000-12-14) accidentally !
    }
    type2[tI <- type2 %in% c("numeric", "integer") ] <- "I"
    if(warnBin && n > 9 && any(tI) &&
       any(iBin <- apply(x[, tI, drop = FALSE], 2,
			 function(v) length(table(v)) == 2)))
	warning(gettextf("binary variable(s) %s treated as interval scaled",
			 pColl(which(tI)[iBin])))

    type2[type2 == "ordered"] <- "O"
    type2[type2 == "factor"]  <- "N"
    if(any(ilog <- type2 == "logical")) {
	if(warnAsym) warning(sprintf(ngettext(sum(ilog),
				 "setting 'logical' variable %s to type 'asymm'",
				 "setting 'logical' variables %s to type 'asymm'"),
			pColl(which(ilog))), domain = NA)
	type2[ilog] <- "A"
    }
    ## Note: We have 2 status codings:  ndyst = (0,1,2) and jdat = (1,2);
    ##       the latter is superfluous in principle

    ## standardize, if necessary, i.e., *iff* all vars are "I"nterval scaled:
    all.I <- all(type2 == "I")
    miss_metric <- missing(metric)
    metric <- match.arg(metric)
    if(all.I && metric != "gower") {
	if(stand) {
	    x <- scale(x, center = TRUE, scale = FALSE) #-> 0-means
	    sx <- colMeans(abs(x), na.rm = TRUE)# can still have NA's
	    if(0 %in% sx) {
		if(warnConst) warning(gettextf(
		    "%s has constant columns %s; these are standardized to 0",
		    sQuote("x"), pColl(which(sx == 0))))
		sx[sx == 0] <- 1
	    }
	    x <- scale(x, center = FALSE, scale = sx)
	}
	jdat <- 2L
	ndyst <- if(metric == "manhattan") 2L else 1L # == diss_kind
    }
    else { ## mixed case or explicit "gower"
	if(!miss_metric && metric != "gower" && !all.I)
	    warning("with mixed variables, metric \"gower\" is used automatically")
        ## FIXME: think of a robust alternative scaling to
        ##        Gower's  (x - min(x)) / (max(x) - min(x))
	colR <- apply(x, 2, range, na.rm = TRUE)
	colmin <- colR[1,]
	sx <- colR[2,] - colmin
	if(any(sx == 0))
	    sx[sx == 0] <- 1
	x <- scale(x, center = colmin, scale = sx)
	jdat <- 1L
	ndyst <- 0L ## diss_kind = "mixed | gower"
        ## weights only used in this "gower" case
        if(length(weights) == 1)
            weights <- rep.int(weights, p)
        else if(length(weights) != p)
            stop("'weights' must be of length p (or 1)")
    }

    ##	type2 <- paste(type2, collapse = "")
    typeCodes <- c('A','S','N','O','I','T')
    ##              1   2   3   4   5   6  --> passed to Fortran below
    iType <- match(type2, typeCodes)# integer
    if(any(ina <- is.na(iType)))
	stop(gettextf("invalid type %s for column numbers %s",
		      type2[ina], pColl(which(ina))))
    if((mdata <- any(inax <- is.na(x)))) { # TRUE if x[] has any NAs (hence p >= 1)
	jtmd <- integer(p)
	jtmd[apply(inax, 2L, any)] <- -1L
	## VALue for MISsing DATa
	valmisdat <- 1.1* max(abs(range(x, na.rm=TRUE)))
	x[inax] <- valmisdat
    }
    ## call Fortran routine
    storage.mode(x) <- "double"
    disv <- .C(cl_daisy, # -> ../src/daisy.c
		     n,
		     p,
		     x,
		     if(mdata) rep(valmisdat, p) else double(1L),
                     as.double(weights),
		     if(mdata) jtmd else integer(1L),
		     jdat,
		     iType,		# vtype
		     ndyst,
		     as.integer(mdata),
		     dis = double(if(n <= 2L) 1L else (n * (n - 1))/2), # at least 1: C has  `disv--`
		     NAOK = TRUE # only to allow "+- Inf"
		     )$dis
    ## adapt Fortran output to S:
    ## convert lower matrix, read by rows, to upper matrix, read by rows.
    if(n) {
        disv[disv == -1] <- NA
        full <- matrix(0, n, n)
        full[!lower.tri(full, diag = TRUE)] <- disv
        disv <- t(full)[lower.tri(full)]
        ## give warning if some dissimilarities are missimg
        if(anyNA(disv)) attr(disv, "NA.message") <-
                            "NA-values in the dissimilarity matrix !"
    } else # n == 0  {pmax(1,.) above}
        disv <- double(0L)
    ## construct S object -- "dist" methods are *there* !
    class(disv) <- dissiCl # see ./0aaa.R
    attr(disv, "Labels") <- dimnames(x)[[1]]
    attr(disv, "Size") <- n
    attr(disv, "Metric") <- if(!ndyst) "mixed" else metric
    if(!ndyst) attr(disv, "Types") <- typeCodes[iType]
    disv
}

print.dissimilarity <-
    function(x, diag = NULL, upper = NULL,
	     digits = getOption("digits"), justify = "none", right = TRUE, ...)
{
    cat("Dissimilarities :\n")
    NextMethod("print")##-> stats:::print.dist(..)
    cat("\n")
    if(!is.null(attr(x, "na.message")))
	cat("Warning : ", attr(x, "NA.message"), "\n")
    cat("Metric : ", attr(x, "Metric"),
	if(!is.null(aT <- attr(x,"Types")))
	paste(";  Types =", paste(aT, collapse=", ")), "\n")
    cat("Number of objects : ", attr(x, "Size"), "\n", sep="")
    invisible(x)
}

summary.dissimilarity <-
    function(object, digits = max(3, getOption("digits") - 2), ...)
    ## 'digits': want a bit higher precision
{
    sx <- summary(as.vector(object), digits = digits, ...)
    at <- attributes(object)
    r <- c(list(summ = sx, n = length(object)), at[names(at) != "class"])
    class(r) <- "summary.dissimilarity"
    r
}

print.summary.dissimilarity <- function(x, ...)
{
    cat(x$n, "dissimilarities, summarized :\n")
    print(x$summ, ...)
    cat("Metric : ", x $ Metric,
	if(!is.null(aT <- x $ Types))
	paste(";  Types =", paste(aT, collapse=", ")), "\n")
    cat("Number of objects : ", x $ Size, "\n", sep="")
    if(!is.null(x $ na.message))
	cat("Warning : ", x $ NA.message, "\n")
    invisible(x)
}
