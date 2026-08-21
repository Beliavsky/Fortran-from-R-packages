#' Parallel Maxima and Minima for mcnodes
#'
#' Returns the parallel maxima and minima of the input values, extended to \samp{mcnode} objects.
#' The same rules as in \code{\link{Ops.mcnode}} are applied for type determination.
#'
#' \samp{pmax} and \samp{pmin} take one or more \samp{mcnode} and/or vectors as arguments and
#' return a \samp{mcnode} of adequate type and size giving the "parallel" maxima or minima.
#' Note that the first element of \samp{...} should be an \samp{mcnode}.
#'
#' @param ... one or more \samp{mcnode}s or \samp{mcnode}s and compatible vectors.
#'   The first element must be an \samp{mcnode}.
#' @param na.rm a logical indicating whether missing values should be removed.
#' @return an \samp{mcnode} of adequate type and dimension.
#' @seealso \code{\link{min}}, \code{\link{Ops.mcnode}}
#' @keywords utilities
#' @name pmin.mcnode
#' @examples
#' ndvar(10)
#' ndunc(21)
#' x <- mcstoc(rnorm, "V")
#' pmin(x, 0)
#' y <- mcdata(rep(c(-1, 1), length=ndunc()), "U")
#' unclass(pmin(x, y))
#' @export
pmin <- function (..., na.rm = FALSE)
 UseMethod("pmin")

#' @rdname pmin.mcnode
#' @export
pmax <- function (..., na.rm = FALSE)
 UseMethod("pmax")

#' @rdname pmin.mcnode
#' @export
pmin.default <- function(..., na.rm = FALSE) base::pmin(..., na.rm = na.rm)

#' @rdname pmin.mcnode
#' @export
pmax.default <- function(..., na.rm = FALSE) base::pmax(..., na.rm = na.rm)

#' @rdname pmin.mcnode
#' @export
pmin.mcnode <- function(..., na.rm = FALSE)
{
    argsd <- list(...)

    largsd <- length(argsd)
	if(largsd == 1) return(argsd)

	outm <- attr(argsd[[1]],"outm")

    typemc <- lapply(argsd, attr, which = "type")

	type <- unlist(typemc)
	type <- if(all(type == "0")) "0" else
		if(all(type %in% c("0","V"))) "V" else
			if(all(type %in% c("0","U"))) "U" else "VU"


# evaluate the minimal common to build the minimal recycling level...
    maxdim1 <- max(unlist(lapply(argsd,function(x) dim(x)[1])))
    maxdim2 <- max(unlist(lapply(argsd,function(x) dim(x)[2])))
    maxdim3 <- max(unlist(lapply(argsd,function(x) dim(x)[3])))

#### A function to deal mcnodes as arguments, mostly copied from Ops.mcnode

    LAFUNC <- function(argsd,typethismc){

			if(!is.null(typethismc)){                                    #mcnode as arguments

				dimm <- dim(argsd)
				if ((typethismc == "V" && dimm[1] != maxdim1) ||
					(typethismc == "U" && dimm[2] != maxdim2) ||
					(typethismc == "VU" && (dimm[1] != maxdim1 || dimm[2] != maxdim2)))
					stop("Nodes of incompatible dimensions") # incompatible dimension

				if(maxdim3 > 1){  #at least one multivariate node as parameter, need recycling on the third dimension
				  if(typethismc == "U") argsd <- apply(argsd, 3, matrix, nrow=maxdim1, ncol=maxdim2, byrow=TRUE)  # recycling U as matrix (maxdim1*maxdim2) x nvariates
					else
					  argsd <- apply(argsd, 3, matrix, nrow = maxdim1, ncol = maxdim2) # recycling 0, V, VU as matrix (maxdim1*maxdim2) x nvariates
				}

				else { dim(argsd) <- NULL    # as vector
					   if(typethismc == "U" && maxdim1 != 1) argsd <- rep(argsd, each = maxdim1)     #recycling U as vector nsv*nsu
						}
			}
			else if(is.array(argsd)) stop("Array prohibited in pmin as parameter. Use an mcnode or a vector instead")
			else if(is.vector(argsd)){
				l <- length(argsd)
				if(l != 1 && l != maxdim1*maxdim2 && l!= maxdim1*maxdim2*maxdim3) stop("The vector size is not compatible
						with the node dimension. Length should be 1 or n=",maxdim1*maxdim2," or n=",maxdim1*maxdim2*maxdim3)
			}
        return(as.vector(argsd))
    }

####

	argsd <- mapply(LAFUNC, argsd=argsd, typethismc=typemc, SIMPLIFY=FALSE)

	# Copied from function pmin
	mmm <- argsd[[1L]]
    has.na <- FALSE
	for (each in argsd[-1L]) {
		if (length(mmm) < (m <- length(each)))
			mmm <- rep(mmm, length.out = m)
		else if (length(each) < (m <- length(mmm)))
			each <- rep(each, length.out = m)
		nas <- cbind(is.na(mmm), is.na(each))
		if (has.na || (has.na <- any(nas))) {
			mmm[nas[, 1L]] <- each[nas[, 1L]]
			each[nas[, 2L]] <- mmm[nas[, 2L]]
		}
		change <- mmm > each
		change <- change & !is.na(change)
		mmm[change] <- each[change]
		if (has.na && !na.rm)
			mmm[nas[, 1L] | nas[, 2L]] <- NA
	}

    mmm <- array(mmm,dim=c(maxdim1, maxdim2, maxdim3))
	class(mmm) <- "mcnode"
    attr(mmm,"type") <- type
    attr(mmm,"outm") <- outm
    return(mmm)
    }

#' @rdname pmin.mcnode
#' @export
pmax.mcnode <- function(..., na.rm = FALSE)
{
    argsd <- list(...)

    largsd <- length(argsd)
	if(largsd == 1) return(argsd)

	outm <- attr(argsd[[1]],"outm")

    typemc <- lapply(argsd, attr, which = "type")

	type <- unlist(typemc)
	type <- if(all(type == "0")) "0" else
		if(all(type %in% c("0","V"))) "V" else
			if(all(type %in% c("0","U"))) "U" else "VU"


# evaluate the minimal common to build the minimal recycling level...
    maxdim1 <- max(unlist(lapply(argsd,function(x) dim(x)[1])))
    maxdim2 <- max(unlist(lapply(argsd,function(x) dim(x)[2])))
    maxdim3 <- max(unlist(lapply(argsd,function(x) dim(x)[3])))

#### A function to deal mcnodes as arguments, mostly copied from Ops.mcnode

    LAFUNC <- function(argsd,typethismc){

			if(!is.null(typethismc)){                                    #mcnode as arguments

				dimm <- dim(argsd)
				if ((typethismc == "V" && dimm[1] != maxdim1) ||
					(typethismc == "U" && dimm[2] != maxdim2) ||
					(typethismc == "VU" && (dimm[1] != maxdim1 || dimm[2] != maxdim2)))
					stop("Nodes of incompatible dimensions") # incompatible dimension

				if(maxdim3 > 1){  #at least one multivariate node as parameter, need recycling on the third dimension
				  if(typethismc == "U") argsd <- apply(argsd, 3, matrix, nrow=maxdim1, ncol=maxdim2, byrow=TRUE)  # recycling U as matrix (maxdim1*maxdim2) x nvariates
					else
					  argsd <- apply(argsd, 3, matrix, nrow = maxdim1, ncol = maxdim2) # recycling 0, V, VU as matrix (maxdim1*maxdim2) x nvariates
				}

				else { dim(argsd) <- NULL    # as vector
					   if(typethismc == "U" && maxdim1 != 1) argsd <- rep(argsd, each = maxdim1)     #recycling U as vector nsv*nsu
						}
			}
			else if(is.array(argsd)) stop("Array prohibited in pmin as parameter. Use an mcnode or a vector instead")
			else if(is.vector(argsd)){
				l <- length(argsd)
				if(l != 1 && l != maxdim1*maxdim2 && l!= maxdim1*maxdim2*maxdim3) stop("The vector size is not compatible
						with the node dimension. Length should be 1 or n=",maxdim1*maxdim2," or n=",maxdim1*maxdim2*maxdim3)
			}
        return(as.vector(argsd))
    }

####

	argsd <- mapply(LAFUNC, argsd=argsd, typethismc=typemc, SIMPLIFY=FALSE)

	# Copied from function pmin
	mmm <- argsd[[1L]]
    has.na <- FALSE
	for (each in argsd[-1L]) {
		if (length(mmm) < (m <- length(each)))
			mmm <- rep(mmm, length.out = m)
		else if (length(each) < (m <- length(mmm)))
			each <- rep(each, length.out = m)
		nas <- cbind(is.na(mmm), is.na(each))
		if (has.na || (has.na <- any(nas))) {
			mmm[nas[, 1L]] <- each[nas[, 1L]]
			each[nas[, 2L]] <- mmm[nas[, 2L]]
		}
		change <- mmm < each
		change <- change & !is.na(change)
		mmm[change] <- each[change]
		if (has.na && !na.rm)
			mmm[nas[, 1L] | nas[, 2L]] <- NA
	}

    mmm <- array(mmm,dim=c(maxdim1, maxdim2, maxdim3))
	class(mmm) <- "mcnode"
    attr(mmm,"type") <- type
    attr(mmm,"outm") <- outm
    return(mmm)
    }
