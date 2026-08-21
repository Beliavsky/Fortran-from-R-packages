#' Monte Carlo Object
#'
#' Creates \samp{mc} objects from \code{\link{mcnode}} or \samp{mc} objects.
#'
#' A \samp{mc} object is a list of \code{\link{mcnode}} objects.
#' \samp{mcnode} objects must be of coherent dimensions.
#'
#' If one of the arguments is a \samp{mc} object, the names of its elements are used.
#' \samp{devname = TRUE} develops the name using the name of the \samp{mc} object as prefix.
#' Finally, names are made unique.
#'
#' @param ... \samp{mcnode} and/or \samp{mc} object(s) to be gathered in a \samp{mc} object,
#'   separated by a comma.
#' @param name vector of characters of the same length as the final \samp{mc} object.
#'   If \samp{NULL}, names are taken from the names of the elements.
#' @param devname develop the name from the name of the \samp{mc} objects, if any.
#' @return An object of class \samp{mc}.
#' @seealso \code{\link{mcnode}}, the basic element of a \samp{mc} object.
#'   To evaluate \samp{mc} objects: \code{\link{mcmodel}}, \code{\link{evalmcmod}},
#'   \code{\link{evalmccut}}.
#'   To study \samp{mc} objects: \code{\link{print.mc}}, \code{\link{summary.mc}},
#'   \code{\link{plot.mc}}, \code{\link{converg}}, \code{\link{hist.mc}},
#'   \code{\link{tornado}}, \code{\link{tornadounc}}.
#' @keywords methods
#' @examples
#' x <- mcstoc(runif)
#' y <- mcdata(3, type = "0")
#' z <- x * y
#' (m <- mc(x, y, z, name = c('n1', 'n2', 'n3')))
#' mc(m, x, devname = TRUE)
#' @export
mc <- function(..., name=NULL, devname=FALSE)
{
# the function list.names is taken from the base table function

	list.names <- function(...) {
        	l <- as.list(substitute(list(...)))[-1]
        	nm <- names(l)
        	fixup <- if (is.null(nm))
            		seq(along = l)
        	else nm == ""
        	dep <- sapply(l[fixup], function(x) if (is.symbol(x)) as.character(x) else "")
        	if (is.null(nm))
            		return(dep)
        	else {
            		nm[fixup] <- dep
            		return(nm)
        	}
    	}

#
# the function make.unique2 makes name unique

	args <- list(...)
  nameori <- names(args)
  nameobj <- unlist(list.names(...))

	if(any(duplicated(nameobj))) nameobj <- make.unique(nameobj,sep="")

  rv <- vector(mode="list",length=0)
  nom <- character(0)

  if(!all(sapply(args,inherits,"mcnode")|sapply(args,inherits,"mc"))) stop("arguments should be mc or mcnode objects")

  for(i in 1:length(args)){   # Should find better                                                  # loop to help memory

	if(is.list(args[[i]])){
      rv <- c(rv,args[[i]])
      if(devname) nom <- c(nom,paste(nameobj[i],names(args[[i]]),sep="."))
      else nom <- c(nom,names(args[[i]]))}
    else {rv <- c(rv,list(args[[i]]))
          nom <- c(nom,nameobj[i])}
    }
  rm(args)

  dimm <- sapply(rv,dim)
	if(!all(dimm[1,] %in% c(1,max(dimm[1,])))) stop("element should be of consistent variability dimensions")
	if(!all(dimm[2,] %in% c(1,max(dimm[2,])))) stop("element should be of consistent uncertainty dimensions")

  # Build the object

  typen <- sapply(rv,attr,which="type")
  type <- ifelse(all(typen < 2),"1D","2D")

  # Build the attributes

	if(!is.null(name)) {
		if(length(name)!=length(rv)) stop("the vector name is not equal to the number of rv")
		names(rv) <- make.unique(name,sep="") }
  else names(rv) <- nom

  class(rv) <- "mc"
  attr(rv,which="type") <- type
	return(rv)
}
