
## This helps to define single-dispatch methods for  asNumeric() :

setOldClass("bigz")#, prototype=as.bigz(integer()))
setOldClass("bigq")#, prototype=as.bigq(integer()))
##                cannot use as.bigz() yet which is only defined in ./bigz.R

## diff() method for these: this is just base::diff.default()
## ---- with 2 lines commented out: '##>>'
.diff.big <- function(x, lag = 1L, differences = 1L, ...)
{
    ismat <- is.matrix(x)
    xlen <- if(ismat) dim(x)[1L] else length(x)
    if (length(lag) != 1L || length(differences) > 1L ||
        lag < 1L || differences < 1L)
	stop("'lag' and 'differences' must be integers >= 1")
    if (lag * differences >= xlen)
	return(x[0L]) # empty, but of proper mode
    ##>>  r <- unclass(x)  # don't want class-specific subset methods
    i1 <- -seq_len(lag)
    if (ismat)
	for (i in seq_len(differences))
	    x <- x[i1, , drop = FALSE] -
                x[-nrow(x):-(nrow(x)-lag+1L), , drop = FALSE]
    else
        for (i in seq_len(differences))
            x <- x[i1] - x[-length(x):-(length(x)-lag+1L)]
    ##>>  class(r) <- oldClass(x)
    x
}
##--> and entries in ../NAMESPACE
