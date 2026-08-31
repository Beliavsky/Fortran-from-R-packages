# HEADER ####################################################
# This is file spam/R/tailhead.R.                           #
# It is part of the R package spam,                         #
#  --> https://CRAN.R-project.org/package=spam              #
#  --> https://CRAN.R-project.org/package=spam64            #
#  --> https://git.math.uzh.ch/reinhard.furrer/spam         #
# by Reinhard Furrer [aut, cre], Florian Gerber [aut],      #
#    Roman Flury [aut], Daniel Gerber [ctb],                #
#    Kaspar Moesinger [ctb]                                 #
# HEADER END ################################################



head.spam <- function(x, n = 6L, m = n, ...) {
  stopifnot(length(n) == 1L, length(m) == 1L)
   n <- if (n < 0L)
        max(dim(x)[1] + n, 0L)
    else min(n, dim(x)[1])
   m <- if (m < 0L)
        max(dim(x)[2] + m, 0L)
    else min(m, dim(x)[2])
    as.matrix(x[seq_len(n), seq_len(m), drop = FALSE])
}

tail.spam <- function (x, n = 6L, m = n, addrownums = TRUE,  ...)
{
    stopifnot(length(n) == 1L, length(m) == 1L)
    nrx <- dim(x)[1]
    ncx <- dim(x)[2]
    n <- if (n < 0L)
        max(nrx + n, 0L)
    else min(n, nrx)
    m <- if (m < 0L)
        max(ncx + m, 0L)
    else min(m, ncx)
    selr <- seq.int(to = nrx, length.out = n)
    selc <- seq.int(to = ncx, length.out = n)
    ans <- as.matrix( x[selr, selc, drop = FALSE])
    if (addrownums) {
#    if (addrownums && is.null(rownames(x))) { # rownames is null by default
#        rownames(ans) <- paste0("[", selr, ",]") # can be used from R2.15
#        colnames(ans) <- paste0("[,", selc, "]")
        rownames(ans) <- paste("[", selr, ",]", sep = "")
        colnames(ans) <- paste("[,", selc, "]", sep = "")
      }
    ans
}


setMethod("head","spam",head.spam)
setMethod("tail","spam",tail.spam)

setMethod("head","spam.chol.NgPeyton", function(x, n = 6L, m = n, ...) head.spam(as.spam(x), n = 6L, m = n, ...))
setMethod("tail","spam.chol.NgPeyton", function(x, n = 6L, m = n, addrownums = TRUE,  ...) tail.spam(as.spam(x), n = 6L, m = n, addrownums = TRUE,  ...))
