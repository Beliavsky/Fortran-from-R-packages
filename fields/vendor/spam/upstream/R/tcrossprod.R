# HEADER ####################################################
# This is file spam/R/tcrossprod.R.                         #
# It is part of the R package spam,                         #
#  --> https://CRAN.R-project.org/package=spam              #
#  --> https://CRAN.R-project.org/package=spam64            #
#  --> https://git.math.uzh.ch/reinhard.furrer/spam         #
# by Reinhard Furrer [aut, cre], Florian Gerber [aut],      #
#    Roman Flury [aut], Daniel Gerber [ctb],                #
#    Kaspar Moesinger [ctb]                                 #
# HEADER END ################################################

     

crossprod.spam <- function(x, y=NULL, ...) {
    dimx <- dim(x)
    if( is.null(y)) {
        if(!is.spam(x)) return(crossprod(x))

        if(dimx[2]==1L) return(matrix( sum(x@entries^2)))

        return( t.spam(x) %*% x)
    }
    if( (!is.spam(x)) & (!is.spam(y))) return(crossprod(x,y))
    
    return( t(x) %*% y)
    
}
tcrossprod.spam <- function(x, y=NULL, ...) {
    dimx <- dim(x)
    if( is.null(y)) {
        if(!is.spam(x)) return(tcrossprod(x))

        if(dimx[2]==1L) return(matrix( sum(x@entries^2)))

        return( x %*% t.spam(x))
    }
    if( (!is.spam(x)) & (!is.spam(y))) return(tcrossprod(x,y))
    
    return( x %*% t(y))
    
}


setMethod("crossprod",signature(x="spam",y="missing"), crossprod.spam)
setMethod("crossprod",signature(x="spam",y="spam"), crossprod.spam)
setMethod("crossprod",signature(x="spam",y="ANY"), crossprod.spam)
setMethod("crossprod",signature(x="ANY",y="spam"), crossprod.spam)
setMethod("tcrossprod",signature(x="spam",y="missing"), tcrossprod.spam)
setMethod("tcrossprod",signature(x="spam",y="spam"), tcrossprod.spam)
setMethod("tcrossprod",signature(x="spam",y="ANY"), tcrossprod.spam)
setMethod("tcrossprod",signature(x="ANY",y="spam"), tcrossprod.spam)
