cvt <- function(object,stopcrit=c("change","maxit"),tol=NULL,
                maxit=100,verbose=FALSE) {
#
# Centroidal Voronoi Tessellation (by Lloyd's algorithm).
#
if(inherits(object,"deldir")) {
    l <- tile.list(object)
    rw <- object$rw
} else if(inherits(object,"tile.list")) {
    l <- object
    rw <- attr(object,"rw")
} else {
    whinge <- paste0("Argument \"object\" must be of class either\n",
                     "  \"deldir\" or \"tile.list\".\n")
    stop(whinge)
}
stopcrit <- match.arg(stopcrit)
if(stopcrit=="change") {
    if(is.null(tol)) tol <- sqrt(.Machine$double.eps)
}
g <- tile.centroids(l)
K <- 0
repeat {
    K <- K+1
    pts <- lapply(l,function(x){x$pt})
    pts <- as.data.frame(matrix(unlist(pts),byrow=TRUE,ncol=2))
    names(pts) <- c("x","y")
    dv2 <- (pts$x - g$x)^2 + (pts$y - g$y)^2
    dm  <- sqrt(max(dv2))
    if(verbose & K%%10 == 0) {
        cat("iteration:",K,"change:",dm,"\n")
    }
    if(stopcrit=="change") {
        if(dm < tol) break
    } else {
        if(K >= maxit) break
    }
    d <- try(deldir(g,rw=rw,round=FALSE))
    if(inherits(d,"try-error")) browser()
    l <- tile.list(d)
    g <- tile.centroids(l)
}
if(verbose & K%%10 != 0) cat("\n")
list(centroids=g,tiles=l)
}
