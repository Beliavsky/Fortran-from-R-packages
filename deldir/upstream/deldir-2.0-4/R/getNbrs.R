getNbrs <- function(object,interior=NULL) {
    ddd  <- object$dirsgs
    id   <- object$summary[["id"]]
    noid <- is.null(id)
    if(noid) id <- 1:nrow(object$summary)
    if(is.null(interior)) {
        npts <- nrow(object$summary)
    } else {
        x <- object$summary[["x"]]
        y <- object$summary[["y"]]
        if(inherits(interior,"list")) { # "interior" is/should be a polygon
            ok <- insidePoly(x,y,interior)
        } else {
            rw   <- object[["rw"]]
            ok <- insideRect(x,y,interior,rw)
        }
        id   <- id[ok]
        npts <- length(id)
    }
    nbrs <- vector("list",npts)
    for(i in seq(along=id)) {
        filter1   <- ddd$ind1 == id[i]
        filter2   <- ddd$ind2 == id[i]
        subset    <- ddd[which(filter1 | filter2),,drop=FALSE]
        nbrs[[i]] <- unname(apply(subset[,c("ind1","ind2")],1,
                                  function(x){x[x!=id[i]]}))
    }
    if(noid) {
        names(nbrs) <- paste0("pt.",id)
    } else {
        names(nbrs) <- id
    }
    nbrs
}
