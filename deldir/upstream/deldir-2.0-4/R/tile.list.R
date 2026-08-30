tile.list <- local({
edgeLengths <- function(x,y) {
    n <- length(x)
    el <- numeric(n) 
    for(i in 1:n) {
        ii <- if(i < n) i+1 else 1
        el[i] <- sqrt((x[i] - x[ii])^2 + (y[i] - y[ii])^2)
    } 
    el 
}

function (object,minEdgeLength=NULL,clipp=NULL) {
  if (!inherits(object, "deldir"))
    stop("Argument \"object\" is not of class \"deldir\".\n")
  rw <- object$rw
  if(is.null(minEdgeLength)) {
    drw <- sqrt((rw[2] - rw[1])^2 + (rw[4] - rw[3])^2)
    minEdgeLength <- drw*sqrt(.Machine$double.eps)
  }
  x.crnrs <- rw[c(1, 2, 2, 1)]
  y.crnrs <- rw[c(3, 3, 4, 4)]
  ddd <- object$dirsgs
  sss <- object$summary
  npts <- nrow(sss)
  x <- sss[["x"]]
  y <- sss[["y"]]
  z <- sss[["z"]]
  id <- sss[["id"]]
  noid <- is.null(id)
  if(noid) id <- 1:nrow(sss)
  noz <- is.null(z)
  i.crnr <- get.cnrind(x, y, rw)
  ind.orig <- object$ind.orig
  rslt <- vector("list",npts)
  for (i in 1:npts) {
    filter1 <- ddd$ind1 == id[i]
    filter2 <- ddd$ind2 == id[i]
    subset  <- ddd[which(filter1 | filter2),,drop=FALSE]
    m <- matrix(unlist(subset[, 1:4]), ncol = 4)
    bp1 <- subset[, 7]
    bp2 <- subset[, 8]
    m1 <- cbind(m[, 1:2, drop = FALSE], 0 + bp1)
    m2 <- cbind(m[, 3:4, drop = FALSE], 0 + bp2)
    m <- rbind(m1, m2)
    pt <- c(x = sss$x[i], y = sss$y[i])
    theta <- atan2(m[, 2] - pt[2], m[, 1] - pt[1])
    theta.0 <- sort(unique(theta))
    mm <- m[match(theta.0, theta),,drop=FALSE]
    xx <- mm[, 1]
    yy <- mm[, 2]
    bp <- as.logical(mm[, 3])
    ii <- i.crnr %in% i
    xx <- c(xx, x.crnrs[ii])
    yy <- c(yy, y.crnrs[ii])
    bp <- c(bp, rep(TRUE, sum(ii)))
    tmp <- list(ptNum = ind.orig[i],
                pt    = pt,
                x     = unname(xx),
                y     = unname(yy),
                bp    = bp,
                area  = sss$dir.area[i])
    tmp    <- acw(tmp)
    bird   <- edgeLengths(tmp$x,tmp$y)
    ok     <- bird >= minEdgeLength
    tmp$x  <- tmp$x[ok]
    tmp$y  <- tmp$y[ok]
    tmp$bp <- tmp$bp[ok]
    rslt[[i]] <-acw(tmp)
    if(!noz) {
        rslt[[i]]["z"] <- z[i]
    }
    if(is.null(clipp)) {
        attr(rslt[[i]],"ncomp") <- 1
    }
}
if(!is.null(clipp)) {
    if(requireNamespace("polyclip",quietly=TRUE)) {
        rslt <- lapply(rslt,doClip,clipp=clipp,rw=rw)
    } else {
        stop("Cannot clip the tiles; package \"polyclip\" not available.\n")
    }
}
    ok <- !sapply(rslt,is.null)
    rslt <- rslt[ok]
    if(noid) {
        id <- paste0("pt.",id)
    }
    names(rslt) <- id[ok]
    class(rslt) <- "tile.list"
    attr(rslt, "rw") <- object$rw
    attr(rslt,"clipp") <- clipp
    return(rslt)
}

})

"[.tile.list" <- function(x,i,...){
    y <- unclass(x)[i]
    class(y) <- "tile.list"
    attr(y,"rw") <- attr(x,"rw")
    y
}
