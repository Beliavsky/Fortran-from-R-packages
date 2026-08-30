plot.tile.list <- function (x, verbose = FALSE, close = FALSE, pch = 1,
                            fillcol = getCol(x,warn=warn), col.pts=NULL,
                            col.lbls=NULL,border=NULL, showpoints = !labelPts,
                            add = FALSE, asp = 1, clipp=NULL, xlab = "x",
                            ylab = "y", main = "", axes=TRUE, warn=TRUE,
                            labelPts=FALSE,adj=NULL,...) {
# Check for use of the defunct argument name "number".
ccc <- match.call()
i   <- match("number",names(ccc))
if(!is.na(i)) {
    if("labelPts" %in% names(ccc)) {
        whinge <- paste0("Both \"labelPts\" and the defunct argument",
                         " \"number\" have been\n  specified.  Do not use",
                         " the defunct argument \"number\".  Use\n",
                         "  \"labelPts\" only.\n")
        stop(whinge)
    }
    whinge <- paste0("The argument name \"number\" is defunct. Please",
                     " use \"labelPts\"\n  instead.\n")
    warning(whinge)
    names(ccc)[i] <- "labelPts"
    return(eval(ccc))
}

# Carry on.
    object <- x
    if (!inherits(object, "tile.list")) 
        stop("Argument \"object\" is not of class tile.list.\n")
    clip  <- !is.null(clipp)
    if(clip) {
        if(!is.null(attr(object,"clipp"))) {
            whinge <- paste0("Argument \"x\" is already clipped.  Re-clip it\n",
                             "  if you want a different clipping polygon.\n")
            stop(whinge)
            if(!requireNamespace("polyclip",quietly=TRUE)) {
                stop("Cannot clip the tiles; package \"polyclip\" not available.\n")
            }
        }
    }
    n     <- length(object)
    rw    <- attr(object, "rw")
    rx    <- rw[1:2]
    ry    <- rw[3:4]
    x.pts <- unlist(lapply(object, function(w) {
        w$pt[1]
    }))
    y.pts <- unlist(lapply(object, function(w) {
        w$pt[2]
    }))
    if (!add) 
        plot(0, 0, type = "n", asp = asp, xlim = rx, ylim = ry, 
            xlab = xlab, ylab = ylab, main = main, axes=axes)
    fillcol <- apply(col2rgb(fillcol, TRUE), 2, function(x) {
        do.call(rgb, as.list(x/255))
    })
    fillcol <- rep(fillcol, length = length(object))
    hexbla <- do.call(rgb, as.list(col2rgb("black", TRUE)/255))
    hexwhi <- do.call(rgb, as.list(col2rgb("white", TRUE)/255))
    if(is.null(col.pts)){
        col.pts <- ifelse(fillcol == hexbla, hexwhi, hexbla)
    } else {
        col.pts <- apply(col2rgb(col.pts, TRUE), 2, function(x) {
            do.call(rgb, as.list(x/255))
        })
        col.pts <- rep(col.pts, length = length(object))
    }
    if(is.null(col.lbls)){
        col.lbls <- ifelse(fillcol == hexbla, hexwhi, hexbla)
    } else {
        col.lbls <- apply(col2rgb(col.lbls, TRUE), 2, function(x) {
            do.call(rgb, as.list(x/255))
        })
        col.lbls <- rep(col.lbls, length = length(object))
    }
    if(is.null(border)) {
        border <- if(all(fillcol == hexbla)) hexwhi else hexbla
    } else if(length(border) > 1) border <- border[1]
    lnwid <- if(all(fillcol == hexbla)) 2 else 1
    ptNms <- names(x)
    Adj <- adj
    if(is.null(Adj)) Adj <- if(showpoints) -1 else 0
    pch <- rep(pch,n)
    pgons <- vector("list",n)
    icol <- 0
    for(i in 1:n) {
        pgon <- if(clip) doClip(object[[i]],clipp,rw) else object[[i]]
        pgons[[i]] <- pgon
        if(is.null(pgon)) next
        icol <- icol+1
        if(is.null(attr(pgon,"ncomp"))) attr(pgon,"ncomp") <- 1
        if(attr(pgon,"ncomp") > 1) {
            pgon <- pgon$tileParts
        } else pgon <- list(pgon)
        for(ii in seq(along=pgon)){
            ptmp <- pgon[[ii]]
            inner <- !any(ptmp$bp)
            polygon(ptmp,col=fillcol[icol],border=NA)
            if (close | inner) { 
                polygon(ptmp,col = NA, border = border, lwd = lnwid)
            } else {
                x <- ptmp$x
                y <- ptmp$y
                ni <- length(x)
                for (j in 1:ni) {
                    jnext <- if (j < ni) j + 1 else 1
                    do.it <- mid.in(x[c(j, jnext)], y[c(j, jnext)], rx, ry)
                    if (do.it) 
                        segments(x[j], y[j], x[jnext], y[jnext],
                                 col = border, lwd = lnwid)
                }
            }
         }
         if(verbose) {
             if(showpoints) points(object[[i]]$pt[1], object[[i]]$pt[2],
                                   pch = pch[i], col = col.pts[i],...)
             if(labelPts) text(object[[i]]$pt[1], object[[i]]$pt[2],
                             labels=ptNms[i], col = col.lbls[i],adj=Adj,...)
             if(i < n) readline(paste("i = ",i,"; Go? ",sep=""))
             if(i == n) cat("i = ",i,"\n",sep="")
         }
    }
    ok <- !sapply(pgons,is.null)
    if(showpoints & !verbose)
        points(x.pts[ok], y.pts[ok], pch = pch[ok], col = col.pts[ok],...)
    if (labelPts & !verbose) 
        text(x.pts[ok], y.pts[ok], labels = ptNms[ok], col = col.lbls[ok],
             adj=Adj,...)
    pgons <- pgons[ok]
    pgons <- if(length(pgons)) pgons else NULL
    invisible(pgons)
}
