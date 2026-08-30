tilePerim0 <- function (object,inclbdry=TRUE) {
    ncomp <- attr(object,"ncomp")
    if(is.null(ncomp)) ncomp <- 1
    tobj <- if(ncomp==1) list(object) else object$tileParts
    peri <- numeric(ncomp)
    for(i in 1:ncomp) {
        tmp <- tobj[[i]]
        x <- tmp[["x"]]
        y <- tmp[["y"]]
        xx <- c(x,x[1])
        yy <- c(y,y[1])
        if(inclbdry) {
            ok <- rep(TRUE,length(x))
        } else {
            bp1 <- tmp[["bp"]]
            bp2 <- c(bp1,bp1[1])
            bpm <- cbind(bp1,bp2[-1])
            ok  <- !apply(bpm,1,all)
        }
        peri[i] <- sum(sqrt(((xx[-1] - x)[ok])^2 + ((yy[-1] - y)[ok])^2))
    }
    perTot <- sum(peri)
    attr(perTot,"perComps") <- peri
    perTot
}
