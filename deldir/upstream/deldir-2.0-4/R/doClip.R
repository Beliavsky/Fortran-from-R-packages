doClip <- function(object,clipp,rw) {
    pgon <- polyclip::polyclip(object,clipp)
    n    <- length(pgon)
    if(n) {
        rslt  <- object[c("ptNum","pt")]
        cmps  <- vector("list",n)
        xold  <- object$x
        yold  <- object$y
        bpold <- object$bp
        for(ii in 1:n) {
            xnew <- pgon[[ii]]$x
            ynew <- pgon[[ii]]$y
            imtch <- findNewInOld(xnew,xold,ynew,yold)
            bp         <- rep(FALSE,length(imtch))
            bp[imtch!=0]  <- bpold[imtch]
            area       <- tileArea(xnew,ynew,rw)
            cmps[[ii]] <- list(x=xnew,y=ynew,bp=bp,area=area)
        }
        if(n==1) {
            rslt <- c(rslt,cmps[[1]])
        } else {
            rslt <- c(rslt,list(tileParts=cmps))
        }
        attr(rslt,"ncomp") <- n
    } else rslt <- NULL
    rslt
}
