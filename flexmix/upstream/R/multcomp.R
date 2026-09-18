#
#  Copyright (C) 2004-2016 Friedrich Leisch
#  $Id: multcomp.R 5298 2025-02-26 14:50:01Z gruen $
#
###*********************************************************

setGeneric("flxglht",
           function(model, linfct, ...) standardGeneric("flxglht"))

setMethod("flxglht", signature(model="flexmix", linfct="character"),
function(model, linfct, ...)
{
    model <- refit(model)
    flxglht(model, linfct, ...)
})

setMethod("flxglht", signature(model="FLXRoptim", linfct="character"),
function(model, linfct, ...)
{
    if (!requireNamespace("multcomp", quietly = TRUE)) {
        stop("install package multcomp to use this function")
    }
    if(length(model@components)>1)
        stop("Can currently handle only models with one response!\n")

    type <- match.arg(linfct, c("zero", "tukey"))

    k <- model@k
    cf <- model@coef
    vc <- model@vcov

    nc <- rownames(model@components[[1]][[1]])
    nc <- as.vector(outer(paste("model.1_Comp.", 1:k, "_coef.", sep=""),
                          nc, paste, sep=""))

    ## FIXME: for zero-inflated models some components are missing
    nc <- nc[nc %in% names(cf)]

    cf <- cf[nc]
    vc <- vc[nc, nc]
    
    nc <- sub("model.1_Comp.(.*)_coef.", "C\\1.", nc)
    names(cf) <- colnames(vc) <- rownames(vc) <- nc
    
    if(type=="zero"){
        linfct <- diag(length(cf))
        rownames(linfct) <- nc
    }
    else{
        k <- length(cf)/length(rownames(model@components[[1]][[1]]))
        p <- length(cf) / k
        tmp <- rep(3, k)
        Ktmp <- multcomp::contrMat(tmp, "Tukey")
        linfct <- kronecker(diag(p), Ktmp)
        colnames(linfct) <- names(cf)
        rownames(linfct) <- 1:nrow(linfct)
        for (i in 1:nrow(linfct))
            rownames(linfct)[i] <-
                paste(colnames(linfct)[linfct[i,] == 1],
                      colnames(linfct)[linfct[i,] == -1], sep = "-")
    }
    
    multcomp::glht(multcomp::parm(coef = cf, vcov = vc), linfct = linfct, ...)
})

