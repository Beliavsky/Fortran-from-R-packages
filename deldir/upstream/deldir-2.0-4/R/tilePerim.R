tilePerim <- function(object,inclbdry=TRUE) {
    if(!inherits(object,"tile.list"))
        stop("Argument \"object\" must be of class \"tile.list\".\n")
    perims   <- lapply(object,tilePerim0,inclbdry=inclbdry)
    perComps <- lapply(perims,function(x){attr(x,"perComps")})
    perims   <- unlist(perims)
    rslt     <- list(perimeters=perims,totalPerim=sum(perims),
                     meanPerim=mean(perims),perComps=perComps)
    rslt
}
