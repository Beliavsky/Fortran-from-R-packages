triMat <- function(object){
    tl <- triang.list(object)
    vl <- lapply(tl,function(x){x[,"ptNum"]})
    vl <- lapply(vl,sort)
    tm <- matrix(unlist(vl),byrow=TRUE,ncol=3)
    return(tm[order(tm[,1]),])
}
