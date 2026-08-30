print.tileInfo <- function(x,digits=4,npl=6,...) {
#
# Matrix of edge lengths:
lel <- lapply(x$indivTiles,function(u){u$edgeLengths})
nel <- sapply(lel,length)
pNs  <- names(x$indivTiles)
m   <- length(nel)
n   <- max(nel)
ld  <- ceiling(log10(max(unlist(lel)))) + digits + 1
fmt <- paste0("%",ld,".",digits,"f")
nms <- paste0("point ",sub("pt.","",pNs),": ")
npad <- max(nchar(nms))
pad  <- paste(rep(" ",npad),collapse="")
cat("\nEdge lengths:\n")
cat("=============\n")
for(i in 1:m) {
    xxx <- sprintf(fmt,lel[[i]])
    irpt <- 0
    repeat{
       buff <- rep(" ",nchar(pad) - nchar(nms[i]))
       buff <- paste(buff,collapse="")
       lhe  <- paste0(nms[i],buff)
       if(irpt==0) cat(lhe ) else cat(pad)
       ibit <- min(npl,length(xxx))
       yyy <- xxx[1:ibit]
       xxx <- xxx[-(1:ibit)]
       cat(yyy,"\n")
       irpt <- 1
       if(length(xxx)==0) break
    }
}

# Table of edgecounts:
tec <- x$tabEdgeCounts
names(attr(tec,"dimnames")) <- ""
mode(tec) <- "character"
tec <- c("  ",tec)
cat("\nTable of edge counts:\n")
cat("=====================\n\n")
print(tec,quote=FALSE)

# Areas:
cat("\nTile areas:\n")
cat("===========\n\n")
print(round(x$Areas,digits=digits))
cat("\n")

# Perimeters:
cat("\nTile perimeters:\n")
cat("================\n\n")
print(round(x$perimeters$perimeters,digits=digits))
cat("\n")

# Neighbours:
cat("\nDelaunay neighbours of tile centres:\n")
cat("=====================================\n\n")
nbrs <- x$nbrs
nms  <- names(nbrs)
nms  <- sub("pt.","",nms)
for(i in 1:length(nbrs)) {
    cat("point ",nms[i]," has neighbours: ",paste(nbrs[[i]],collapse=" "),"\n",sep="")
}
cat("\n")

invisible()
}
