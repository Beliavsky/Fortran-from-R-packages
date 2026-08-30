lawSummary <- function(object) {
#
# Function to produce a summary of a Dirichlet (Voronoi)
# tessellation in terms of parameters relevant to Lewis's law
# and Aboav-Weaire's law.  Note that "law" in the function name
# corresponds to "Lewis-Aboav-Weaire.
#
# The parameters of interest are:
# * the areas of each of the interior Dirichlet tiles
# * the number of edges of each of the interior Dirichlet tiles
# * the number of edges of all neighbouring tiles of
#   each of the interior Dirichlet tiles.
#
# This function was created at the request of Kai Xu
# (Fisheries College, Jimei University, Xiamen, Fujian, China 361021).
#

dnbrs <- function(dsgs) {
# Delaunay neighbours.
    iii <- dsgs[,c("ind1","ind2")]
    uind <- with(iii,sort(unique(c(ind1,ind2))))
    rslt <- lapply(uind,function(i,m){sort(c(m[m[,1]==i,2],m[m[,2]==i,1]))},m=iii)
    names(rslt) <- uind
    rslt
}

dsgs <- object$dirsgs
nbrs <- dnbrs(dsgs)

# Layer 1; tiles whose edges have vertices on the boundary.
ex1   <- apply(dsgs[,c("bp1","bp2")],1,any)
dout1 <- dsgs[ex1,]
iout1 <- as.character(unique(c(dout1[,"ind1"],dout1[,"ind2"])))

# Layer 2: tiles having vertices that are Delaunay neighbours
# of vertices of tiles in Layer 1.
nb1 <- as.character(unique(unlist(nbrs[iout1])))

# Layer 3: tiles having vertices that are Delaunay neighbours
# of vertices of tiles in Layer 2.
iout2 <- setdiff(nb1,iout1)
nb2   <- as.character(unique(unlist(nbrs[iout2])))
iout3 <- setdiff(nb2,union(iout1,iout2))

# Keepers.
smry     <- object$summary
nms.all  <- rownames(object$summary)
iout.12  <- union(iout1,iout2)
iout.123 <- union(iout.12,iout3)
nms.12   <- setdiff(nms.all,iout.12)
nms.123  <- setdiff(nms.all,iout.123)
if(!length(nms.123)) return(NULL)

smry.123 <- smry[nms.123,]
smry.12 <- smry[nms.12,]

tile.areas        <- smry.123$dir.area
names(tile.areas) <- nms.123
tile.tags         <- smry.123$z
if(!is.null(tile.tags)) {
    names(tile.tags)  <- nms.123
}
num.edges         <- smry.12$n.tside
names(num.edges)  <- nms.12
i.12              <- as.numeric(nms.12)
i.123             <- as.numeric(nms.123)
nbrs.12           <- lapply(nbrs,function(x,iok){intersect(x,iok)},iok=i.12)
nbrs.123          <- lapply(nbrs.12[i.123],as.character)
num.nbr.edges     <- lapply(nbrs.123,function(k,x){x[k]},x=num.edges)
totnum.nbr.edges  <- sapply(num.nbr.edges,sum)
num.edges         <- num.edges[nms.123]

tl        <- tile.list(object)[i.123]
tv        <- lapply(tl,function(x){data.frame(x=x$x,y=x$y)})
names(tv) <- nms.123
rslt <- list(tile.vertices=tv,tile.areas=tile.areas)
if(!is.null(tile.tags)) {
    rslt <- c(rslt,list(tile.tags=tile.tags))
}
rslt <- c(rslt,list(num.edges=num.edges,num.nbr.edges=num.nbr.edges,
                    totnum.nbr.edges=totnum.nbr.edges))
attr(rslt,"i1") <- as.numeric(iout1)
attr(rslt,"i2") <- as.numeric(iout2)
attr(rslt,"i3") <- as.numeric(iout3)
attr(rslt,"i.kept") <- i.123
rslt
}
