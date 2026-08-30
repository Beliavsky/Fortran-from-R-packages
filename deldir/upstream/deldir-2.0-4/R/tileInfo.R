tileInfo <- function(object,bndry=FALSE,clipp=NULL) {
#
# Function to provide a summary of information about the tiles
# of a Dirichlet (Voronoi) tessellation.
#
# First check that we really are looking at a tessellation/triangulation.
    if (!inherits(object, "deldir")) 
        stop("Argument \"object\" is not of class \"deldir\".\n")
# List the (possibly clipped) tiles in the tessellation.
    tl <- tile.list(object,clipp=clipp)

# Required info:
# * for each tile, the number of edges
# * for each tile, a vector of the lengths of the edges
# * for each tile, the vector of indices of the Delaunay neighbours
#   of the centre of the tile
# * a tabulation of the numbers of edges of tiles
# * a vector of all lengths of edges (with repetitions)
# * a vector of lengths of _unique_ edges
# * the area of each tile

    getEdges <- function(tile){
        ncomp <- attr(tile,"ncomp")
        if(ncomp==1) {
            cmps <- list(tile)
        } else {
            cmps <- tile$tileParts
        }
        edges <- vector("list",ncomp)
        for(i in 1:ncomp) {
            tP <- cmps[[i]]
            x  <- tP$x
            y  <- tP$y
            x1 <- c(x,x[1])
            y1 <- c(y,y[1])
            hedges <- cbind(x1[-1],x,y1[-1],y)
            edges[[i]] <- t(apply(hedges,1,function(x){
                                o <- order(x[1:2],x[3:4])
                                c(x[1:2][o],x[3:4][o])
                                }))
       }
       edges <- do.call(rbind,edges)
       ledge <- apply(edges,1,function(x){sqrt((x[1]-x[2])^2 + (x[3]-x[4])^2)})
       areas <- sapply(cmps,function(x){x$area})
       list(edges=edges,edgeLengths=ledge,area=sum(areas),ptNum=tile$ptNum)
    }

    cnob <- function(tile) {
    # Check not on boundary.
        if(attr(tile,"ncomp")==1) {
            cmps <- list(tile)
        } else {
            cmps <- tile$tileParts
        }
        all(sapply(cmps,function(x){!any(x$bp)}))
    }

    ok <- sapply(tl,if(bndry) function(x){TRUE} else cnob)
    if(sum(ok)==0) {
        whinge <- paste0("All tiles are boundary tiles.  To get a non-vacuous\n",
                         "  result, set bndry=TRUE.\n")
        stop(whinge)
    }
    xxx   <- lapply(tl[ok],getEdges)
    ptNms <- object$summary[["id"]]
    if(is.null(ptNms)) {
        ptNums <- unname(sapply(xxx,function(x){x$ptNum}))
        ptNms  <- paste0("pt.",1:nrow(object$summary))
    }
    nms <- ptNms[ptNms %in% names(ok)][ok]
    names(xxx) <- nms

# Extract and tabulate the edge counts.
    allnedge <- sapply(xxx,function(x){nrow(x$edges)})
    tabnedge <- table(allnedge)

# Extract and combine the edge lengths.  Note that there
# will be duplication since many edges are edges of *two* tiles.
    all.lengths <- unname(unlist(lapply(xxx,function(x){x$edgeLengths})))

# Extract the tile areas into a single vector.
    areas <- sapply(xxx,function(x){x$area})

# The lengths of unique edges.
    all.edges <- do.call(rbind,lapply(xxx,function(x){x$edges}))
    edupe     <- duplicated(all.edges)
    ue        <- all.lengths[!edupe]

# Tile perimeters.
  perims <- tilePerim(tl[ok],inclbdry=bndry)
  names(perims$perimeters) <- nms

# Delaunay neighbours of tile centres.
nbrs   <- getNbrs(object,interior=clipp)
nms.ok <- intersect(nms,names(nbrs))
nbrs   <- nbrs[nms.ok]

# Pack up and go home.
    rslt <- list(indivTiles=xxx,allEdgeCounts=allnedge,tabEdgeCounts=tabnedge,
                 allEdgeLengths=all.lengths,Areas=areas,uniqueEdgeLengths=ue,
                 perimeters=perims,nbrs=nbrs)
    class(rslt) <- "tileInfo"
    rslt
}
