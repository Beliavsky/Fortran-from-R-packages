#  Copyright (C) 2022 Martin Maechler
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  https://www.R-project.org/Licenses/


### Compute pam-consistent  medoids  from a data partition,
### such as from a cutree(<hierarchical clustering>)

medoids <- function(x, clustering, diss = inherits(x, "dist"), USE.NAMES=FALSE, ...) {
    ## as for pam(), 'x'  may *either* be a numeric matrix (n x p) *or* a dissimilarity
    diss <- as.logical(diss)
    if(is.na(diss)) stop("'diss' must be one of {TRUE, FALSE}")

    ## split data 'x' according to clustering
    if(diss) {
        ## FIXME!
        stop("'diss = TRUE' is not implemented yet.  Please write to maintainer(\"cluster\")")

        n <- attr(x, "Size")
    } else {
	x <- data.matrix(x)# dropping "automatic rownames" compatibly with daisy()
	if(!is.numeric(x)) stop("'x' is not a numeric dataframe or matrix.")
        n <- NROW(x)
        ##
        xCl <- split.data.frame(x, clustering) # a list of data sets
    }
    ## k <- length(xCl)
    ids <- split(seq_len(n), clustering) # how the 1:n are split in k
    ##  medoids = integer indices of obs., using pam(*, k=1) within each cluster
    idCl <- lapply(xCl, function(d) pam(d, k=1L, ...)$id.med)
    ## select the medoid original 1:n indices:
    vapply(seq_along(ids), function(j) ids[[j]] [idCl[[j]]], 1L)
}

