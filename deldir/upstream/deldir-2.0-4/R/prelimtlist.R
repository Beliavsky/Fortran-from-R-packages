prelimtlist <- function(object) {
#
# prelimtlist <--> "preliminary triangle list"
# The ("preliminary") matrix produced by this function may contain
# rows which are indices of the vertices of triangles which are
# formed by the union of three contiguous Delaunay triangles and
# which are not themselves Delaunay triangles.
#
stopifnot(inherits(object, "deldir"))
    a <- object$delsgs[, 5]
    b <- object$delsgs[, 6]
    prelist <- matrix(integer(0), 0, 3)
    for (i in seq(length.out=nrow(object$summary))) {
        jj <- c(b[a == i], a[b == i])
        jj <- sort(unique(jj))
        jj <- jj[jj > i]
        if (length(jj) > 0) {
            for (j in jj) {
                kk <- c(b[a == j], a[b == j])
                kk <- kk[(kk %in% jj) & (kk > j)]
                if (length(kk) > 0) {
                    for (k in kk) {
                        prelist <- rbind(prelist, c(i, j, k))
                    }
                }
            }
        }
    }
    prelist
}
