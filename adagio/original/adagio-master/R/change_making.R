##
##  c h a n g e _ m a k i n g . R
##


change_making <- function(items, value) {
    stopifnot(is.numeric(items), is.numeric(value))
    if (any(items != floor(items)))
        stop("Argument 'items' must be a vector of integers.")
    if (length(value) != 1 || value != floor(value))
        stop("Argument 'value' must be a single integer.")
    
    n <- length(items)
    obj <- rep(1, n)
    sol <- lp( direction = "min",
               objective.in = obj,
               const.mat = matrix(items, nrow = 1),
               const.dir = "==",
               const.rhs = value,
               all.int = TRUE )

    if (sol$objval == 0) {
        return(list(count = 0, solution = c()))
    } else
        return(list(count = sol$objval, solution = sol$solution))
}
