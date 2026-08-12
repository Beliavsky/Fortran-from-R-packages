##
## s e t c o v e r . R  Set cover problem
##

setcover <- function(Sets, weights = rep(1, nrow(Sets))) {
    stopifnot(is.numeric(Sets), is.numeric(weights))
    if (!is.matrix(Sets) || !all(Sets %in% c(0, 1))) {
        stop("Argument 'Sets' must be a matrix of zeros and ones.")
    }
    n = nrow(Sets)  # number of sets
    m = ncol(Sets)  # size of universe
    if (any(colSums(Sets) == 0)) {
        stop("Not all elements covered by sets in 'Sets'.")
    }

    sol <- lp(direction = "min",
              objective.in = weights,
              const.mat = Sets,
              const.dir = rep(">=", m),
              const.rhs = rep(1, m),
              transpose.constraints = FALSE,
              all.bin =TRUE)

    sets <- which(sol$solution == 1)
    return(list(sets = sets, objective = sol$objval))
}
