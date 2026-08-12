# w <- c( 40,  60, 30, 40, 20, 5)
# p <- c(110, 150, 70, 80, 30, 5)
# cap <- 90  # c(65, 85)

mknapsack <- function(w, p, cap) {
    stopifnot(is.numeric(p), is.numeric(w), is.numeric(cap))
    if (any(w <= 0) || any(p <= 0))
        stop("Profits 'w' and weights 'p' must be vectors of positive numbers.")

    m <- length(cap)
    n <- length(p)  # == length(w)
    if (length(w) != n)
        stop("Profits 'p' and weights 'w' must be vectors of equal length.")

    if (n == 1) {
        inds <- which(cap >= w)
        if (length(inds) == 0) {
            return(list(ksack = 0, val = 0))
        } else {
            ind <- inds[1]
            return(list(ksack = ind, val = p, bs = 0))
        }
    }

    if (m == 1) {
        sol <- lp(direction = "max",
                  objective.in = p,
                  const.mat = matrix(w, nrow = 1),
                  const.dir = "<=",
                  const.rhs = cap,
                  all.bin = TRUE)
        return(list(ksack = sol$solution, val = sol$objval, bs = 0))
    }

    obj <- rep(p, m)
    cm1 <- matrix(0, nrow = m, ncol = m * n)
    for (k in 1:m) {
        cm1[k, ((k-1)*n+1):(k*n)] <- w
    }
    cm2 <- diag(1, n)
    for (k in 2:m) {
        cm2 <- cbind(cm2, diag(1, n))
    }
    cm <- rbind(cm1, cm2)

    sol <- lp(direction = "max",
              objective.in = obj,
              const.mat = cm,
              const.dir = rep("<=", n + m),
              const.rhs = c(cap, rep(1, n)),
              all.bin = TRUE)

    sls = sol$solution
    my = sls[1:n]
    for (k in 2:m) {
        my = my + k * sls[((k-1)*n+1):(k*n)]
    }

    return(list(ksack = my, val = sol$objval, bs = 0))
}

