if (interactive() || any(!c("package:highs", "package:tinytest") %in% search())) {
    library("tinytest")
    library("highs")
}


test_unconstrained_qp <- function() {
    L <- c(0, -1, -3)
    Q <- rbind(c(2, 0.0, -1), c(0, 0.2, 0), c(-1, 0.0, 2))
    s <- highs_solve(Q = Q, L = L)
    expect_equal(s[["objective_value"]], -5.5)
}


test_constrained_qp <- function() {
    L <- c(0, -1, -3)
    Q <- rbind(c(2, 0.0, -1), c(0, 0.2, 0), c(-1, 0.0, 2))
    A <- rbind(c(1, 0, 1))
    s <- highs_solve(Q = Q, L = L, lower = 0, A = A, rhs = 2)
    expect_equal(s[["objective_value"]], -5.25)
}


test_unconstrained_qp()
test_constrained_qp()
