if (interactive() || any(!c("package:highs", "package:tinytest") %in% search())) {
    library("tinytest")
    library("highs")
}


test_lp <- function() {
    L <- c(2, 4, 3)
    A <- matrix(c(3, 4, 2, 2, 1, 2, 1, 3, 2), nrow = 3, byrow = TRUE)
    rhs <- c(60, 40, 80)
    s <- highs_solve(L = L, lower = 0, A = A, rhs = rhs, maximum = TRUE)
    expect_equal(s[["objective_value"]], 230 / 3)

    s <- highs_solve(L = L, lower = 0, maximum = TRUE)
    expect_equal(s[["primal_solution"]], c(0, 0, 0))
}


test_milp <- function() {
    L <- c(3, 1, 3)
    A <- rbind(c(-1,  2,  1), c( 0,  4, -3), c( 1, -3,  2))
    rhs <- c(4, 2, 3)
    lower <- c(-Inf, 0, 2)
    upper <- c(4, 100, Inf)
    types <- c("I", "C", "I")
    s <- highs_solve(L = L, lower = lower, upper = upper, A = A, rhs = rhs,
                     types = types, maximum = TRUE)
    expect_equal(s[["objective_value"]], 23.5)
    isol <- s[["primal_solution"]][types == "I"]
    expect_equal(isol, round(isol))
}


test_milp2 <- function() {
    L <- c(3, 1, 3)
    A <- rbind(c(-1, 2, 1), c(0, 4, -3), c(1, -3, 2))
    rhs <- c(4, 2, 3)
    types <- c("I", "C", "I")
    s <- highs_solve(L = L, lower = 0, A = A, rhs = rhs, types = types,
                     maximum = TRUE, start = c(NA, NA, 3))
    expect_equal(s[["objective_value"]], 26.75)
}


test_lp()
test_milp()
test_milp2()
