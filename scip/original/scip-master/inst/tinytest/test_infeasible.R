library(scip)

## -----------------------------------------------------------------------
## Test 1: Infeasible — contradictory constraints
## x + y <= 1 AND x + y >= 3
## From SCIP's test philosophy: verify status detection
## -----------------------------------------------------------------------
A <- matrix(c(1, 1, 1, 1), nrow = 2, byrow = TRUE)
b <- c(1, 3)
res <- scip_solve(
    obj   = c(1, 1),
    A     = A,
    b     = b,
    sense = c("<=", ">="),
    control = list(verbose = FALSE)
)
expect_equal(res$status, "infeasible")
expect_true(is.na(res$objval))
expect_null(res$x)
expect_equal(res$sol_count, 0L)

## -----------------------------------------------------------------------
## Test 2: Infeasible — from SCIP's parsing.c infeasible_constraint test
## lhs > rhs: 10 <= 5*x <= 0 (impossible for x >= 0)
## -----------------------------------------------------------------------
A2 <- matrix(c(5), nrow = 1)
b2_upper <- 0
b2_lower <- 10
## lhs=10, rhs=0: this means 10 <= 5x <= 0, infeasible
## We use sense "<=", but need lhs > rhs. Use both <= and >= with contradictory bounds.
## Actually: 5x >= 10 AND 5x <= 0 → x >= 2 AND x <= 0 → infeasible
A2 <- matrix(c(5, 5), nrow = 2, byrow = TRUE)
b2 <- c(0, 10)
res2 <- scip_solve(
    obj   = c(1),
    A     = A2,
    b     = b2,
    sense = c("<=", ">="),
    control = list(verbose = FALSE)
)
expect_equal(res2$status, "infeasible")

## -----------------------------------------------------------------------
## Test 3: Infeasible binary — no feasible binary assignment
## x1 + x2 >= 3, x1, x2 binary (max sum = 2, can't reach 3)
## -----------------------------------------------------------------------
A3 <- matrix(c(1, 1), nrow = 1)
b3 <- 3
res3 <- scip_solve(
    obj   = c(1, 1),
    A     = A3,
    b     = b3,
    sense = ">=",
    vtype = "B",
    control = list(verbose = FALSE)
)
expect_equal(res3$status, "infeasible")
