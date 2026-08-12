library(scip)

## -----------------------------------------------------------------------
## Test 1: SCIP issue/3047.lp — Pure LP from SCIP's test instances
## min 25*X0 s.t.
##   8*X0 + 4*X1 + 7*X2 + 7*X4 == 25200
##   28*X0 + 39*X1 - 63*X2 - 13*X4 <= 3200
##   -8*X0 + 21*X1 - 7*X2 + 43*X4 <= 4800
## Bounds: X0 in [0,1000], X1 in [0,3000], X2 in [0,2000], X4 in [0,2600]
## -----------------------------------------------------------------------
A <- matrix(c(
    8,   4,   7,   7,
   28,  39, -63, -13,
   -8,  21,  -7,  43
), nrow = 3, byrow = TRUE)
b <- c(25200, 3200, 4800)
res <- scip_solve(
    obj   = c(25, 0, 0, 0),
    A     = A,
    b     = b,
    sense = c("==", "<=", "<="),
    lb    = c(0, 0, 0, 0),
    ub    = c(1000, 3000, 2000, 2600),
    control = list(verbose = FALSE)
)
expect_equal(res$status, "optimal")
expect_true(!is.na(res$objval))
## Verify feasibility: equality constraint should hold
expect_equal(sum(c(8, 4, 7, 7) * res$x), 25200, tolerance = 1e-6)

## -----------------------------------------------------------------------
## Test 2: Simple LP (x + y <= 2 from SCIP's cons.c test, but continuous)
## min -x - y s.t. x + y <= 2, 0 <= x,y <= 2
## Optimal: x=y=1 -> objval=-2  (actually x+y=2, many optima)
## -----------------------------------------------------------------------
A2 <- matrix(c(1, 1), nrow = 1)
b2 <- 2
res2 <- scip_solve(
    obj   = c(-1, -1),
    A     = A2,
    b     = b2,
    sense = "<=",
    lb    = c(0, 0),
    ub    = c(2, 2),
    control = list(verbose = FALSE)
)
expect_equal(res2$status, "optimal")
expect_equal(res2$objval, -2, tolerance = 1e-8)
## x + y should equal 2 at optimum
expect_equal(sum(res2$x), 2, tolerance = 1e-8)

## -----------------------------------------------------------------------
## Test 3: LP with equality and >= constraints
## min x1 + x2 s.t. x1 + x2 == 5, x1 >= 2
## Optimal: x1=2, x2=3, objval=5
## -----------------------------------------------------------------------
A3 <- matrix(c(1, 1, 1, 0), nrow = 2, byrow = TRUE)
b3 <- c(5, 2)
res3 <- scip_solve(
    obj   = c(1, 1),
    A     = A3,
    b     = b3,
    sense = c("==", ">="),
    control = list(verbose = FALSE)
)
expect_equal(res3$status, "optimal")
expect_equal(res3$objval, 5, tolerance = 1e-8)

## -----------------------------------------------------------------------
## Test 4: LP with free variables (lb = -Inf)
## Inspired by SCIP's parsing test: free constraint with unbounded vars
## min x s.t. x + y >= 4, y <= 1, x free, y free
## Optimal: y=1, x=3, objval=3
## -----------------------------------------------------------------------
A4 <- matrix(c(1, 1, 0, 1), nrow = 2, byrow = TRUE)
b4 <- c(4, 1)
res4 <- scip_solve(
    obj   = c(1, 0),
    A     = A4,
    b     = b4,
    sense = c(">=", "<="),
    lb    = c(-Inf, -Inf),
    ub    = c(Inf, Inf),
    control = list(verbose = FALSE)
)
expect_equal(res4$status, "optimal")
expect_equal(res4$objval, 3, tolerance = 1e-8)
