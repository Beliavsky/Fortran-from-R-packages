library(scip)

## -----------------------------------------------------------------------
## Test 1: Integer program from SCIP's cons.c unit test
## min -x - y s.t. x + y <= 2, x,y integer in [0,2]
## Optimal: x=2, y=0 or x=0,y=2 or x=1,y=1 — all give objval=-2
## -----------------------------------------------------------------------
A <- matrix(c(1, 1), nrow = 1)
b <- 2
res <- scip_solve(
    obj   = c(-1, -1),
    A     = A,
    b     = b,
    sense = "<=",
    vtype = "I",
    lb    = c(0, 0),
    ub    = c(2, 2),
    control = list(verbose = FALSE)
)
expect_equal(res$status, "optimal")
expect_equal(res$objval, -2)
expect_equal(sum(res$x), 2)
expect_true(all(res$x == floor(res$x)))  # integer solution

## -----------------------------------------------------------------------
## Test 2: Knapsack from SCIP's solveknapsackexactly.c test_general
## 6 items, weights=[7,2,7,5,1,3], capacity=13, all profits=1
## Optimal: take items {1,3,4,5} (0-indexed), value=4
## As MIP: max sum(x_i) s.t. 7x0+2x1+7x2+5x3+x4+3x5 <= 13, binary
## -----------------------------------------------------------------------
A_ks <- matrix(c(7, 2, 7, 5, 1, 3), nrow = 1)
b_ks <- 13
res_ks <- scip_solve(
    obj   = c(-1, -1, -1, -1, -1, -1),
    A     = A_ks,
    b     = b_ks,
    sense = "<=",
    vtype = "B",
    control = list(verbose = FALSE)
)
expect_equal(res_ks$status, "optimal")
expect_equal(res_ks$objval, -4)
## Items 0 and 2 (weights 7 each) should NOT both be taken
expect_true(res_ks$x[1] + res_ks$x[3] <= 1)
## Total weight should be <= 13
expect_true(sum(c(7, 2, 7, 5, 1, 3) * res_ks$x) <= 13)

## -----------------------------------------------------------------------
## Test 3: Knapsack from solveknapsackexactly.c test3 (equal weights)
## 4 items, weights=[2,2,2,2], capacity=4, profits=[1,2,3,4]
## Optimal: take items 2,3 (0-indexed), value=7
## -----------------------------------------------------------------------
A_ks2 <- matrix(c(2, 2, 2, 2), nrow = 1)
b_ks2 <- 4
res_ks2 <- scip_solve(
    obj   = c(-1, -2, -3, -4),
    A     = A_ks2,
    b     = b_ks2,
    sense = "<=",
    vtype = "B",
    control = list(verbose = FALSE)
)
expect_equal(res_ks2$status, "optimal")
expect_equal(res_ks2$objval, -7)
expect_equal(res_ks2$x, c(0, 0, 1, 1))

## -----------------------------------------------------------------------
## Test 4: Knapsack from solveknapsackexactly.c test_greedy1
## 3 items, weights=[1,2,1], capacity=3, profits=[3,2,0.5]
## Optimal: take items 0,1, value=5
## -----------------------------------------------------------------------
A_ks3 <- matrix(c(1, 2, 1), nrow = 1)
b_ks3 <- 3
res_ks3 <- scip_solve(
    obj   = c(-3, -2, -0.5),
    A     = A_ks3,
    b     = b_ks3,
    sense = "<=",
    vtype = "B",
    control = list(verbose = FALSE)
)
expect_equal(res_ks3$status, "optimal")
expect_equal(res_ks3$objval, -5)
expect_equal(res_ks3$x[1], 1)
expect_equal(res_ks3$x[2], 1)

## -----------------------------------------------------------------------
## Test 5: Knapsack from solveknapsackexactly.c test1 (all redundant)
## 2 items, weights=[2,1], capacity=1, profits=[1,-1]
## Optimal: take nothing, value=0
## -----------------------------------------------------------------------
A_ks4 <- matrix(c(2, 1), nrow = 1)
b_ks4 <- 1
res_ks4 <- scip_solve(
    obj   = c(-1, 1),  # profit[0]=1 but too heavy, profit[1]=-1 so don't take
    A     = A_ks4,
    b     = b_ks4,
    sense = "<=",
    vtype = "B",
    control = list(verbose = FALSE)
)
expect_equal(res_ks4$status, "optimal")
expect_equal(res_ks4$objval, 0, tolerance = 1e-8)
expect_equal(res_ks4$x, c(0, 0))

## -----------------------------------------------------------------------
## Test 6: Mixed-integer problem
## min -x1, x1 continuous in [0,2], x2 binary, x1 + x2 <= 1.5
## Optimal: x2=0, x1=1.5, objval=-1.5
## -----------------------------------------------------------------------
A_mi <- matrix(c(1, 1), nrow = 1)
b_mi <- 1.5
res_mi <- scip_solve(
    obj   = c(-1, 0),
    A     = A_mi,
    b     = b_mi,
    sense = "<=",
    vtype = c("C", "B"),
    lb    = c(0, 0),
    ub    = c(2, 1),
    control = list(verbose = FALSE)
)
expect_equal(res_mi$status, "optimal")
expect_equal(res_mi$x[1], 1.5, tolerance = 1e-8)
expect_equal(res_mi$x[2], 0)

## -----------------------------------------------------------------------
## Test 7: flugpl.mps — Classic airline MIP from SCIP test instances
## Known optimal: 1201500
## 18 vars (11 integer), 18 constraints
## We verify just the optimal value (problem is complex to encode by hand,
## so this tests that SCIP handles MIP correctly)
## -----------------------------------------------------------------------
## Encoding flugpl manually is complex; skip for now — covered by
## the simpler MIP tests above which are derived from SCIP's own C tests.
