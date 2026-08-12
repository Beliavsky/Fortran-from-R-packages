library(scip)

## -----------------------------------------------------------------------
## Test 1: Basic LP via model API (from SCIP's cons.c)
## min -x - y s.t. x + y <= 2, x,y in [0,2]
## -----------------------------------------------------------------------
m <- scip_model("test_lp")
v1 <- scip_add_var(m, obj = -1, lb = 0, ub = 2, name = "x")
v2 <- scip_add_var(m, obj = -1, lb = 0, ub = 2, name = "y")
scip_add_linear_cons(m, vars = c(v1, v2), coefs = c(1, 1), rhs = 2)
scip_optimize(m)

expect_equal(scip_get_status(m), "optimal")
sol <- scip_get_solution(m)
expect_equal(sol$objval, -2, tolerance = 1e-8)
expect_equal(sum(sol$x), 2, tolerance = 1e-8)

## -----------------------------------------------------------------------
## Test 2: Integer knapsack from SCIP's solveknapsackexactly.c test_general
## 6 items, weights=[7,2,7,5,1,3], capacity=13, all profits=1
## -----------------------------------------------------------------------
m2 <- scip_model("knapsack")
first <- scip_add_vars(m2,
    obj = rep(-1, 6), lb = 0, ub = 1, vtype = "B",
    names = paste0("item", 1:6))
scip_add_linear_cons(m2, vars = first:(first + 5),
    coefs = c(7, 2, 7, 5, 1, 3), rhs = 13)
scip_optimize(m2)

expect_equal(scip_get_status(m2), "optimal")
expect_equal(scip_get_objval(m2), -4)
sol2 <- scip_get_solution(m2)
expect_true(sum(c(7, 2, 7, 5, 1, 3) * sol2$x) <= 13)
expect_true(scip_get_nsols(m2) >= 1L)

## -----------------------------------------------------------------------
## Test 3: Maximize objective sense
## -----------------------------------------------------------------------
m3 <- scip_model("maximize")
scip_set_objective_sense(m3, "maximize")
scip_add_vars(m3, obj = c(1, 1), lb = 0, ub = 2)
scip_add_linear_cons(m3, vars = 1:2, coefs = c(1, 1), rhs = 3)
scip_optimize(m3)

expect_equal(scip_get_status(m3), "optimal")
expect_equal(scip_get_objval(m3), 3, tolerance = 1e-8)

## -----------------------------------------------------------------------
## Test 4: SOS1 constraint (from SCIP SOS test philosophy)
## At most one of x,y nonzero; max x + y
## -----------------------------------------------------------------------
m4 <- scip_model("sos1")
scip_set_objective_sense(m4, "maximize")
scip_add_vars(m4, obj = c(1, 1), lb = 0, ub = c(3, 5))
scip_add_sos1_cons(m4, vars = 1:2)
scip_optimize(m4)

expect_equal(scip_get_status(m4), "optimal")
sol4 <- scip_get_solution(m4)
expect_equal(sol4$objval, 5, tolerance = 1e-8)
expect_equal(sol4$x[1], 0, tolerance = 1e-8)
expect_equal(sol4$x[2], 5, tolerance = 1e-8)

## -----------------------------------------------------------------------
## Test 5: Quadratic constraint
## min x s.t. x^2 <= 4, x in [-2, 2] => x = -2
## -----------------------------------------------------------------------
m5 <- scip_model("quadratic")
scip_add_var(m5, obj = 1, lb = -2, ub = 2, name = "x")
scip_add_quadratic_cons(m5,
    quadvars1 = 1L, quadvars2 = 1L, quadcoefs = 1.0,
    rhs = 4)
scip_optimize(m5)

expect_equal(scip_get_status(m5), "optimal")
expect_equal(scip_get_objval(m5), -2, tolerance = 1e-8)

## -----------------------------------------------------------------------
## Test 6: Indicator constraint
## z binary, x continuous in [0, 10]
## z = 1 => x <= 1
## min -x - 100*z
## When z=0: x=10, obj=-10; when z=1: x=1, obj=-101 => z=1 optimal
## -----------------------------------------------------------------------
m6 <- scip_model("indicator")
scip_add_var(m6, obj = -100, lb = 0, ub = 1, vtype = "B", name = "z")
scip_add_var(m6, obj = -1, lb = 0, ub = 10, name = "x")
scip_add_indicator_cons(m6, binvar = 1L, vars = 2L, coefs = 1.0, rhs = 1.0)
scip_optimize(m6)

expect_equal(scip_get_status(m6), "optimal")
sol6 <- scip_get_solution(m6)
expect_equal(sol6$x[1], 1)    # z = 1
expect_equal(sol6$x[2], 1, tolerance = 1e-8)  # x = 1
expect_equal(sol6$objval, -101, tolerance = 1e-8)

## -----------------------------------------------------------------------
## Test 7: Set parameters
## -----------------------------------------------------------------------
m7 <- scip_model("params")
scip_add_vars(m7, obj = c(-1, -1), lb = 0, ub = 10, vtype = "I")
scip_add_linear_cons(m7, vars = 1:2, coefs = c(1, 1), rhs = 15)
scip_set_param(m7, "limits/time", 100.0)
scip_set_param(m7, "display/verblevel", 0L)
scip_optimize(m7)

expect_equal(scip_get_status(m7), "optimal")

## -----------------------------------------------------------------------
## Test 8: Solution pool access
## -----------------------------------------------------------------------
m8 <- scip_model("solpool")
scip_add_vars(m8, obj = rep(-1, 4), lb = 0, ub = 1, vtype = "B")
scip_add_linear_cons(m8, vars = 1:4, coefs = rep(1, 4), rhs = 2)
scip_optimize(m8)

nsols <- scip_get_nsols(m8)
expect_true(nsols >= 1L)
sol_best <- scip_get_sol(m8, 1)
expect_equal(sol_best$objval, scip_get_objval(m8))

## -----------------------------------------------------------------------
## Test 9: Get info
## -----------------------------------------------------------------------
info <- scip_get_info(m8)
expect_true(info$solve_time >= 0)
expect_true(info$nodes >= 0)
expect_true(info$sol_count >= 1L)

## -----------------------------------------------------------------------
## Test 10: Explicit free
## -----------------------------------------------------------------------
m9 <- scip_model("to_free")
scip_add_var(m9, obj = 1, lb = 0, ub = 1)
scip_model_free(m9)
expect_error(scip_get_status(m9))  # should error on freed model
