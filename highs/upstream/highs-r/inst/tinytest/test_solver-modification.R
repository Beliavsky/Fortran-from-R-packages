if (interactive() || any(!c("package:highs", "package:tinytest") %in% search())) {
    library("tinytest")
    library("highs")
}


test_add_vars <- function() {
    solver <- example_solver("LP")
    ncol_before <- hi_solver_get_num_col(solver)
    hi_solver_add_vars(solver, lower = 0, upper = 10)
    expect_equal(hi_solver_get_num_col(solver), ncol_before + 1L)
}


test_add_rows <- function() {
    solver <- example_solver("LP")
    nrow_before <- hi_solver_get_num_row(solver)
    ncol <- hi_solver_get_num_col(solver)
    hi_solver_add_rows(solver,
        lhs = 0, rhs = 100,
        start = 0L, index = 0L:(ncol - 1L),
        value = rep(1.0, ncol))
    expect_equal(hi_solver_get_num_row(solver), nrow_before + 1L)
}


test_change_variable_bounds <- function() {
    solver <- example_solver("LP")
    hi_solver_change_variable_bounds(solver, idx = 0L, lower = 5.0, upper = 50.0)
    hi_solver_run(solver)
    sol <- hi_solver_get_solution(solver)
    expect_true(sol$col_value[1] >= 5.0 || abs(sol$col_value[1] - 5.0) < 1e-8)
}


test_change_constraint_bounds <- function() {
    solver <- example_solver("LP")
    hi_solver_change_constraint_bounds(solver, idx = 0L, lhs = -Inf, rhs = 10.0)
    status <- hi_solver_run(solver)
    expect_equal(status, 0L)
}


test_solver_new_options <- function() {
    ## Options added in HiGHS 1.14.0 -- validates that R/constants.R was
    ## regenerated from v1.14.0's HighsOptions.h.  Restricted to non-PDLP
    ## options because inst/patch_highs.R strips PDLP from the CRAN build.
    opts <- highs_available_solver_options()
    expect_true("iis_time_limit" %in% opts[["option"]])
    expect_true("mip_allow_cut_separation_at_nodes" %in% opts[["option"]])
}


test_options_roundtrip <- function() {
    solver <- example_solver("LP")
    hi_solver_set_option(solver, "time_limit", "99.5")
    expect_equal(hi_solver_get_dbl_option(solver, "time_limit"), 99.5)

    hi_solver_set_option(solver, "output_flag", "FALSE")
    expect_false(hi_solver_get_bool_option(solver, "output_flag"))

    hi_solver_set_option(solver, "threads", "2")
    expect_equal(hi_solver_get_int_option(solver, "threads"), 2L)

    hi_solver_set_option(solver, "solver", "simplex")
    expect_equal(hi_solver_get_str_option(solver, "solver"), "simplex")
}


test_infeasible_model <- function() {
    # Create an infeasible LP: x >= 10 and x <= 5
    L <- 1
    A <- matrix(1, nrow = 1)
    s <- highs_solve(L = L, lower = 10, upper = 5, A = A, lhs = 0, rhs = 100)
    expect_true(s$status != 7L)  # Not optimal
}


test_add_vars()
test_add_rows()
test_change_variable_bounds()
test_change_constraint_bounds()
test_solver_new_options()
test_options_roundtrip()
test_infeasible_model()
