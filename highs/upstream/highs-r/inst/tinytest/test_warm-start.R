if (interactive() || any(!c("package:highs", "package:tinytest") %in% search())) {
    library("tinytest")
    library("highs")
}


test_basis_roundtrip <- function() {
    solver <- example_solver("LP")
    hi_solver_run(solver)
    basis <- hi_solver_get_basis(solver)
    expect_true(basis$valid)
    expect_equal(length(basis$col_status), hi_solver_get_num_col(solver))
    expect_equal(length(basis$row_status), hi_solver_get_num_row(solver))

    # Clear solver state and re-set basis
    hi_solver_clear_solver(solver)
    status <- hi_solver_set_basis(solver, basis$col_status, basis$row_status)
    expect_equal(status, 0L)

    # Re-solve should converge in 0 simplex iterations
    hi_solver_run(solver)
    info <- hi_solver_info(solver)
    expect_equal(info$simplex_iteration_count, 0L)
}


test_solution_warm_start <- function() {
    solver <- example_solver("LP")
    hi_solver_run(solver)
    sol <- hi_solver_get_solution(solver)
    obj_val <- hi_solver_info(solver)$objective_function_value

    # Clear solver and set solution
    hi_solver_clear_solver(solver)
    status <- hi_solver_set_solution(solver,
        col_value = sol$col_value,
        row_value = sol$row_value,
        col_dual = sol$col_dual,
        row_dual = sol$row_dual)
    expect_equal(status, 0L)

    # Re-solve and verify same objective
    hi_solver_run(solver)
    info <- hi_solver_info(solver)
    expect_equal(info$objective_function_value, obj_val, tolerance = 1e-8)
}


test_clear_basis <- function() {
    solver <- example_solver("LP")
    hi_solver_run(solver)

    # Clear basis (invalidate it)
    status <- hi_solver_clear_basis(solver)
    expect_equal(status, 0L)

    # Verify basis is now invalid
    basis <- hi_solver_get_basis(solver)
    expect_false(basis$valid)
}


test_basis_status_values <- function() {
    solver <- example_solver("LP")
    hi_solver_run(solver)
    basis <- hi_solver_get_basis(solver)

    # All status values should be in valid range 0-4
    expect_true(all(basis$col_status >= 0L & basis$col_status <= 4L))
    expect_true(all(basis$row_status >= 0L & basis$row_status <= 4L))
}


test_basis_roundtrip()
test_solution_warm_start()
test_clear_basis()
test_basis_status_values()
