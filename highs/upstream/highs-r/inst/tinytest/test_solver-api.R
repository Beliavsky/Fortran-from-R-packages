if (interactive() || any(!c("package:highs", "package:tinytest") %in% search())) {
    library("tinytest")
    library("highs")
}
library("checkmate")


test_new_solver <- function() {
    solver <- example_solver()

    expect_class(solver, "externalptr")
    expect_true(!is.null(solver))
}


test_solver_sense <- function() {
    solver <- example_solver()

    hi_solver_set_sense(solver, TRUE)
    expect_true(hi_solver_get_sense(solver))

    hi_solver_set_sense(solver, FALSE)
    expect_false(hi_solver_get_sense(solver))
}


test_solver_offset <- function() {
    solver <- example_solver()

    offset_value <- 5.5
    hi_solver_set_offset(solver, offset_value)
    expect_true(TRUE)
}


test_solver_integrality <- function() {
    solver <- example_solver()

    expect_silent(hi_solver_set_integrality(solver, 0, 1))
    expect_silent(hi_solver_set_integrality(solver, 1, 0))
}


test_solver_objective <- function() {
    solver <- example_solver()

    expect_equal(hi_solver_set_objective(solver, 0, 2.5), 0)
    expect_equal(hi_solver_set_objective(solver, 1, -1.0), 0)
}


test_solver_variable_bounds <- function() {
    solver <- example_solver()

    expect_equal(hi_solver_set_variable_bounds(solver, 0, 1.0, 5.0), 0)
    expect_equal(hi_solver_set_variable_bounds(solver, 1, 0.0, 10.0), 0)
    expect_equal(hi_solver_set_variable_bounds(solver, 0, 5.0, 15.0), 0)
}


test_solver_constraint_bounds <- function() {
    solver <- example_solver()

    expect_silent(hi_solver_set_constraint_bounds(solver, 0, -Inf, 10.0))
    expect_silent(hi_solver_set_constraint_bounds(solver, 0, 0.0, Inf))
}


test_solver_coefficients <- function() {
    solver <- example_solver()

    expect_silent(hi_solver_set_coeff(solver, 0, 0, 2.5))
    expect_silent(hi_solver_set_coeff(solver, 0, 1, 1.5))
}


test_solver_write_model <- function() {
    solver <- example_solver()

    temp_file <- tempfile(fileext = ".mps")
    expect_silent(hi_solver_write_model(solver, temp_file))
    expect_true(file.exists(temp_file))

    unlink(temp_file)
}

test_solver_start <- function() {
    solver <- example_solver()

    expect_silent(hi_solver_set_start(solver, c(0, NA, NA)))
}

test_new_solver()
test_solver_sense()
test_solver_offset()
test_solver_integrality()
test_solver_objective()
test_solver_variable_bounds()
test_solver_constraint_bounds()
test_solver_coefficients()
test_solver_write_model()
test_solver_start()
