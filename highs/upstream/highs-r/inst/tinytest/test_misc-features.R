if (interactive() || any(!c("package:highs", "package:tinytest") %in% search())) {
    library("tinytest")
    library("highs")
}


test_version <- function() {
    solver <- example_solver("LP")
    v <- hi_solver_version(solver)
    expect_true(is.character(v))
    expect_true(nchar(v) > 0)
    expect_true(grepl("^[0-9]+\\.[0-9]+\\.[0-9]+", v))
}


test_run_time <- function() {
    solver <- example_solver("LP")
    hi_solver_run(solver)
    rt <- hi_solver_get_run_time(solver)
    expect_true(is.numeric(rt))
    expect_true(rt >= 0)
}


test_presolve <- function() {
    solver <- example_solver("LP")
    status <- hi_solver_presolve(solver)
    expect_equal(status, 0L)
}


test_read_model <- function() {
    # Write a model, then read it back into a fresh solver
    solver <- example_solver("LP")
    temp_file <- tempfile(fileext = ".mps")
    hi_solver_write_model(solver, temp_file)

    solver2 <- hi_new_solver(highs_model(L = 0))
    status <- hi_solver_read_model(solver2, temp_file)
    expect_equal(status, 0L)
    expect_equal(hi_solver_get_num_col(solver2), hi_solver_get_num_col(solver))
    expect_equal(hi_solver_get_num_row(solver2), hi_solver_get_num_row(solver))
    unlink(temp_file)
}


test_version()
test_run_time()
test_presolve()
test_read_model()
