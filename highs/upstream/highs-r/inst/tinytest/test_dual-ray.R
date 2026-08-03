if (interactive() || any(!c("package:highs", "package:tinytest") %in% search())) {
    library("tinytest")
    library("highs")
}


# A dual (Farkas) ray certifies primal infeasibility and is available only after
# an infeasible LP solved by simplex. Build a trivially infeasible LP: x >= 0 and
# x <= -1.
infeasible_solver <- function() {
    model <- highs_model(L = 0, lower = 0, upper = Inf,
                         A = matrix(1, nrow = 1L), lhs = -Inf, rhs = -1)
    solver <- hi_new_solver(model)
    hi_solver_set_option(solver, "solver", "simplex")
    hi_solver_set_option(solver, "presolve", "off")
    hi_solver_run(solver)
    solver
}


test_dual_ray_on_infeasible_lp <- function() {
    solver <- infeasible_solver()
    ray <- hi_solver_get_dual_ray(solver)
    expect_equal(ray$status, 0L)
    expect_true(ray$has_dual_ray)
    expect_equal(length(ray$dual_ray), hi_solver_get_num_row(solver))
    expect_true(is.numeric(ray$dual_ray))
}


test_no_dual_ray_on_feasible_lp <- function() {
    solver <- example_solver("LP")
    hi_solver_set_option(solver, "solver", "simplex")
    hi_solver_set_option(solver, "presolve", "off")
    hi_solver_run(solver)
    ray <- hi_solver_get_dual_ray(solver)
    expect_equal(ray$status, 0L)
    expect_false(ray$has_dual_ray)
    expect_null(ray$dual_ray)
}


test_primal_ray_structure <- function() {
    # Unbounded LP: minimize -x, x >= 0, single trivially-satisfied row.
    model <- highs_model(L = -1, lower = 0, upper = Inf,
                         A = matrix(0, nrow = 1L), lhs = -Inf, rhs = Inf)
    solver <- hi_new_solver(model)
    hi_solver_set_option(solver, "solver", "simplex")
    hi_solver_set_option(solver, "presolve", "off")
    hi_solver_run(solver)
    ray <- hi_solver_get_primal_ray(solver)
    # Always a well-formed list; the ray itself is present iff HiGHS flagged one.
    expect_equal(ray$status, 0L)
    expect_true(is.logical(ray$has_primal_ray))
    if (isTRUE(ray$has_primal_ray)) {
        expect_equal(length(ray$primal_ray), hi_solver_get_num_col(solver))
    } else {
        expect_null(ray$primal_ray)
    }
}


test_dual_ray_on_infeasible_lp()
test_no_dual_ray_on_feasible_lp()
test_primal_ray_structure()
