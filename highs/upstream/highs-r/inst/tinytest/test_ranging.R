if (interactive() || any(!c("package:highs", "package:tinytest") %in% search())) {
    library("tinytest")
    library("highs")
}


test_ranging_lp <- function() {
    solver <- example_solver("LP")
    hi_solver_run(solver)
    ranging <- hi_solver_get_ranging(solver)
    expect_true(ranging$valid)
    ncol <- hi_solver_get_num_col(solver)
    nrow <- hi_solver_get_num_row(solver)
    # Cost ranging vectors have length ncol + nrow (includes row duals)
    expect_equal(length(ranging$col_cost_up$value), ncol + nrow)
    expect_equal(length(ranging$col_cost_dn$value), ncol + nrow)
    expect_equal(length(ranging$col_bound_up$value), ncol)
    expect_equal(length(ranging$col_bound_dn$value), ncol)
    expect_equal(length(ranging$row_bound_up$value), nrow)
    expect_equal(length(ranging$row_bound_dn$value), nrow)
    expect_true(all(c("value", "objective", "in_var", "ou_var") %in%
                    names(ranging$col_cost_up)))
}


test_ranging_lp()
