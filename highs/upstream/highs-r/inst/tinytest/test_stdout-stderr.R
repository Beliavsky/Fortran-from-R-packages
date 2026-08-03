if (interactive() || any(!c("package:highs", "package:tinytest") %in% search())) {
  library("tinytest")
  library("highs")
}

test_hi_solver_run_output_capturable <- function() {
  model <- highs_model(
    L = 0, lower = 0, upper = Inf,
    A = matrix(1, nrow = 1L), lhs = -Inf, rhs = -1
  )
  solver <- hi_new_solver(model)
  hi_solver_set_option(solver, "solver", "simplex")
  hi_solver_set_option(solver, "presolve", FALSE)

  output <- capture.output(hi_solver_run(solver))

  expect_true(is.character(output),
    info = "hi_solver_run output should be capturable as character"
  )
}

test_highs_solve_output_capturable <- function() {
  L <- c(0, 0)
  lower <- c(0, 0)
  upper <- c(10, 10)
  A <- matrix(c(1, 1), nrow = 1)
  lhs <- 1
  rhs <- 15

  output <- capture.output(
    result <- highs_solve(
      L = L, lower = lower, upper = upper,
      A = A, lhs = lhs, rhs = rhs
    )
  )

  expect_true(is.character(output),
    info = "highs_solve output should be capturable as character"
  )

  expect_true(is.list(result),
    info = "highs_solve should return a list"
  )
  expect_true(result$info$valid,
    info = "highs_solve should return valid result information"
  )
}

test_solver_with_log_control <- function() {
  model <- highs_model(
    L = c(1, 1), lower = c(0, 0), upper = c(10, 10),
    A = matrix(c(1, 1), nrow = 1),
    lhs = 1, rhs = 15
  )
  solver <- hi_new_solver(model)

  hi_solver_set_option(solver, "log_to_console", FALSE)

  output <- capture.output(result <- hi_solver_run(solver))

  expect_true(length(output) < 10L,
    info = "With log_to_console=off, output should be minimal"
  )

  expect_equal(result, 0L)
}
