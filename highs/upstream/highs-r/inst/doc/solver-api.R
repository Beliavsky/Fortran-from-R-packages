## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----create-------------------------------------------------------------------
library(highs)

model <- highs_model(
  L = c(2, 4, 3),
  lower = 0,
  A = matrix(c(3, 4, 2, 2, 1, 2, 1, 3, 2), nrow = 3, byrow = TRUE),
  rhs = c(60, 40, 80),
  maximum = TRUE
)
solver <- hi_new_solver(model)

## ----solve--------------------------------------------------------------------
hi_solver_run(solver)
hi_solver_status_message(solver)
sol <- hi_solver_get_solution(solver)
sol$col_value

## ----info---------------------------------------------------------------------
info <- hi_solver_info(solver)
info$objective_function_value
info$simplex_iteration_count
info$primal_solution_status

## ----modify-------------------------------------------------------------------
# Change the objective coefficient for variable 1 (0-based indexing)
hi_solver_set_objective(solver, index = 0L, coeff = 10.0)

# Change variable bounds
hi_solver_set_variable_bounds(solver, index = 0L, lower = 1.0, upper = 20.0)

# Change constraint bounds
hi_solver_set_constraint_bounds(solver, index = 0L, lower = -Inf, upper = 50.0)

# Re-solve
hi_solver_run(solver)
hi_solver_get_solution(solver)$col_value

## ----add----------------------------------------------------------------------
solver <- hi_new_solver(model)

# Add a new variable with bounds [0, 15]
hi_solver_add_vars(solver, lower = 0, upper = 15)
cat("Columns after add:", hi_solver_get_num_col(solver), "\n")

# Add a new constraint: x1 + x2 + x3 + x4 <= 100
hi_solver_add_rows(solver,
  lhs = -Inf, rhs = 100,
  start = 0L,
  index = 0L:3L,
  value = rep(1.0, 4)
)
cat("Rows after add:", hi_solver_get_num_row(solver), "\n")

## ----query--------------------------------------------------------------------
solver <- hi_new_solver(model)
hi_solver_get_lp_costs(solver)
hi_solver_get_num_col(solver)
hi_solver_get_num_row(solver)

## ----wrapper------------------------------------------------------------------
hw <- highs_solver(model)
hw$solve()
hw$status_message()
hw$solution()$col_value
hw$info()$objective_function_value

# Modify and re-solve
hw$L(1, 10.0)     # set cost of variable 1 to 10
hw$vbounds(1, 1, 20)  # set bounds of variable 1 to [1, 20]
hw$solve()
hw$solution()$col_value

## ----io-----------------------------------------------------------------------
solver <- hi_new_solver(model)
temp <- tempfile(fileext = ".mps")
hi_solver_write_model(solver, temp)

# Read into a fresh solver
solver2 <- hi_new_solver(highs_model(L = 0))
hi_solver_read_model(solver2, temp)
hi_solver_run(solver2)
hi_solver_get_solution(solver2)$col_value
unlink(temp)

## ----version------------------------------------------------------------------
solver <- hi_new_solver(model)
hi_solver_version(solver)
hi_solver_run(solver)
hi_solver_get_run_time(solver)

