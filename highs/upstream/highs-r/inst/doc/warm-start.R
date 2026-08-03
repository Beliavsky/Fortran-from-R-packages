## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----basis-roundtrip----------------------------------------------------------
library(highs)

model <- highs_model(
  L = c(2, 4, 3),
  lower = 0,
  A = matrix(c(3, 4, 2, 2, 1, 2, 1, 3, 2), nrow = 3, byrow = TRUE),
  rhs = c(60, 40, 80),
  maximum = TRUE
)
solver <- hi_new_solver(model)

# Solve the original problem
hi_solver_run(solver)
info1 <- hi_solver_info(solver)
cat("First solve:", info1$simplex_iteration_count, "iterations\n")

# Save the basis
basis <- hi_solver_get_basis(solver)
cat("Basis valid:", basis$valid, "\n")
cat("Column statuses:", basis$col_status, "\n")
cat("Row statuses:", basis$row_status, "\n")

## ----basis-restore------------------------------------------------------------
hi_solver_clear_solver(solver)
hi_solver_set_basis(solver, basis$col_status, basis$row_status)
hi_solver_run(solver)
info2 <- hi_solver_info(solver)
cat("Warm-start solve:", info2$simplex_iteration_count, "iterations\n")
cat("Same objective:", info1$objective_function_value == info2$objective_function_value, "\n")

## ----perturbation-------------------------------------------------------------
solver <- hi_new_solver(model)
hi_solver_run(solver)
obj_original <- hi_solver_info(solver)$objective_function_value
cat("Original objective:", obj_original, "\n")

# Save basis before modification
basis <- hi_solver_get_basis(solver)

# Tighten a constraint: rhs from 60 to 50
hi_solver_change_constraint_bounds(solver, idx = 0L, lhs = -Inf, rhs = 50)

# Warm-start from the saved basis
hi_solver_set_basis(solver, basis$col_status, basis$row_status)
hi_solver_run(solver)
info <- hi_solver_info(solver)
cat("After perturbation:", info$objective_function_value,
    "(", info$simplex_iteration_count, "iterations)\n")

## ----solution-warmstart-------------------------------------------------------
solver <- hi_new_solver(model)
hi_solver_run(solver)
sol <- hi_solver_get_solution(solver)

# Clear and warm-start from solution
hi_solver_clear_solver(solver)
hi_solver_set_solution(
  solver,
  col_value = sol$col_value,
  row_value = sol$row_value,
  col_dual  = sol$col_dual,
  row_dual  = sol$row_dual
)
hi_solver_run(solver)
hi_solver_info(solver)$objective_function_value

## ----sparse-solution----------------------------------------------------------
solver <- hi_new_solver(model)

# Set only the non-zero entries (0-based column indices)
hi_solver_set_sparse_solution(solver, index = c(0L, 1L), value = c(5.0, 10.0))
hi_solver_run(solver)
hi_solver_get_solution(solver)$col_value

## ----clear-basis--------------------------------------------------------------
solver <- hi_new_solver(model)
hi_solver_run(solver)
cat("Basis valid after solve:", hi_solver_get_basis(solver)$valid, "\n")

hi_solver_clear_basis(solver)
cat("Basis valid after clear:", hi_solver_get_basis(solver)$valid, "\n")

## ----closure------------------------------------------------------------------
hw <- highs_solver(model)
hw$solve()
basis <- hw$get_basis()
cat("Basis valid:", basis$valid, "\n")

# Perturb and warm-start
hw$cbounds(1, -Inf, 50)
hw$set_basis(basis$col_status, basis$row_status)
hw$solve()
hw$info()$simplex_iteration_count

## ----ranging------------------------------------------------------------------
solver <- hi_new_solver(model)
hi_solver_run(solver)
ranging <- hi_solver_get_ranging(solver)
cat("Ranging valid:", ranging$valid, "\n")

## ----cost-ranging-------------------------------------------------------------
cost_up <- ranging$col_cost_up
data.frame(
  variable = seq_along(cost_up$value),
  max_increase = cost_up$value,
  new_objective = cost_up$objective
)

## ----bound-ranging------------------------------------------------------------
bound_up <- ranging$row_bound_up
data.frame(
  constraint = seq_along(bound_up$value),
  max_increase = bound_up$value,
  new_objective = bound_up$objective
)

## ----presolve-----------------------------------------------------------------
solver <- hi_new_solver(model)
hi_solver_presolve(solver)

