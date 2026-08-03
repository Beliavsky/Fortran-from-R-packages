## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")

## ----lp-----------------------------------------------------------------------
library(highs)

L <- c(2, 4, 3)
A <- matrix(c(3, 4, 2,
              2, 1, 2,
              1, 3, 2), nrow = 3, byrow = TRUE)
rhs <- c(60, 40, 80)

sol <- highs_solve(L = L, lower = 0, A = A, rhs = rhs, maximum = TRUE)
sol$objective_value
sol$primal_solution

## ----milp---------------------------------------------------------------------
L <- c(3, 1, 3)
A <- rbind(c(-1,  2,  1),
           c( 0,  4, -3),
           c( 1, -3,  2))
rhs <- c(4, 2, 3)
lower <- c(-Inf, 0, 2)
upper <- c(4, 100, Inf)
types <- c("I", "C", "I")

sol <- highs_solve(L = L, lower = lower, upper = upper,
                   A = A, rhs = rhs, types = types, maximum = TRUE)
sol$objective_value
sol$primal_solution

## ----qp-----------------------------------------------------------------------
Q <- matrix(c(8, 2, 2,
              2, 6, 0,
              2, 0, 4), nrow = 3)
L <- c(-14, -6, -12)

sol <- highs_solve(Q = Q, L = L, lower = 0)
sol$objective_value
sol$primal_solution

## ----control------------------------------------------------------------------
sol <- highs_solve(
  L = c(2, 4, 3), lower = 0,
  A = matrix(c(3, 4, 2, 2, 1, 2, 1, 3, 2), nrow = 3, byrow = TRUE),
  rhs = c(60, 40, 80), maximum = TRUE,
  control = highs_control(
    threads = 1L,
    time_limit = 60,
    log_to_console = FALSE
  )
)
sol$status_message

## ----options, eval = FALSE----------------------------------------------------
# highs_available_solver_options()

## ----sparse-------------------------------------------------------------------
library(Matrix)
A_sparse <- Matrix(c(3, 4, 2, 2, 1, 2, 1, 3, 2),
                   nrow = 3, byrow = TRUE, sparse = TRUE)
sol <- highs_solve(L = c(2, 4, 3), lower = 0,
                   A = A_sparse, rhs = c(60, 40, 80), maximum = TRUE)
sol$objective_value

