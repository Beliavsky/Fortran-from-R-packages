## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, echo = FALSE------------------------------------------------------
library(piqp)

## -----------------------------------------------------------------------------
P <- matrix(c(6, 0, 0, 4), nrow = 2)
c <- c(-1, -4)
A <- matrix(c(1, -2), nrow = 1)
b <- 1
G <- matrix(c(1, 2, -1, 0), nrow = 2)
h_u <- c(0.2, -1)
x_l <- c(-1, -Inf)  ## 2 variables
x_u <- c(1, Inf)    ## 2 variables

## -----------------------------------------------------------------------------
sol <- solve_piqp(P, c, A, b, G, h_u = h_u, x_l = x_l, x_u = x_u, backend = "auto")
cat(sprintf("(Solution status, description): = (%d, %s)\n",
            sol$status, sol$info$status_desc))
cat(sprintf("Objective: %f, solution: (x1, x2) = (%f, %f)\n", sol$info$primal_obj, sol$x[1], sol$x[2]))

## -----------------------------------------------------------------------------
status_description(sol$status)

## -----------------------------------------------------------------------------
model <- piqp(P, c, A, b, G, h_u = h_u, x_l = x_l, x_u = x_u)
sol2 <- solve(model)
identical(sol, sol2)

## -----------------------------------------------------------------------------
get_dims(model)

## -----------------------------------------------------------------------------
update(model, c = c(-2, -3), h_u = c(0.1, -1.5))
sol3 <- solve(model)
cat(sprintf("Status: %s\n", sol3$info$status_desc))
cat(sprintf("Objective: %f, solution: (x1, x2) = (%f, %f)\n",
            sol3$info$primal_obj, sol3$x[1], sol3$x[2]))

## -----------------------------------------------------------------------------
update(model, A = matrix(c(1, -1), nrow = 1), b = 0,
       x_l = c(-2, -2), x_u = c(2, 2))
sol4 <- solve(model)
cat(sprintf("Status: %s\n", sol4$info$status_desc))
cat(sprintf("Objective: %f, solution: (x1, x2) = (%f, %f)\n",
            sol4$info$primal_obj, sol4$x[1], sol4$x[2]))

## ----error = TRUE-------------------------------------------------------------
try({
update(model, b = c(5, 2))
})

## -----------------------------------------------------------------------------
update_settings(model, new_settings = list(verbose = FALSE, max_iter = 100L))

## -----------------------------------------------------------------------------
sparse_sol <- solve_piqp(P, c, A, b, G, h_u = h_u, x_l = x_l, x_u = x_u, backend = "sparse")
str(sparse_sol)

## -----------------------------------------------------------------------------
P <- matrix(2 * c(3, 0, 0, 2), nrow = 2, ncol = 2)
c <- c(-1, -4)
A <- matrix(c(1, -2), ncol = 2)
b <- 0
x_l <- rep(-1.0, 2)
x_u <- rep(1.0, 2)
sol <- solve_piqp(P = P, c = c, A = A, b = b, x_l = x_l, x_u = x_u)
cat(sprintf("(Solution status, description): = (%d, %s)\n",
            sol$status, sol$info$status_desc))
cat(sprintf("Objective: %f, solution: (x1, x2) = (%f, %f)\n", sol$info$primal_obj, sol$x[1], sol$x[2]))

## -----------------------------------------------------------------------------
G <- diag(2)
h_u <- c(1, 1)
sol <- solve_piqp(P = P, c = c, A = A, b = b, G = G, h_u = h_u,
                  x_l = c(-1, -1), x_u = c(Inf, Inf))
cat(sprintf("(Solution status, description): = (%d, %s)\n",
            sol$status, sol$info$status_desc))
cat(sprintf("Objective: %f, solution: (x1, x2) = (%f, %f)\n", sol$info$primal_obj, sol$x[1], sol$x[2]))

## -----------------------------------------------------------------------------
G <- Matrix::Matrix(c(1, 0, -1, 0, 0, 1, 0, -1), byrow = TRUE,
                    nrow = 4, sparse = TRUE)
h_u <- rep(1, 4)

sol <- solve_piqp(A = A, b = b, c = c, P = P, G = G, h_u = h_u)
cat(sprintf("(Solution status, description): = (%d, %s)\n",
            sol$status, status_description(sol$status)))
cat(sprintf("Objective: %f, solution: (x1, x2) = (%f, %f)\n", sol$info$primal_obj, sol$x[1], sol$x[2]))

## -----------------------------------------------------------------------------
s <- solve_piqp(P = P, c = c, A = A, b = b, G = G, h_u = h_u,
          settings = list(max_iter = 3)) ## Reduced number of iterations
cat(sprintf("(Solution status, description): = (%d, %s)\n",
            s$status, s$info$status_desc))
cat(sprintf("Objective: %f, solution: (x1, x2) = (%f, %f)\n", s$info$primal_obj, s$x[1], s$x[2]))

