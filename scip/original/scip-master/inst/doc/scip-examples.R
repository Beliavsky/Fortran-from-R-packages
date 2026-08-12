## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
library(scip)

## ----production-lp------------------------------------------------------------
A <- matrix(c(1, 2,
              2, 1), nrow = 2, byrow = TRUE)
b <- c(6, 8)
## scip_solve minimizes, so negate the objective for maximization
res <- scip_solve(obj = c(-5, -4), A = A, b = b,
                  sense = c("<=", "<="),
                  control = list(verbose = FALSE))
res$status
-res$objval  # negate back to get the maximized profit
res$x

## ----knapsack-----------------------------------------------------------------
A <- matrix(c(7, 2, 7, 5, 1, 3), nrow = 1)
res <- scip_solve(obj = c(-1, -1, -1, -1, -1, -1),
                  A = A, b = 13, sense = "<=",
                  vtype = "B",
                  control = list(verbose = FALSE))
res$status
items_packed <- which(res$x == 1)
items_packed
total_weight <- sum(c(7, 2, 7, 5, 1, 3)[items_packed])
total_weight

## ----nqueens-function---------------------------------------------------------
solve_nqueens <- function(n) {
    m <- scip_model(paste0(n, "-queens"))
    scip_set_objective_sense(m, "maximize")

    ## Create n*n binary variables x[i,j] with objective 1
    ## Variable (i,j) is stored at index (i-1)*n + j
    idx <- function(i, j) (i - 1L) * n + j
    scip_add_vars(m, obj = rep(1, n * n), lb = 0, ub = 1, vtype = "B",
                  names = sprintf("x%d_%d", rep(1:n, each = n), rep(1:n, n)))

    ## Row constraints: exactly one queen per row
    for (i in 1:n) {
        scip_add_linear_cons(m, vars = idx(i, 1:n), coefs = rep(1, n),
                             lhs = 1, rhs = 1, name = sprintf("row_%d", i))
    }

    ## Column constraints: exactly one queen per column
    for (j in 1:n) {
        scip_add_linear_cons(m, vars = idx(1:n, j), coefs = rep(1, n),
                             lhs = 1, rhs = 1, name = sprintf("col_%d", j))
    }

    ## Diagonal constraints: at most one queen per diagonal
    ## Down-right diagonals
    for (d in (-(n - 2)):(n - 2)) {
        rows <- max(1, 1 + d):min(n, n + d)
        cols <- rows - d
        if (length(rows) >= 2) {
            scip_add_linear_cons(m, vars = idx(rows, cols),
                                 coefs = rep(1, length(rows)),
                                 rhs = 1, name = sprintf("diag_down_%d", d))
        }
    }
    ## Down-left diagonals (anti-diagonals)
    for (s in 3:(2 * n - 1)) {
        rows <- max(1, s - n):min(n, s - 1)
        cols <- s - rows
        if (length(rows) >= 2) {
            scip_add_linear_cons(m, vars = idx(rows, cols),
                                 coefs = rep(1, length(rows)),
                                 rhs = 1, name = sprintf("diag_up_%d", s))
        }
    }

    scip_optimize(m)
    sol <- scip_get_solution(m)

    ## Extract queen positions
    x_mat <- matrix(round(sol$x), nrow = n, byrow = TRUE)
    positions <- which(x_mat == 1, arr.ind = TRUE)
    colnames(positions) <- c("row", "col")
    list(status = scip_get_status(m),
         n_queens = as.integer(sol$objval),
         positions = positions[order(positions[, "row"]), ],
         board = x_mat)
}

## ----nqueens-solve------------------------------------------------------------
result <- solve_nqueens(8)
result$status
result$n_queens
result$positions

## ----nqueens-board------------------------------------------------------------
board_str <- apply(result$board, 1, function(row) {
    paste(ifelse(row == 1, "Q", "."), collapse = " ")
})
cat(board_str, sep = "\n")

## ----circlepacking------------------------------------------------------------
radii <- c(1.0, 1.5, 0.8)
n <- length(radii)
H <- 5.0  # fixed height

m <- scip_model("circlepacking")

## Variables: x_i, y_i for each circle, plus W
x_idx <- integer(n)
y_idx <- integer(n)
for (i in seq_len(n)) {
    x_idx[i] <- scip_add_var(m, obj = 0, lb = radii[i], ub = 100,
                              name = sprintf("x%d", i))
    y_idx[i] <- scip_add_var(m, obj = 0, lb = radii[i], ub = H - radii[i],
                              name = sprintf("y%d", i))
}
w_idx <- scip_add_var(m, obj = 1, lb = 0, ub = 100, name = "W")

## x_i <= W - r_i  =>  x_i - W <= -r_i
for (i in seq_len(n)) {
    scip_add_linear_cons(m, vars = c(x_idx[i], w_idx), coefs = c(1, -1),
                         rhs = -radii[i],
                         name = sprintf("fit_x_%d", i))
}

## Non-overlap: (x_i - x_j)^2 + (y_i - y_j)^2 >= (r_i + r_j)^2
## Expand: x_i^2 - 2*x_i*x_j + x_j^2 + y_i^2 - 2*y_i*y_j + y_j^2 >= (r_i+r_j)^2
for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
        min_dist_sq <- (radii[i] + radii[j])^2
        scip_add_quadratic_cons(m,
            quadvars1 = c(x_idx[i], x_idx[i], x_idx[j], y_idx[i], y_idx[i], y_idx[j]),
            quadvars2 = c(x_idx[i], x_idx[j], x_idx[j], y_idx[i], y_idx[j], y_idx[j]),
            quadcoefs = c(1, -2, 1, 1, -2, 1),
            lhs = min_dist_sq,
            name = sprintf("nooverlap_%d_%d", i, j))
    }
}

scip_optimize(m)
cat("Status:", scip_get_status(m), "\n")
sol <- scip_get_solution(m)
W_opt <- sol$x[w_idx]
cat(sprintf("Minimum width: %.3f (height fixed at %.1f)\n", W_opt, H))
for (i in seq_len(n)) {
    cat(sprintf("  Circle %d (r=%.1f): center = (%.3f, %.3f)\n",
                i, radii[i], sol$x[x_idx[i]], sol$x[y_idx[i]]))
}

## ----facility-location--------------------------------------------------------
## Problem data
p <- 3  # warehouses
q <- 4  # customers
fixed_cost <- c(100, 150, 120)
capacity <- c(50, 60, 40)
demand <- c(20, 25, 15, 30)
## Shipping cost matrix (warehouses x customers)
ship_cost <- matrix(c(
    4, 8, 5, 6,
    6, 3, 7, 4,
    5, 5, 4, 8
), nrow = p, byrow = TRUE)

m <- scip_model("facility_location")

## Binary variables y_i: open warehouse i?
y_idx <- integer(p)
for (i in 1:p) {
    y_idx[i] <- scip_add_var(m, obj = fixed_cost[i], vtype = "B",
                              name = sprintf("y%d", i))
}

## Continuous variables x_ij: amount shipped from i to j
x_idx <- matrix(0L, nrow = p, ncol = q)
for (i in 1:p) {
    for (j in 1:q) {
        x_idx[i, j] <- scip_add_var(m, obj = ship_cost[i, j],
                                     lb = 0, ub = max(capacity),
                                     name = sprintf("x%d_%d", i, j))
    }
}

## Demand constraints: sum_i x_ij >= d_j
for (j in 1:q) {
    scip_add_linear_cons(m, vars = x_idx[, j], coefs = rep(1, p),
                         lhs = demand[j],
                         name = sprintf("demand_%d", j))
}

## Capacity constraints: sum_j x_ij - k_i * y_i <= 0
for (i in 1:p) {
    scip_add_linear_cons(m,
                         vars = c(x_idx[i, ], y_idx[i]),
                         coefs = c(rep(1, q), -capacity[i]),
                         rhs = 0,
                         name = sprintf("capacity_%d", i))
}

scip_optimize(m)
sol <- scip_get_solution(m)
cat("Status:", scip_get_status(m), "\n")
cat("Total cost:", sol$objval, "\n")

open_warehouses <- which(round(sol$x[y_idx]) == 1)
cat("Open warehouses:", open_warehouses, "\n")

shipments <- matrix(sol$x[x_idx], nrow = p)
cat("\nShipment plan (rows=warehouses, cols=customers):\n")
rownames(shipments) <- paste0("W", 1:p)
colnames(shipments) <- paste0("C", 1:q)
print(round(shipments, 1))

## ----indicator----------------------------------------------------------------
revenue <- c(12, 8, 15)
max_prod <- c(10, 20, 8)
setup_cost <- 50

m <- scip_model("setup_production")
scip_set_objective_sense(m, "maximize")

z_idx <- integer(3)
x_idx <- integer(3)
for (i in 1:3) {
    z_idx[i] <- scip_add_var(m, obj = -setup_cost, vtype = "B",
                              name = sprintf("z%d", i))
    x_idx[i] <- scip_add_var(m, obj = revenue[i], lb = 0, ub = max_prod[i],
                              name = sprintf("x%d", i))
    ## x_i <= max_prod[i] * z_i  =>  x_i - max_prod[i] * z_i <= 0
    scip_add_linear_cons(m, vars = c(x_idx[i], z_idx[i]),
                         coefs = c(1, -max_prod[i]),
                         rhs = 0,
                         name = sprintf("link_%d", i))
}

scip_optimize(m)
sol <- scip_get_solution(m)
cat("Status:", scip_get_status(m), "\n")
cat("Profit:", sol$objval, "\n")
for (i in 1:3) {
    on <- round(sol$x[z_idx[i]])
    prod <- sol$x[x_idx[i]]
    cat(sprintf("  Machine %d: %s, produce %.0f units\n",
                i, if (on) "ON" else "OFF", prod))
}

## ----controls-----------------------------------------------------------------
## One-shot: time limit and gap tolerance
ctrl <- scip_control(verbose = FALSE, time_limit = 60, gap_limit = 0.01)

## Model-building: set SCIP parameters directly
m <- scip_model("tuning_example")
scip_set_param(m, "display/verblevel", 0L)   # suppress output
scip_set_param(m, "limits/time", 60.0)        # 60-second time limit
scip_set_param(m, "limits/gap", 0.01)         # 1% optimality gap

