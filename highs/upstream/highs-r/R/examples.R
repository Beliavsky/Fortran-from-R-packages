#' Generate Example Optimization Models
#'
#' Creates example optimization models for different problem types:
#' - Linear Programming (LP)
#' - Mixed Integer Linear Programming (MILP)
#' - Quadratic Programming (QP)
#'
#' @param op_type Character string specifying the type of optimization model.
#'        Must be one of "LP", "MILP", or "QP".
#'
#' @return A HiGHS model object configured according to the specified type:
#'   - LP: Maximization problem with 3 variables and 3 constraints
#'   - MILP: Maximization problem with mixed integer and continuous variables
#'   - QP: Problem with quadratic objective function
#'
#' @examples
#' model <- example_model("LP")
#' model <- example_model("MILP")
#' model <- example_model("QP")
#'
#' @export
example_model <- function(op_type = c("LP", "MILP", "QP")) {
  op_type <- match.arg(op_type)
  if (op_type == "LP") {
    L <- c(2, 4, 3)
    A <- matrix(c(3, 4, 2, 2, 1, 2, 1, 3, 2), nrow = 3, byrow = TRUE)
    rhs <- c(60, 40, 80)
    highs_model(L = L, lower = 0, A = A, rhs = rhs, maximum = TRUE)
  } else if (op_type == "MILP") {
    L <- c(3, 1, 3)
    A <- rbind(c(-1, 2, 1), c(0, 4, -3), c(1, -3, 2))
    rhs <- c(4, 2, 3)
    lower <- c(-Inf, 0, 2)
    upper <- c(4, 100, Inf)
    types <- c("I", "C", "I")
    highs_model(
      L = L,
      lower = lower,
      upper = upper,
      A = A,
      rhs = rhs,
      types = types,
      maximum = TRUE
    )
  } else {
    L <- c(0, -1, -3)
    Q <- rbind(c(2, 0.0, -1), c(0, 0.2, 0), c(-1, 0.0, 2))
    A <- cbind(1, 0, 1)
    highs_model(Q = Q, L = L, lower = 0, A = A, rhs = 2)
  }
}


#' Create a HiGHS Solver Object
#'
#' Creates and solves an example optimization model using the HiGHS solver.
#' This is a convenience wrapper that combines model creation and solving
#' in a single function call.
#'
#' @param op_type Character string specifying the type of optimization model.
#'        Must be one of "LP", "MILP", or "QP".
#'
#' @return An object of class \code{"highs_solver"}.
#'
#' @examples
#' solver <- example_solver("LP")
#' solver <- example_solver("MILP")
#' solver <- example_solver("QP")
#'
#' @export
example_solver <- function(op_type = c("LP", "MILP", "QP")) {
  hi_new_solver(example_model(op_type))
}
