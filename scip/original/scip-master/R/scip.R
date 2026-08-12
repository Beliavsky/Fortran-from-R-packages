#' Solve a linear or mixed-integer program using SCIP
#'
#' One-shot interface to the SCIP solver. Formulates and solves:
#' \deqn{\min_{x} \; obj' x}
#' subject to constraint rows defined by \code{A}, \code{b}, \code{sense},
#' with variable types \code{vtype} and bounds \code{lb}, \code{ub}.
#'
#' @param obj Numeric vector of length \code{n}; objective coefficients.
#' @param A Constraint matrix (\code{m x n}). Can be a dense matrix,
#'   \code{dgCMatrix}, or \code{simple_triplet_matrix}.
#' @param b Numeric vector of length \code{m}; constraint right-hand side.
#' @param sense Character vector of length \code{m}; constraint sense.
#'   Each element must be \code{"<="}, \code{">="}, or \code{"=="}.
#' @param vtype Character; variable types. Either a single value applied to all
#'   variables, or a vector of length \code{n}. Values: \code{"C"} (continuous),
#'   \code{"B"} (binary), \code{"I"} (integer). Default \code{"C"}.
#' @param lb Numeric; lower bounds for variables. Single value or vector of
#'   length \code{n}. Default \code{0}.
#' @param ub Numeric; upper bounds for variables. Single value or vector of
#'   length \code{n}. Default \code{Inf}.
#' @param control A list of solver parameters, typically from \code{\link{scip_control}}.
#' @return A named list with components:
#' \describe{
#'   \item{status}{Character; solver status (e.g., "optimal", "infeasible", "unbounded").}
#'   \item{objval}{Numeric; optimal objective value (or \code{NA} if no solution).}
#'   \item{x}{Numeric vector; primal solution (or \code{NULL} if no solution).}
#'   \item{sol_count}{Integer; number of solutions found.}
#'   \item{gap}{Numeric; relative optimality gap.}
#'   \item{info}{List with additional solver information (solve_time, iterations, nodes).}
#' }
#' @export
scip_solve <- function(obj, A, b, sense,
                       vtype = "C", lb = 0, ub = Inf,
                       control = list()) {
    ## Validate dimensions
    n <- length(obj)
    if (is.matrix(A) || inherits(A, "sparseMatrix") || inherits(A, "simple_triplet_matrix")) {
        csc <- make_csc_matrix(A)
    } else {
        stop("A must be a matrix, dgCMatrix, or simple_triplet_matrix")
    }
    m <- csc$nrow
    stopifnot(csc$ncol == n)
    stopifnot(length(b) == m)
    stopifnot(length(sense) == m)
    stopifnot(all(sense %in% c("<=", ">=", "==")))

    ## Expand scalar vtype/lb/ub to vectors
    if (length(vtype) == 1L) vtype <- rep(vtype, n)
    if (length(lb) == 1L) lb <- rep(lb, n)
    if (length(ub) == 1L) ub <- rep(ub, n)
    stopifnot(length(vtype) == n, length(lb) == n, length(ub) == n)
    stopifnot(all(vtype %in% c("C", "B", "I")))

    ## Fix bounds for binary variables
    is_binary <- vtype == "B"
    lb[is_binary] <- pmax(lb[is_binary], 0)
    ub[is_binary] <- pmin(ub[is_binary], 1)

    ## Accept either a scip_control object or a list of arguments
    if (inherits(control, "scip_control")) {
        ctrl <- control
    } else {
        ctrl <- do.call(scip_control, control)
    }

    ## Convert sense to lhs/rhs for SCIP's ranged constraint format
    lhs <- ifelse(sense == "<=", -Inf, b)
    rhs <- ifelse(sense == ">=", Inf, b)

    ## Call C
    result <- .Call(
        R_scip_solve,
        as.double(obj),
        as.integer(csc$i),
        as.integer(csc$p),
        as.double(csc$x),
        as.integer(m),
        as.integer(n),
        as.double(lhs),
        as.double(rhs),
        vtype,
        as.double(lb),
        as.double(ub),
        ctrl,
        PACKAGE = "scip"
    )

    result
}
