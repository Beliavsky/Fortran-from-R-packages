#' Create a SCIP model
#'
#' Creates a new SCIP optimization model for incremental problem construction.
#'
#' @param name Character; problem name. Default \code{"scip_model"}.
#' @return An external pointer representing the SCIP model.
#' @export
scip_model <- function(name = "scip_model") {
    .Call(R_scip_model_create, as.character(name), PACKAGE = "scip")
}

#' Add a variable to a SCIP model
#'
#' @param model A SCIP model (external pointer from \code{\link{scip_model}}).
#' @param obj Numeric; objective coefficient.
#' @param lb Numeric; lower bound. Default \code{0}.
#' @param ub Numeric; upper bound. Default \code{Inf}.
#' @param vtype Character; variable type: \code{"C"} (continuous),
#'   \code{"B"} (binary), or \code{"I"} (integer). Default \code{"C"}.
#' @param name Character; variable name. Default auto-generated.
#' @return Integer; 1-based variable index.
#' @export
scip_add_var <- function(model, obj, lb = 0, ub = Inf,
                         vtype = "C", name = NULL) {
    if (vtype == "B") { lb <- max(lb, 0); ub <- min(ub, 1) }
    if (is.null(name)) {
        nv <- .Call(R_scip_model_get_nvars, model, PACKAGE = "scip")
        name <- paste0("x", nv + 1L)
    }
    .Call(R_scip_model_add_var, model,
          as.double(obj), as.double(lb), as.double(ub),
          as.character(vtype), as.character(name),
          PACKAGE = "scip")
}

#' Add multiple variables to a SCIP model
#'
#' @param model A SCIP model.
#' @param obj Numeric vector; objective coefficients.
#' @param lb Numeric; lower bounds (scalar or vector). Default \code{0}.
#' @param ub Numeric; upper bounds (scalar or vector). Default \code{Inf}.
#' @param vtype Character; variable types (scalar or vector). Default \code{"C"}.
#' @param names Character vector; variable names. Default auto-generated.
#' @return Integer; 1-based index of first variable added.
#' @export
scip_add_vars <- function(model, obj, lb = 0, ub = Inf,
                          vtype = "C", names = NULL) {
    n <- length(obj)
    if (length(lb) == 1L) lb <- rep(lb, n)
    if (length(ub) == 1L) ub <- rep(ub, n)
    if (length(vtype) == 1L) vtype <- rep(vtype, n)
    ## Fix bounds for binary variables
    is_binary <- vtype == "B"
    lb[is_binary] <- pmax(lb[is_binary], 0)
    ub[is_binary] <- pmin(ub[is_binary], 1)
    if (is.null(names)) {
        nv <- .Call(R_scip_model_get_nvars, model, PACKAGE = "scip")
        names <- paste0("x", seq(nv + 1L, nv + n))
    }
    stopifnot(length(lb) == n, length(ub) == n,
              length(vtype) == n, length(names) == n)
    .Call(R_scip_model_add_vars, model,
          as.double(obj), as.double(lb), as.double(ub),
          as.character(vtype), as.character(names),
          PACKAGE = "scip")
}

#' Add a linear constraint to a SCIP model
#'
#' Adds \code{lhs <= sum(coefs * x[vars]) <= rhs}.
#'
#' @param model A SCIP model.
#' @param vars Integer vector; 1-based variable indices.
#' @param coefs Numeric vector; coefficients (same length as \code{vars}).
#' @param lhs Numeric; left-hand side. Default \code{-Inf}.
#' @param rhs Numeric; right-hand side. Default \code{Inf}.
#' @param name Character; constraint name. Default auto-generated.
#' @return Integer; 1-based constraint index.
#' @export
scip_add_linear_cons <- function(model, vars, coefs, lhs = -Inf, rhs = Inf,
                                 name = NULL) {
    if (is.null(name)) {
        nc <- .Call(R_scip_model_get_nconss, model, PACKAGE = "scip")
        name <- paste0("c", nc + 1L)
    }
    .Call(R_scip_model_add_linear_cons, model,
          as.integer(vars), as.double(coefs),
          as.double(lhs), as.double(rhs), as.character(name),
          PACKAGE = "scip")
}

#' Add a quadratic constraint to a SCIP model
#'
#' Adds \code{lhs <= linexpr + quadexpr <= rhs} where
#' \code{quadexpr = sum(quadcoefs[k] * x[quadvars1[k]] * x[quadvars2[k]])}.
#'
#' @param model A SCIP model.
#' @param linvars Integer vector; 1-based variable indices for linear part.
#' @param lincoefs Numeric vector; linear coefficients.
#' @param quadvars1,quadvars2 Integer vectors; 1-based variable indices for
#'   quadratic terms.
#' @param quadcoefs Numeric vector; quadratic coefficients.
#' @param lhs,rhs Numeric; constraint bounds.
#' @param name Character; constraint name.
#' @return Integer; 1-based constraint index.
#' @export
scip_add_quadratic_cons <- function(model,
                                    linvars = integer(0),
                                    lincoefs = double(0),
                                    quadvars1 = integer(0),
                                    quadvars2 = integer(0),
                                    quadcoefs = double(0),
                                    lhs = -Inf, rhs = Inf,
                                    name = NULL) {
    if (is.null(name)) {
        nc <- .Call(R_scip_model_get_nconss, model, PACKAGE = "scip")
        name <- paste0("qc", nc + 1L)
    }
    .Call(R_scip_model_add_quadratic_cons, model,
          as.integer(linvars), as.double(lincoefs),
          as.integer(quadvars1), as.integer(quadvars2), as.double(quadcoefs),
          as.double(lhs), as.double(rhs), as.character(name),
          PACKAGE = "scip")
}

#' Add a SOS1 constraint to a SCIP model
#'
#' At most one variable in the set can be nonzero.
#'
#' @param model A SCIP model.
#' @param vars Integer vector; 1-based variable indices.
#' @param weights Numeric vector; weights determining branching order.
#' @param name Character; constraint name.
#' @return Integer; 1-based constraint index.
#' @export
scip_add_sos1_cons <- function(model, vars, weights = NULL, name = NULL) {
    if (is.null(weights)) weights <- seq_along(vars)
    if (is.null(name)) {
        nc <- .Call(R_scip_model_get_nconss, model, PACKAGE = "scip")
        name <- paste0("sos1_", nc + 1L)
    }
    .Call(R_scip_model_add_sos1_cons, model,
          as.integer(vars), as.double(weights), as.character(name),
          PACKAGE = "scip")
}

#' Add a SOS2 constraint to a SCIP model
#'
#' At most two adjacent variables in the set can be nonzero.
#'
#' @param model A SCIP model.
#' @param vars Integer vector; 1-based variable indices.
#' @param weights Numeric vector; weights determining adjacency order.
#' @param name Character; constraint name.
#' @return Integer; 1-based constraint index.
#' @export
scip_add_sos2_cons <- function(model, vars, weights = NULL, name = NULL) {
    if (is.null(weights)) weights <- seq_along(vars)
    if (is.null(name)) {
        nc <- .Call(R_scip_model_get_nconss, model, PACKAGE = "scip")
        name <- paste0("sos2_", nc + 1L)
    }
    .Call(R_scip_model_add_sos2_cons, model,
          as.integer(vars), as.double(weights), as.character(name),
          PACKAGE = "scip")
}

#' Add an indicator constraint to a SCIP model
#'
#' If \code{binvar = 1} then \code{sum(coefs * x[vars]) <= rhs}.
#'
#' @param model A SCIP model.
#' @param binvar Integer; 1-based index of the binary indicator variable.
#' @param vars Integer vector; 1-based variable indices.
#' @param coefs Numeric vector; coefficients.
#' @param rhs Numeric; right-hand side.
#' @param name Character; constraint name.
#' @return Integer; 1-based constraint index.
#' @export
scip_add_indicator_cons <- function(model, binvar, vars, coefs, rhs,
                                    name = NULL) {
    if (is.null(name)) {
        nc <- .Call(R_scip_model_get_nconss, model, PACKAGE = "scip")
        name <- paste0("ind_", nc + 1L)
    }
    .Call(R_scip_model_add_indicator_cons, model,
          as.integer(binvar), as.integer(vars), as.double(coefs),
          as.double(rhs), as.character(name),
          PACKAGE = "scip")
}

#' Set a SCIP parameter
#'
#' @param model A SCIP model.
#' @param name Character; SCIP parameter name (e.g., \code{"limits/time"}).
#' @param value The parameter value (type is auto-detected by SCIP).
#' @return Invisible \code{NULL}.
#' @export
scip_set_param <- function(model, name, value) {
    .Call(R_scip_model_set_param, model,
          as.character(name), value,
          PACKAGE = "scip")
    invisible(NULL)
}

#' Set objective sense
#'
#' @param model A SCIP model.
#' @param sense Character; \code{"minimize"} (default) or \code{"maximize"}.
#' @return Invisible \code{NULL}.
#' @export
scip_set_objective_sense <- function(model, sense = "minimize") {
    .Call(R_scip_model_set_objective_sense, model,
          as.character(sense),
          PACKAGE = "scip")
    invisible(NULL)
}

#' Solve a SCIP model
#'
#' @param model A SCIP model.
#' @return Invisible \code{NULL}. Use \code{\link{scip_get_status}} and
#'   \code{\link{scip_get_solution}} to retrieve results.
#' @export
scip_optimize <- function(model) {
    .Call(R_scip_model_optimize, model, PACKAGE = "scip")
    invisible(NULL)
}

#' Get solver status
#'
#' @param model A SCIP model (after \code{\link{scip_optimize}}).
#' @return Character; status string (e.g., \code{"optimal"}, \code{"infeasible"}).
#' @export
scip_get_status <- function(model) {
    .Call(R_scip_model_get_status, model, PACKAGE = "scip")
}

#' Get the best solution
#'
#' @param model A SCIP model (after \code{\link{scip_optimize}}).
#' @return A list with \code{objval} and \code{x}.
#' @export
scip_get_solution <- function(model) {
    .Call(R_scip_model_get_solution, model, PACKAGE = "scip")
}

#' Get objective value of best solution
#'
#' @param model A SCIP model (after \code{\link{scip_optimize}}).
#' @return Numeric; objective value, or \code{NA} if no solution.
#' @export
scip_get_objval <- function(model) {
    sol <- .Call(R_scip_model_get_solution, model, PACKAGE = "scip")
    sol$objval
}

#' Get number of solutions found
#'
#' @param model A SCIP model (after \code{\link{scip_optimize}}).
#' @return Integer.
#' @export
scip_get_nsols <- function(model) {
    .Call(R_scip_model_get_nsols, model, PACKAGE = "scip")
}

#' Get the k-th solution from the solution pool
#'
#' @param model A SCIP model.
#' @param k Integer; 1-based solution index (1 = best).
#' @return A list with \code{objval} and \code{x}.
#' @export
scip_get_sol <- function(model, k) {
    .Call(R_scip_model_get_sol, model, as.integer(k), PACKAGE = "scip")
}

#' Get solver information
#'
#' @param model A SCIP model (after \code{\link{scip_optimize}}).
#' @return A list with \code{solve_time}, \code{nodes}, \code{iterations},
#'   \code{gap}, \code{sol_count}.
#' @export
scip_get_info <- function(model) {
    .Call(R_scip_model_get_info, model, PACKAGE = "scip")
}

#' Free a SCIP model
#'
#' Explicitly frees the SCIP model and all associated memory.
#' The model is also freed automatically when garbage collected.
#'
#' @param model A SCIP model.
#' @return Invisible \code{NULL}.
#' @export
scip_model_free <- function(model) {
    .Call(R_scip_model_free, model, PACKAGE = "scip")
    invisible(NULL)
}
