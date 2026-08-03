#' Create a new solver instance.
#'
#' This function initializes a new Highs solver instance using the provided model pointer.
#'
#' @param model an object of class \code{"highs_model"}
#'
#' @return A new solver instance.
#'
#' @examples
#' model <- example_model()
#' solver <- hi_new_solver(model)
#' @export
hi_new_solver <- function(model) {
  checkmate::assert_class(model, classes = "highs_model")
  msg <- capture.output(result <- new_solver(model))
  if (result$status == -1L) {
    stop(msg)
  } else if (result$status == 1L) {
    warning(msg)
  }
  solver <- result$solver
  if (is.null(solver)) {
    stop(sprintf("%s\nSolver could not be created check your model!", msg))
  }
  class(solver) <- c("highs_solver", class(solver))
  solver
}


#' Get the optimization sense of the solver instance.
#'
#' This function returns the optimization sense (e.g., minimization or maximization) of the provided solver instance.
#'
#' @param solver An object of class "highs_solver" representing the solver instance.
#'
#' @return The optimization sense of the solver instance.
#'
#' @examples
#' solver <- example_solver()
#' hi_solver_get_sense(solver)
#'
#' @export
hi_solver_get_sense <- function(solver) {
  checkmate::assert_class(solver, classes = "highs_solver")
  as.logical(solver_get_sense(solver))
}

#' Set the optimization sense of the solver instance.
#'
#' This function updates the optimization sense for the given solver instance. Use TRUE for maximization and FALSE for minimization.
#'
#' @param solver An object of class "highs_solver".
#' @param maximum A boolean indicating whether to set maximization (TRUE) or minimization (FALSE).
#'
#' @return The updated solver instance with the new optimization sense.
#'
#' @examples
#' solver <- example_solver()
#' hi_solver_set_sense(solver, TRUE)
#'
#' @export
hi_solver_set_sense <- function(solver, maximum) {
  checkmate::assert_class(solver, classes = "highs_solver")
  checkmate::assert_logical(maximum, len = 1, any.missing = FALSE)
  solver_set_sense(solver, maximum)
}


#' Set the objective offset for the solver.
#'
#' This function sets the objective offset in the solver instance.
#'
#' @param solver An object of class "highs_solver".
#' @param ext_offset A numeric value representing the offset.
#'
#' @return The solver instance with the updated offset.
#'
#' @examples
#' solver <- example_solver()
#' hi_solver_set_offset(solver, 5.0)
#'
#' @export
hi_solver_set_offset <- function(solver, ext_offset) {
  checkmate::assert_class(solver, classes = "highs_solver")
  checkmate::assert_numeric(ext_offset, len = 1, any.missing = FALSE)
  solver_set_offset(solver, as.double(ext_offset))
}


#' Set integrality for a set of variables in the solver.
#'
#' This function defines whether a variable is categorized as integral or continuous.
#'
#' @param solver An object of class "highs_solver".
#' @param index An integer vector specifying the variable index.
#' @param type An integer vector representing the integrality type.
#'
#' @return The solver instance with updated integrality settings.
#'
#' @examples
#' solver <- example_solver()
#' hi_solver_set_integrality(solver, 1, 1)
#'
#' @export
hi_solver_set_integrality <- function(solver, index, type) {
  checkmate::assert_class(solver, classes = "highs_solver")
  checkmate::assert_integerish(
    index,
    lower = 0,
    any.missing = FALSE,
    min.len = 1
  )
  checkmate::assert_integerish(
    type,
    lower = 0,
    any.missing = FALSE,
    len = length(index)
  )
  solver_set_integrality(solver, as.integer(index), as.integer(type))
}


#' Set the objective coefficient for a set of variables.
#'
#' This function assigns a coefficient to a variable in the objective function.
#'
#' @param solver An object of class "highs_solver".
#' @param index The variable index.
#' @param coeff A numeric value representing the objective coefficient.
#'
#' @return The solver instance with the updated objective.
#'
#' @examples
#' solver <- example_solver()
#' hi_solver_set_objective(solver, 2, 3.5)
#'
#' @export
hi_solver_set_objective <- function(solver, index, coeff) {
  checkmate::assert_class(solver, classes = "highs_solver")
  checkmate::assert_integerish(
    index,
    lower = 0,
    any.missing = FALSE,
    min.len = 1
  )
  checkmate::assert_numeric(coeff, len = length(index), any.missing = FALSE)
  solver_set_objective(solver, as.integer(index), as.double(coeff))
}


#' Set variable bounds for a set of variables.
#'
#' This function sets the lower and upper bounds for a set of variables.
#'
#' @param solver An object of class "highs_solver".
#' @param index The variable index.
#' @param lower The lower bound.
#' @param upper The upper bound.
#'
#' @return The solver instance with updated variable bounds.
#'
#' @examples
#' solver <- example_solver()
#' hi_solver_set_variable_bounds(solver, 2, 0, 10)
#'
#' @export
hi_solver_set_variable_bounds <- function(solver, index, lower, upper) {
  checkmate::assert_class(solver, classes = "highs_solver")
  checkmate::assert_integerish(
    index,
    lower = 0,
    any.missing = FALSE,
    min.len = 1
  )
  checkmate::assert_numeric(lower, len = length(index), any.missing = FALSE)
  checkmate::assert_numeric(upper, len = length(index), any.missing = FALSE)
  solver_set_variable_bounds(
    solver,
    as.integer(index),
    as.double(lower),
    as.double(upper)
  )
}


#' Set constraint bounds for a given constraint.
#'
#' This function sets the lower and upper bounds for a specific constraint.
#'
#' @param solver An object of class "highs_solver".
#' @param index The constraint index.
#' @param lower The lower bound.
#' @param upper The upper bound.
#'
#' @return The solver instance with updated constraint bounds.
#'
#' @examples
#' solver <- example_solver()
#' hi_solver_set_constraint_bounds(solver, 1, -Inf, 100)
#'
#' @export
hi_solver_set_constraint_bounds <- function(solver, index, lower, upper) {
  checkmate::assert_class(solver, classes = "highs_solver")
  checkmate::assert_integerish(index, lower = 0, any.missing = FALSE, len = 1)
  checkmate::assert_numeric(lower, len = 1, any.missing = FALSE)
  checkmate::assert_numeric(upper, len = 1, any.missing = FALSE)
  solver_set_constraint_bounds(
    solver,
    as.integer(index),
    as.double(lower),
    as.double(upper)
  )
}


#' Set a coefficient in the constraint matrix.
#'
#' This function assigns a coefficient value to a specific entry in the constraint matrix.
#'
#' @param solver An object of class "highs_solver".
#' @param row The row index.
#' @param col The column index.
#' @param val The coefficient value.
#'
#' @return The solver instance with the updated coefficient.
#'
#' @examples
#' solver <- example_solver()
#' hi_solver_set_coeff(solver, 1, 1, 4.2)
#'
#' @export
hi_solver_set_coeff <- function(solver, row, col, val) {
  checkmate::assert_class(solver, classes = "highs_solver")
  checkmate::assert_integerish(row, lower = 0, any.missing = FALSE, len = 1)
  checkmate::assert_integerish(col, lower = 0, any.missing = FALSE, len = 1)
  checkmate::assert_numeric(val, len = 1, any.missing = FALSE)
  solver_set_coeff(solver, as.integer(row), as.integer(col), as.double(val))
}


#' Add Variables to the Solver
#'
#' This function adds new variables to the solver with specified bounds.
#'
#' @param solver An object of class "highs_solver".
#' @param lower A numeric vector of lower bounds for the new variables.
#' @param upper A numeric vector of upper bounds for the new variables.
#'
#' @return The solver instance with the new variables added.
#'
#' @examples
#' solver <- example_solver()
#' hi_solver_add_vars(solver, lower = c(0, 0, 0), upper = c(10, 10, 10))
#'
#' @export
hi_solver_add_vars <- function(solver, lower, upper) {
  checkmate::assert_class(solver, classes = "highs_solver")
  checkmate::assert_numeric(lower, any.missing = FALSE)
  checkmate::assert_numeric(upper, any.missing = FALSE)
  checkmate::assert_true(length(lower) == length(upper))
  solver_add_vars(solver, as.double(lower), as.double(upper))
}


#' Clear All Solver Data
#'
#' This function clears all data from the solver instance, including the model and solution.
#'
#' @param solver An object of class "highs_solver".
#'
#' @return The cleared solver instance.
#'
#' @examples
#' solver <- example_solver()
#' hi_solver_clear(solver)
#'
#' @export
hi_solver_clear <- function(solver) {
  checkmate::assert_class(solver, classes = "highs_solver")
  solver_clear(solver)
}


#' Clear Model Data
#'
#' This function clears only the model data from the solver instance.
#'
#' @param solver An object of class "highs_solver".
#'
#' @return The solver instance with cleared model data.
#'
#' @examples
#' solver <- example_solver()
#' hi_solver_clear_model(solver)
#'
#' @export
hi_solver_clear_model <- function(solver) {
  checkmate::assert_class(solver, classes = "highs_solver")
  solver_clear_model(solver)
}


#' Clear Solver State
#'
#' This function clears the internal solver state while preserving the model.
#'
#' @param solver An object of class "highs_solver".
#'
#' @return The solver instance with cleared solver state.
#'
#' @examples
#' solver <- example_solver()
#' hi_solver_clear_solver(solver)
#'
#' @export
hi_solver_clear_solver <- function(solver) {
  checkmate::assert_class(solver, classes = "highs_solver")
  solver_clear_solver(solver)
}

#' Set Starting Solution for the Solver
#'
#' This function sets the starting solution for the solver.
#'
#' @param solver An object of class "highs_solver".
#' @param val A numeric vector of values for the variables.
#' Missing values can be used to specify that particular
#' variables should have their starting values calculated automatically.
#'
#' @return An integer value indicating the status.
#'
#' @examples
#' solver <- example_solver()
#' hi_solver_set_start(solver, c(0, NA, NA))
#'
#' @export
hi_solver_set_start <- function(solver, val) {
  checkmate::assert_class(solver, classes = "highs_solver")
  checkmate::assert_numeric(val, any.missing = TRUE, all.missing = FALSE)
  idx <- which(!is.na(val))
  solver_set_solution_vec(solver, idx - 1L, val[idx])
}

#' Run the Solver
#'
#' Executes the optimization solver on the current model.
#'
#' Solver output is controlled by the `output_flag` option:
#' use `hi_solver_set_option(solver, "output_flag", "off")` to suppress logging.
#'
#' @param solver An object of class "highs_solver".
#' @param verbose a logical if TRUE, prints solver log messages to the console.
#'
#' @return An integer status code (-1 = error, 0 = success, 1 = warning).
#'
#' @examples
#' solver <- example_solver()
#' hi_solver_run(solver)
#'
#' @export
hi_solver_run <- function(solver, verbose = FALSE) {
  checkmate::assert_class(solver, classes = "highs_solver")
  msgs <- capture.output(status <- solver_run(solver))
  if (isTRUE(verbose)) {
    writeLines(msgs)
  }
  status
}


# Get Model from Solver
#
# This function retrieves the current optimization model from the solver.
#
# @param solver An object of class "highs_solver".
#
# @return The current optimization model.
#
# @examples
# solver <- example_solver()
# model <- hi_solver_get_model(solver)
#
hi_solver_get_model <- function(solver) {
  checkmate::assert_class(solver, classes = "highs_solver")
  solver_get_model(solver)
}


#' Get Number of Variables
#'
#' This function returns the number of variables (columns) in the optimization model.
#'
#' @param solver An object of class "highs_solver".
#'
#' @return An integer representing the number of variables.
#'
#' @examples
#' solver <- example_solver()
#' n_vars <- hi_solver_get_num_col(solver)
#'
#' @export
hi_solver_get_num_col <- function(solver) {
  checkmate::assert_class(solver, classes = "highs_solver")
  solver_get_num_col(solver)
}


#' Get Number of Constraints
#'
#' This function returns the number of constraints (rows) in the optimization model.
#'
#' @param solver An object of class "highs_solver".
#'
#' @return An integer representing the number of constraints.
#'
#' @examples
#' solver <- example_solver()
#' n_constraints <- hi_solver_get_num_row(solver)
#'
#' @export
hi_solver_get_num_row <- function(solver) {
  checkmate::assert_class(solver, classes = "highs_solver")
  solver_get_num_row(solver)
}


#' Write Model to File
#'
#' This function writes the current optimization model to a file.
#'
#' @param solver An object of class "highs_solver".
#' @param filename A character string specifying the output file path.
#'
#' @return Invisible NULL.
#'
#' @examples
#' solver <- example_solver()
#' model_file <- tempfile(fileext = ".mps")
#' hi_solver_write_model(solver, model_file)
#'
#' @export
hi_solver_write_model <- function(solver, filename) {
  checkmate::assert_class(solver, classes = "highs_solver")
  solver_write_model(solver, filename)
}


#' Write Basis to File
#'
#' This function writes the current basis information to a file.
#'
#' @param solver An object of class "highs_solver".
#' @param filename A character string specifying the output file path.
#'
#' @return Invisible NULL.
#'
#' @examples
#' solver <- example_solver()
#' basis_file <- tempfile(fileext = ".txt")
#' hi_solver_write_basis(solver, basis_file)
#'
#' @export
hi_solver_write_basis <- function(solver, filename) {
  checkmate::assert_class(solver, classes = "highs_solver")
  solver_write_basis(solver, filename)
}


#' Get Solver Status Message
#'
#' This function returns a human-readable message describing the current solver status.
#'
#' @param solver An object of class "highs_solver".
#'
#' @return A character string containing the status message.
#'
#' @examples
#' solver <- example_solver()
#' hi_solver_run(solver)
#' message <- hi_solver_status_message(solver)
#'
#' @export
hi_solver_status_message <- function(solver) {
  checkmate::assert_class(solver, classes = "highs_solver")
  solver_status_message(solver)
}


#' Get Solver Status
#'
#' This function returns the current status of the solver.
#'
#' @param solver An object of class "highs_solver".
#'
#' @return A status code indicating the solver state.
#'
#' @examples
#' solver <- example_solver()
#' hi_solver_run(solver)
#' status <- hi_solver_status(solver)
#'
#' @export
hi_solver_status <- function(solver) {
  checkmate::assert_class(solver, classes = "highs_solver")
  solver_status(solver)
}


#' Get Solver Infinity Value
#'
#' This function returns the value that the solver uses to represent infinity.
#'
#' @return A numeric value representing infinity in the solver.
#'
#' @examples
#' inf <- hi_solver_infinity()
#'
#' @export
hi_solver_infinity <- function() {
  solver_infinity()
}


#' Reset Global Scheduler
#'
#' This function resets the global scheduler used by the solver.
#'
#' @param blocking A logical value indicating whether to wait for completion.
#'
#' @return Invisible NULL.
#'
#' @examples
#' hi_reset_global_scheduler(TRUE)
#'
#' @export
hi_reset_global_scheduler <- function(blocking) {
  checkmate::assert_logical(blocking, len = 1, any.missing = FALSE)
  reset_global_scheduler(blocking)
}


#' Get Solver Information
#'
#' This function retrieves detailed information about the solver's state and performance.
#'
#' @param solver An object of class "highs_solver".
#'
#' @return A list containing solver information.
#'
#' @examples
#' solver <- example_solver()
#' info <- hi_solver_info(solver)
#'
#' @export
hi_solver_info <- function(solver) {
  checkmate::assert_class(solver, classes = "highs_solver")
  solver_info(solver)
}


#' Get Solution
#'
#' This function retrieves the solution from the solver after optimization.
#'
#' @param solver An object of class "highs_solver".
#'
#' @return A list containing the solution information.
#'
#' @examples
#' solver <- example_solver()
#' hi_solver_run(solver)
#' solution <- hi_solver_get_solution(solver)
#'
#' @export
hi_solver_get_solution <- function(solver) {
  checkmate::assert_class(solver, classes = "highs_solver")
  solver_get_solution(solver)
}


#' Set Solution for Warm-Start
#'
#' Sets a complete primal and/or dual solution on the solver for warm-starting.
#'
#' @param solver An object of class "highs_solver".
#' @param col_value Numeric vector of primal column values.
#' @param row_value Numeric vector of primal row values.
#' @param col_dual Numeric vector of dual column values.
#' @param row_dual Numeric vector of dual row values.
#' @param value_valid Logical, whether the primal values are valid.
#' @param dual_valid Logical, whether the dual values are valid.
#'
#' @return An integer status code (0 = OK).
#'
#' @export
hi_solver_set_solution <- function(solver, col_value, row_value,
                                   col_dual, row_dual,
                                   value_valid = TRUE, dual_valid = TRUE) {
  checkmate::assert_class(solver, classes = "highs_solver")
  solver_set_solution_obj(solver, value_valid, dual_valid,
                          col_value, col_dual, row_value, row_dual)
}


#' Set Sparse Solution for Warm-Start
#'
#' Sets a sparse primal solution (by index/value pairs) for warm-starting.
#'
#' @param solver An object of class "highs_solver".
#' @param index Integer vector of column indices (0-based).
#' @param value Numeric vector of values for those columns.
#'
#' @return An integer status code (0 = OK).
#'
#' @export
hi_solver_set_sparse_solution <- function(solver, index, value) {
  checkmate::assert_class(solver, classes = "highs_solver")
  checkmate::assert_integerish(index, any.missing = FALSE)
  checkmate::assert_numeric(value, len = length(index), any.missing = FALSE)
  solver_set_solution_vec(solver, as.integer(index), as.double(value))
}


#' Get the Dual Ray (Farkas Infeasibility Certificate)
#'
#' For a primal-infeasible LP, returns the dual unbounded ray: a Farkas
#' certificate of primal infeasibility, with one entry per constraint (row).
#'
#' A dual ray is available only after an LP detected to be **infeasible** and
#' solved by the **simplex** method (the interior-point solver does not produce
#' a ray). Disable presolve, or rely on HiGHS's postsolve, to recover the ray on
#' the original model. Requesting the ray may solve an auxiliary LP.
#'
#' @param solver An object of class "highs_solver".
#'
#' @return A list with `status` (integer, 0 = OK), `has_dual_ray` (logical), and
#'   `dual_ray` (numeric vector of length `n_row`, or `NULL` when no ray exists).
#'
#' @export
hi_solver_get_dual_ray <- function(solver) {
  checkmate::assert_class(solver, classes = "highs_solver")
  solver_get_dual_ray(solver)
}


#' Get the Primal Ray (Unboundedness Certificate)
#'
#' For a primal-unbounded LP (i.e. dual-infeasible), returns the primal unbounded
#' ray: a certificate of unboundedness, with one entry per variable (column).
#'
#' As with the dual ray, this is a simplex-method certificate; requesting it may
#' solve an auxiliary LP.
#'
#' @param solver An object of class "highs_solver".
#'
#' @return A list with `status` (integer, 0 = OK), `has_primal_ray` (logical), and
#'   `primal_ray` (numeric vector of length `n_col`, or `NULL` when no ray exists).
#'
#' @export
hi_solver_get_primal_ray <- function(solver) {
  checkmate::assert_class(solver, classes = "highs_solver")
  solver_get_primal_ray(solver)
}


#' Get Basis
#'
#' Retrieves the current simplex basis from the solver.
#'
#' @param solver An object of class "highs_solver".
#'
#' @return A list with components:
#' \describe{
#'   \item{valid}{Logical, whether the basis is valid.}
#'   \item{col_status}{Integer vector of column basis statuses
#'     (0=Lower, 1=Basic, 2=Upper, 3=Zero, 4=Nonbasic).}
#'   \item{row_status}{Integer vector of row basis statuses.}
#' }
#'
#' @export
hi_solver_get_basis <- function(solver) {
  checkmate::assert_class(solver, classes = "highs_solver")
  solver_get_basis(solver)
}


#' Set Basis for Warm-Start
#'
#' Sets a simplex basis on the solver for warm-starting. Use basis
#' status values: 0=Lower, 1=Basic, 2=Upper, 3=Zero, 4=Nonbasic.
#'
#' @param solver An object of class "highs_solver".
#' @param col_status Integer vector of column basis statuses.
#' @param row_status Integer vector of row basis statuses.
#'
#' @return An integer status code (0 = OK).
#'
#' @export
hi_solver_set_basis <- function(solver, col_status, row_status) {
  checkmate::assert_class(solver, classes = "highs_solver")
  checkmate::assert_integerish(col_status, any.missing = FALSE)
  checkmate::assert_integerish(row_status, any.missing = FALSE)
  solver_set_basis(solver, as.integer(col_status), as.integer(row_status))
}


#' Clear Basis
#'
#' Invalidates the solver's current basis, allowing presolve to run
#' on the next call to \code{hi_solver_run}.
#'
#' @param solver An object of class "highs_solver".
#'
#' @return An integer status code (0 = OK).
#'
#' @export
hi_solver_clear_basis <- function(solver) {
  checkmate::assert_class(solver, classes = "highs_solver")
  solver_clear_basis(solver)
}


#' Get Ranging (Sensitivity Analysis)
#'
#' Computes ranging information for the solved LP. The solver must have
#' an optimal LP solution before calling this function.
#'
#' @param solver An object of class "highs_solver".
#'
#' @return A list with components:
#' \describe{
#'   \item{valid}{Logical, whether the ranging data is valid.}
#'   \item{col_cost_up, col_cost_dn}{Lists with value, objective, in_var, ou_var
#'     vectors for cost ranging (length num_col).}
#'   \item{col_bound_up, col_bound_dn}{Lists for column bound ranging.}
#'   \item{row_bound_up, row_bound_dn}{Lists for row bound ranging (length num_row).}
#' }
#'
#' @export
hi_solver_get_ranging <- function(solver) {
  checkmate::assert_class(solver, classes = "highs_solver")
  solver_get_ranging(solver)
}


#' Get HiGHS Version
#'
#' Returns the version string of the bundled HiGHS solver.
#'
#' @param solver An object of class "highs_solver".
#'
#' @return A character string with the HiGHS version.
#'
#' @export
hi_solver_version <- function(solver) {
  checkmate::assert_class(solver, classes = "highs_solver")
  solver_version(solver)
}


#' Get Run Time
#'
#' Returns the elapsed time of the last solver run.
#'
#' @param solver An object of class "highs_solver".
#'
#' @return A numeric value (seconds).
#'
#' @export
hi_solver_get_run_time <- function(solver) {
  checkmate::assert_class(solver, classes = "highs_solver")
  solver_get_run_time(solver)
}


#' Presolve
#'
#' Runs presolve on the current model without solving.
#'
#' @param solver An object of class "highs_solver".
#'
#' @return An integer status code (0 = OK).
#'
#' @export
hi_solver_presolve <- function(solver) {
  checkmate::assert_class(solver, classes = "highs_solver")
  solver_presolve(solver)
}


#' Postsolve
#'
#' Performs postsolve using a given solution from the presolved model.
#'
#' @param solver An object of class "highs_solver".
#' @param col_value Numeric vector of column values.
#' @param col_dual Numeric vector of column duals.
#' @param row_value Numeric vector of row values.
#' @param row_dual Numeric vector of row duals.
#'
#' @return An integer status code (0 = OK).
#'
#' @export
hi_solver_postsolve <- function(solver, col_value, col_dual, row_value, row_dual) {
  checkmate::assert_class(solver, classes = "highs_solver")
  solver_postsolve_solution(solver, col_value, col_dual, row_value, row_dual)
}


#' Read Model from File
#'
#' Reads a model from a file (MPS or LP format).
#'
#' @param solver An object of class "highs_solver".
#' @param filename Path to the model file.
#'
#' @return An integer status code (0 = OK).
#'
#' @export
hi_solver_read_model <- function(solver, filename) {
  checkmate::assert_class(solver, classes = "highs_solver")
  checkmate::assert_file_exists(filename)
  solver_read_model(solver, filename)
}


#' Get Boolean Option Value
#'
#' This function retrieves the value of a boolean solver option.
#'
#' @param solver An object of class "highs_solver".
#' @param key A character string specifying the option name.
#'
#' @return A logical value.
#'
#' @examples
#' solver <- example_solver()
#' value <- hi_solver_get_bool_option(solver, "mip_detect_symmetry")
#'
#' @export
hi_solver_get_bool_option <- function(solver, key) {
  checkmate::assert_class(solver, classes = "highs_solver")
  checkmate::assert_character(key, len = 1, any.missing = FALSE)
  solver_get_bool_option(solver, key)
}


#' Get Integer Option Value
#'
#' This function retrieves the value of an integer solver option.
#'
#' @param solver An object of class "highs_solver".
#' @param key A character string specifying the option name.
#'
#' @return An integer value.
#'
#' @examples
#' solver <- example_solver()
#' value <- hi_solver_get_int_option(solver, "log_dev_level")
#'
#' @export
hi_solver_get_int_option <- function(solver, key) {
  checkmate::assert_class(solver, classes = "highs_solver")
  checkmate::assert_character(key, len = 1, any.missing = FALSE)
  solver_get_int_option(solver, key)
}


#' Get Double Option Value
#'
#' This function retrieves the value of a double precision solver option.
#'
#' @param solver An object of class "highs_solver".
#' @param key A character string specifying the option name.
#'
#' @return A numeric value.
#'
#' @examples
#' solver <- example_solver()
#' value <- hi_solver_get_dbl_option(solver, "time_limit")
#'
#' @export
hi_solver_get_dbl_option <- function(solver, key) {
  checkmate::assert_class(solver, classes = "highs_solver")
  checkmate::assert_character(key, len = 1, any.missing = FALSE)
  solver_get_dbl_option(solver, key)
}


#' Get String Option Value
#'
#' This function retrieves the value of a string solver option.
#'
#' @param solver An object of class "highs_solver".
#' @param key A character string specifying the option name.
#'
#' @return A character string.
#'
#' @examples
#' solver <- example_solver()
#' value <- hi_solver_get_str_option(solver, "solver")
#'
#' @export
hi_solver_get_str_option <- function(solver, key) {
  checkmate::assert_class(solver, classes = "highs_solver")
  checkmate::assert_character(key, len = 1, any.missing = FALSE)
  solver_get_str_option(solver, key)
}


#' Change Variable Bounds
#'
#' This function updates the bounds of an existing variable in the model.
#'
#' @param solver An object of class "highs_solver".
#' @param idx An integer specifying the variable index.
#' @param lower The new lower bound.
#' @param upper The new upper bound.
#'
#' @return The solver instance with updated bounds.
#'
#' @examples
#' solver <- example_solver()
#' hi_solver_change_variable_bounds(solver, 1, 0, 10)
#'
#' @export
hi_solver_change_variable_bounds <- function(solver, idx, lower, upper) {
  checkmate::assert_class(solver, classes = "highs_solver")
  checkmate::assert_integerish(
    idx,
    lower = 0,
    any.missing = FALSE,
    min.len = 1
  )
  checkmate::assert_numeric(lower, len = length(idx), any.missing = FALSE)
  checkmate::assert_numeric(upper, len = length(idx), any.missing = FALSE)
  solver_change_variable_bounds(
    solver,
    as.integer(idx),
    as.double(lower),
    as.double(upper)
  )
}


#' Change Constraint Bounds
#'
#' This function updates the bounds of an existing constraint in the model.
#'
#' @param solver An object of class "highs_solver".
#' @param idx An integer vector specifying the constraint indices. Note that idx starts at 0.
#' @param lhs The new left-hand side bound.
#' @param rhs The new right-hand side bound.
#'
#' @details
#' This is a low-level wrapper for `highs.changeRowsBounds` that only does basic
#' checks to prevent segfaults, but otherwise behaves exactly like
#' `highs.changeRowsBounds`.
#'
#' @return The solver instance with updated constraint bounds.
#'
#' @examples
#' solver <- example_solver()
#' hi_solver_change_constraint_bounds(solver, 1, -Inf, 100)
#'
#' @export
hi_solver_change_constraint_bounds <- function(solver, idx, lhs, rhs) {
  checkmate::assert_class(solver, classes = "highs_solver")
  checkmate::assert_integerish(
    idx,
    lower = 0,
    any.missing = FALSE,
    min.len = 1
  )
  checkmate::assert_numeric(lhs, len = length(idx), any.missing = FALSE)
  checkmate::assert_numeric(rhs, len = length(idx), any.missing = FALSE)
  solver_change_constraint_bounds(
    solver,
    as.integer(idx),
    as.double(lhs),
    as.double(rhs)
  )
}


#' Add Constraints to Model
#'
#' This function adds new constraints (rows) to the optimization model.
#'
#' @param solver An object of class "highs_solver".
#' @param lhs A numeric vector of left-hand side bounds.
#' @param rhs A numeric vector of right-hand side bounds.
#' @param start An integer vector of starting positions in the sparse matrix.
#' @param index An integer vector of column indices.
#' @param value A numeric vector of coefficient values.
#'
#' @return The solver instance with new constraints added.
#'
#' @examples
#' solver <- example_solver()
#' hi_solver_add_rows(solver, c(-Inf), c(10), c(0, 2), c(0, 1), c(1, 2))
#'
#' @export
hi_solver_add_rows <- function(solver, lhs, rhs, start, index, value) {
  checkmate::assert_class(solver, classes = "highs_solver")
  checkmate::assert_numeric(lhs, min.len = 1, any.missing = FALSE)
  checkmate::assert_numeric(rhs, len = length(lhs), any.missing = FALSE)
  checkmate::assert_integerish(
    start,
    lower = 0,
    any.missing = FALSE,
    min.len = 1
  )
  checkmate::assert_integerish(
    index,
    lower = 0,
    any.missing = FALSE,
    min.len = 1
  )
  checkmate::assert_numeric(
    value,
    len = length(index),
    any.missing = FALSE,
    min.len = 1
  )
  solver_add_rows(
    solver,
    as.double(lhs),
    as.double(rhs),
    as.integer(start),
    as.integer(index),
    as.double(value)
  )
}


#' Add Variables to Model
#'
#' This function adds new variables (columns) to the optimization model.
#'
#' @param solver An object of class "highs_solver".
#' @param costs A numeric vector of objective coefficients.
#' @param lower A numeric vector giving the lower bounds of the new variables.
#' @param upper A numeric vector giving the upper bounds of the new variables.
#' @param start An integer vector of starting positions in the sparse matrix.
#' @param index An integer vector of row indices.
#' @param value A numeric vector of coefficient values.
#'
#' @return The solver instance with new variables added.
#'
#' @examples
#' solver <- example_solver()
#' hi_solver_add_cols(solver, c(1), c(0), c(10), c(0, 1), c(0), c(2))
#'
#' @export
hi_solver_add_cols <- function(
  solver,
  costs,
  lower,
  upper,
  start,
  index,
  value
) {
  checkmate::assert_class(solver, classes = "highs_solver")
  checkmate::assert_numeric(costs, min.len = 1, any.missing = FALSE)
  checkmate::assert_numeric(lower, len = length(costs), any.missing = FALSE)
  checkmate::assert_numeric(upper, len = length(costs), any.missing = FALSE)
  checkmate::assert_integerish(
    start,
    lower = 0,
    any.missing = FALSE,
    min.len = 1
  )
  checkmate::assert_integerish(
    index,
    lower = 0,
    any.missing = FALSE,
    min.len = 1
  )
  checkmate::assert_numeric(
    value,
    len = length(index),
    any.missing = FALSE,
    min.len = 1
  )
  solver_add_cols(
    solver,
    as.double(costs),
    as.double(lower),
    as.double(upper),
    as.integer(start),
    as.integer(index),
    as.double(value)
  )
}


#' Get Objective Coefficients
#'
#' This function retrieves the objective coefficients of the linear program.
#'
#' @param solver An object of class "highs_solver".
#'
#' @return A numeric vector of objective coefficients.
#'
#' @examples
#' solver <- example_solver()
#' costs <- hi_solver_get_lp_costs(solver)
#'
#' @export
hi_solver_get_lp_costs <- function(solver) {
  checkmate::assert_class(solver, classes = "highs_solver")
  solver_get_lp_costs(solver)
}


#' Get Variable Bounds
#'
#' This function retrieves the bounds for all variables in the model.
#'
#' @param solver An object of class "highs_solver".
#'
#' @return A list containing lower and upper bounds for all variables.
#'
#' @examples
#' solver <- example_solver()
#' bounds <- hi_solver_get_variable_bounds(solver)
#'
#' @export
hi_solver_get_variable_bounds <- function(solver) {
  checkmate::assert_class(solver, classes = "highs_solver")
  solver_get_variable_bounds(solver)
}


#' Get Constraint Bounds
#'
#' This function retrieves the bounds for all constraints in the model.
#'
#' @param solver An object of class "highs_solver".
#'
#' @return A list containing lower and upper bounds for all constraints.
#'
#' @examples
#' solver <- example_solver()
#' bounds <- hi_solver_get_constraint_bounds(solver)
#'
#' @export
hi_solver_get_constraint_bounds <- function(solver) {
  checkmate::assert_class(solver, classes = "highs_solver")
  solver_get_constraint_bounds(solver)
}


#' Get Constraint Matrix
#'
#' This function retrieves the constraint matrix of the optimization model.
#'
#' @param solver An object of class "highs_solver".
#'
#' @return A sparse matrix representing the constraints.
#'
#' @examples
#' solver <- example_solver()
#' matrix <- hi_solver_get_constraint_matrix(solver)
#'
#' @export
hi_solver_get_constraint_matrix <- function(solver) {
  checkmate::assert_class(solver, classes = "highs_solver")
  solver_get_constraint_matrix(solver)
}


#' Get Variable Types
#'
#' This function retrieves the type (continuous, integer, etc.) of all variables.
#'
#' @param solver An object of class "highs_solver".
#'
#' @return A character vector of variable types.
#'
#' @examples
#' solver <- example_solver()
#' types <- hi_solver_get_vartype(solver)
#'
#' @export
hi_solver_get_vartype <- function(solver) {
  checkmate::assert_class(solver, classes = "highs_solver")
  solver_get_vartype(solver)
}
