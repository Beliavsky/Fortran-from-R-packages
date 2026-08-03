#' @import checkmate
#' @importFrom Rcpp evalCpp
#' @importFrom stats setNames
#' @importFrom utils modifyList capture.output tail str
#' @useDynLib highs, .registration=TRUE
NULL


highs_globals <- local({
  globals <- list(threads = NA_integer_)
  function(key, value) {
    if (missing(key)) {
      return(globals)
    }
    if (missing(value)) {
      globals[[key]]
    } else {
      globals[[key]] <<- value
    }
  }
})


highs_infinity <- function() {
  if (is.null(highs_globals("Inf"))) {
    highs_globals("Inf", solver_infinity())
  }
  highs_globals("Inf")
}


highs_variable_types <- function() {
  c("C", "I", "SC", "SI", "II")
}


get_number_of_threads <- function() {
  highs_globals("threads")
}


set_number_of_threads <- function(threads) {
  if (!isTRUE(highs_globals("threads") == threads)) {
    reset_global_scheduler(TRUE)
    highs_globals("threads", threads)
  }
}


csc_to_matrix <- function(
  start,
  index,
  value,
  nrow = max(index + 1L),
  ncol = length(start) - 1L
) {
  stopifnot(length(index) == length(value))
  ind <- index + 1L
  M <- matrix(0, nrow, ncol)
  for (i in seq_along(index)) {
    row_id <- ind[i]
    col_id <- min(which(start >= i) - 1L)
    M[row_id, col_id] <- value[i]
  }
  M
}


model_set_hessian <- function(model, x) {
  if (is_dgc(x)) {
    model_set_hessian_(
      model,
      "square",
      dim = x@Dim[1L],
      start = x@p,
      index = x@i,
      value = x@x
    )
  } else if (is_csc(x)) {
    model_set_hessian_(
      model,
      "square",
      dim = x@dimension[1L],
      start = x@ia - 1L,
      index = x@ja - 1L,
      value = x@ra
    )
  } else if (is_dgr(x)) {
    stop("dgRMatrix is not supported for the hessian")
  } else if (is_csr(x)) {
    stop("matrix.csr is not supported for the hessian")
  } else if (is_stm(x)) {
    ind <- order(x$j, x$i)
    start <- c(0L, cumsum(tabulate(x$j[ind], x$ncol)))
    model_set_hessian_(
      model,
      "square",
      dim = x[["nrow"]],
      start = start,
      index = x$i[ind] - 1L,
      value = x$v[ind]
    )
  } else if (is_dmat(x)) {
    ind <- which(x != 0, arr.ind = TRUE)
    start <- c(0L, cumsum(tabulate(ind[, 2L], NCOL(x))))
    model_set_hessian_(
      model,
      "square",
      dim = NROW(x),
      start = start,
      index = ind[, 1] - 1L,
      value = x[ind]
    )
  } else {
    stop(sprintf("unkown class %s", deparse(class(x))))
  }
}


model_set_constraint_matrix <- function(model, x) {
  if (is_dgc(x)) {
    model_set_constraint_matrix_(
      model,
      "colwise",
      start = x@p,
      index = x@i,
      value = x@x
    )
  } else if (is_csc(x)) {
    model_set_constraint_matrix_(
      model,
      "colwise",
      start = x@ia - 1L,
      index = x@ja - 1L,
      value = x@ra
    )
  } else if (is_dgr(x)) {
    model_set_constraint_matrix_(
      model,
      "rowwise",
      start = x@p,
      index = x@j,
      value = x@x
    )
  } else if (is_csr(x)) {
    model_set_constraint_matrix_(
      model,
      "rowwise",
      start = x@ia - 1L,
      index = x@ja - 1L,
      value = x@ra
    )
  } else if (is_stm(x)) {
    ind <- order(x$j, x$i)
    start <- c(0L, cumsum(tabulate(x$j[ind], x$ncol)))
    model_set_constraint_matrix_(
      model,
      "colwise",
      start = start,
      index = x$i[ind] - 1L,
      value = x$v[ind]
    )
  } else if (is_dmat(x)) {
    ind <- which(x != 0, arr.ind = TRUE)
    start <- c(0L, cumsum(tabulate(ind[, 2L], NCOL(x))))
    model_set_constraint_matrix_(
      model,
      "colwise",
      start = start,
      index = ind[, 1] - 1L,
      value = x[ind]
    )
  } else {
    stop(sprintf("unkown class %s", deparse(class(x))))
  }
}


number_or_rows <- function(x) {
  if (is_dgc(x) || is_dgr(x)) {
    x@Dim[1L]
  } else if (is_csc(x) || is_csr(x)) {
    x@dimension[1L]
  } else if (is_stm(x)) {
    x[["nrow"]]
  } else if (is_dmat(x)) {
    NROW(x)
  } else {
    stop(sprintf("unkown class %s", deparse(class(x))))
  }
}


allowed_matrix_classes <- function() {
  c(
    "matrix",
    "simple_triplet_matrix",
    "dgCMatrix",
    "matrix.csc",
    "dgRMatrix",
    "matrix.csr"
  )
}


#' Create a Highs Model
#'
#' Solve linear and quadratic mixed integer optimization problems.
#'
#' @param Q a numeric symmetric matrix giving the quadratic part of the objective.
#' @param L a numeric vector giving the linear part of the objective function.
#' @param lower a numeric vector giving the lower bounds of the variables.
#' @param upper a numeric vector giving the upper bounds of the variables.
#' @param A a numeric matrix giving the linear part of the constraints. Rows are
#'   constraints, and columns are decision variables.
#' @param lhs a numeric vector giving the left hand-side of the linear constraints.
#' @param rhs a numeric vector giving the right hand-side of the linear constraints.
#' @param types a integer vector or character vector giving the variable types.
#'      \code{'C'} or \code{'1'} for continuous,
#'      \code{'I'} or \code{'2'} for integer,
#'      \code{'SC'} or \code{'3'} for semi continuous,
#'      \code{'SI'} or \code{'4'} for semi integer and
#'      \code{'II'} or \code{'5'} for implicit integer.
#' @param maximum a logical if \code{TRUE} the solver searches for a maximum,
#'                if \code{FALSE} the solver searches for a minimum.
#' @param offset a numeric value giving the offset (default is \code{0}).
#'
#' @return A an object of class \code{highs_model}.
#'
#' @examples
#' library("highs")
#' # Minimize:
#' #  x_0 +  x_1 + 3
#' # Subject to:
#' #               x_1 <=  7
#' #  5 <=  x_0 + 2x_1 <= 15
#' #  6 <= 3x_0 + 2x_1
#' #  0 <= x_0 <= 4
#' #  1 <= x_1
#' A <- rbind(c(0, 1), c(1, 2), c(3, 2))
#' m <- highs_model(L = c(1.0, 1), lower = c(0, 1), upper = c(4, Inf),
#'                  A = A, lhs = c(-Inf, 5, 6), rhs = c(7, 15, Inf),
#'                  offset = 3)
#' m
#'
#' # Minimize:
#' #  -x_2 - 3x_3 + (1/2) * (2 x_1^2 - 2 x_1x_3 + 0.2 x_2^2 + 2 x_3^2)
#' # Subject to:
#' #  x_1 + x_3 <= 2
#' #  0 <= x
#' L <- c(0, -1, -3)
#' Q <- rbind(c(2, 0.0, -1), c(0, 0.2, 0), c(-1, 0.0, 2))
#' A <- cbind(1, 0, 1)
#' m <- highs_model(Q = Q, L = L, lower = 0, A = A, rhs = 2)
#' m
#' @export
highs_model <- function(
  Q = NULL,
  L,
  lower,
  upper,
  A = NULL,
  lhs = NULL,
  rhs = NULL,
  types = rep.int(1L, length(L)),
  maximum = FALSE,
  offset = 0
) {
  assert_numeric(L, any.missing = FALSE)
  assert_multi_class(Q, allowed_matrix_classes(), null.ok = TRUE)
  assert_multi_class(A, allowed_matrix_classes(), null.ok = TRUE)

  nvars <- length(L)

  INF <- highs_infinity()
  if (missing(lower) || length(lower) == 0L) {
    lower <- rep.int(-INF, nvars)
  } else if (length(lower) == 1L) {
    lower <- rep.int(lower, nvars)
  }
  assert_numeric(lower, any.missing = FALSE, len = nvars)
  if (missing(upper) || length(upper) == 0L) {
    upper <- rep.int(INF, nvars)
  } else if (length(upper) == 1L) {
    upper <- rep.int(upper, nvars)
  }
  assert_numeric(upper, any.missing = FALSE, len = nvars)

  if (NROW(A) == 0L) {
    A <- lhs <- rhs <- NULL
    ncons <- 0L
  } else {
    ncons <- number_or_rows(A)
  }

  if (ncons > 0L) {
    if (length(lhs) == 0L) {
      lhs <- rep.int(-INF, ncons)
    } else {
      lhs <- replace(lhs, lhs == -Inf, -INF)
    }
    assert_numeric(lhs, any.missing = FALSE, len = ncons)
    if (length(rhs) == 0L) {
      rhs <- rep.int(INF, ncons)
    } else {
      rhs <- replace(rhs, rhs == Inf, INF)
    }
    assert_numeric(rhs, any.missing = FALSE, len = ncons)
  }

  all_types_continuous <- NA
  if (length(types) == 0L) {
    types <- rep.int(0L, nvars)
    all_types_continuous <- TRUE
  } else {
    if (is.character(types)) {
      types <- match(types, highs_variable_types()) - 1L
    } else {
      types <- types - 1L
    }
    assert_integerish(
      types,
      lower = 0,
      upper = 4L,
      any.missing = FALSE,
      len = nvars
    )
  }

  model <- new_model()
  model_set_ncol(model, nvars)
  model_set_nrow(model, ncons)
  model_set_sense(model, maximum)
  model_set_objective(model, L)
  if (!is.null(Q)) {
    model_set_hessian(model, Q)
  }

  if (!isTRUE(all_types_continuous) || any(types != 0)) {
    model_set_vartype(model, as.integer(types))
  }

  lower <- replace(lower, lower == -Inf, -INF)
  upper <- replace(upper, upper == Inf, INF)
  model_set_lower(model, lower)
  model_set_upper(model, upper)
  if (ncons > 0L) {
    model_set_constraint_matrix(model, A)
    model_set_lhs(model, lhs)
    model_set_rhs(model, rhs)
  }
  if (offset != 0) {
    model_set_offset(model, offset)
  }
  class(model) <- c("highs_model", class(model))
  return(model)
}


#' Write a Highs Model to a File
#'
#' Write an highs model to file.
#'
#' @param model an object of class \code{highs_model}.
#' @param file a character string giving the filename.
#'
#' @examples
#' model <- example_model()
#' model_file <- tempfile(fileext = ".mps")
#' highs_write_model(model, model_file)
#'
#' @export
highs_write_model <- function(model, file) {
  checkmate::assert_string(file)
  checkmate::assert_class(model, classes = "highs_model")
  checkmate::assert_directory_exists(dirname(file), access = "w")
  solver <- hi_new_solver(model)
  solver_write_model(solver, file)
}


#' Solve an Optimization Problems
#'
#' Solve linear and quadratic mixed integer optimization problems.
#'
#' @param Q a numeric symmetric matrix giving the quadratic part of the objective.
#' @param L a numeric vector giving the linear part of the objective function.
#' @param lower a numeric vector giving the lower bounds of the variables.
#' @param upper a numeric vector giving the upper bounds of the variables.
#' @param A a numeric matrix giving the linear part of the constraints. Rows are
#'   constraints, and columns are decision variables.
#' @param lhs a numeric vector giving the left hand-side of the linear constraints.
#' @param rhs a numeric vector giving the right hand-side of the linear constraints.
#' @param types a integer vector or character vector giving the variable types.
#'      \code{'C'} or \code{'1'} for continuous,
#'      \code{'I'} or \code{'2'} for integer,
#'      \code{'SC'} or \code{'3'} for semi continuous,
#'      \code{'SI'} or \code{'4'} for semi integer and
#'      \code{'II'} or \code{'5'} for implicit integer.
#' @param maximum a logical if \code{TRUE} the solver searches for a maximum,
#'                if \code{FALSE} the solver searches for a minimum.
#' @param offset a numeric value giving the offset (default is \code{0}).
#' @param start a numeric vector giving the start solution.
#' Missing values can be used to specify that particular
#' variables should have their starting values calculated automatically.
#' @param control a list giving additional options for the solver,
#'                see \link{highs_available_solver_options} or the \code{README} file
#'                for a list of all available options.
#'
#' @return A \code{list} containing the result provided by the solver,
#'  containing the following named objects:
#' \item{\code{primal_solution}}{a numeric vector giving the primal solution.}
#' \item{\code{objective_value}}{a numeric giving the objective value.}
#' \item{\code{status}}{an integer giving the status code}
#' \item{\code{status_message}}{a character string giving the status message (explanation of the \code{status_code}).}
#' \item{\code{solver_msg}}{a list giving the original (not canonicalized) solver message.}
#' \item{\code{info}}{a list giving additional information provided by the solver.}
#'
#' Additional information on can be found in the \code{README} file.
#'
#' @examples
#' library("highs")
#' # Minimize:
#' #  x_0 +  x_1 + 3
#' # Subject to:
#' #               x_1 <=  7
#' #  5 <=  x_0 + 2x_1 <= 15
#' #  6 <= 3x_0 + 2x_1
#' #  0 <= x_0 <= 4
#' #  1 <= x_1
#' A <- rbind(c(0, 1), c(1, 2), c(3, 2))
#' s <- highs_solve(L = c(1.0, 1), lower = c(0, 1), upper = c(4, Inf),
#'                  A = A, lhs = c(-Inf, 5, 6), rhs = c(7, 15, Inf),
#'                  offset = 3)
#' s[["objective_value"]]
#' s[["primal_solution"]]
#'
#' # Minimize:
#' #  -x_2 - 3x_3 + (1/2) * (2 x_1^2 - 2 x_1x_3 + 0.2 x_2^2 + 2 x_3^2)
#' # Subject to:
#' #  x_1 + x_3 <= 2
#' #  0 <= x
#' L <- c(0, -1, -3)
#' Q <- rbind(c(2, 0.0, -1), c(0, 0.2, 0), c(-1, 0.0, 2))
#' A <- cbind(1, 0, 1)
#' s <- highs_solve(Q = Q, L = L, lower = 0, A = A, rhs = 2)
#' s[["objective_value"]]
#' s[["primal_solution"]]
#' @export
highs_solve <- function(
  Q = NULL,
  L,
  lower,
  upper,
  A = NULL,
  lhs = NULL,
  rhs = NULL,
  types = rep.int(1L, length(L)),
  maximum = FALSE,
  offset = 0,
  start = NULL,
  control = highs_control()
) {
  checkmate::assert_list(control)
  if (!inherits(control, "highs_control")) {
    control <- do.call(highs_control, control)
  }
  verbose <- control$log_to_console %||% FALSE
  model <- highs_model(
    Q = Q,
    L = L,
    lower = lower,
    upper = upper,
    A = A,
    lhs = lhs,
    rhs = rhs,
    types = types,
    maximum = maximum,
    offset = offset
  )

  set_number_of_threads(control$threads)
  solver <- hi_new_solver(model)
  hi_solver_set_options(solver, control)
  if (!is.null(start)) {
    hi_solver_set_start(solver, start)
  }

  run_status <- hi_solver_run(solver, verbose = verbose)
  status <- solver_status(solver)
  status_message <- solver_status_message(solver)

  solution <- solver_get_solution(solver)
  info <- solver_info(solver)
  list(
    primal_solution = solution[["col_value"]],
    objective_value = info[["objective_function_value"]],
    status = status,
    status_message = status_message,
    solver_msg = solution,
    info = info
  )
}


#' Highs Control
#'
#' @param threads an integer giving the number of threads to be used.
#' @param time_limit a double giving the time limit.
#' @param log_to_console a logical giving if the output should be shown in
#'                       the console.
#' @param ... other arguments supported by the \code{HiGHS} solver.
#'
#' @examples
#' control <- highs_control()
#' @export
highs_control <- function(
  threads = 1L,
  time_limit = Inf,
  log_to_console = FALSE,
  ...
) {
  checkmate::assert_integerish(
    threads,
    len = 1L,
    lower = 1L,
    any.missing = FALSE
  )
  checkmate::assert_double(time_limit, len = 1L, any.missing = FALSE)
  checkmate::assert_logical(log_to_console, len = 1L, any.missing = FALSE)
  control <- c(as.list(environment()), list(...))
  default_control <- list(parallel = "off", solver = "choose")
  control <- modifyList(default_control, control)
  if (is.infinite(control[["time_limit"]])) {
    control[["time_limit"]] <- NULL
  }
  if (isTRUE(control$solver == "pdlp")) {
    stop("solver 'pdlp' is not supported in the R highs interface.")
  }
  control$parallel <- if (isTRUE(threads > 1)) "on" else "off"
  structure(control, class = "highs_control")
}


#' @noRd
#' @export
print.highs_control <- function(x, ...) {
  writeLines("highs_control")
  out <- capture.output(str(unclass(x), ...))
  writeLines(out[-1L])
}


solver_get_types <- function(solver) {
  vtypes <- solver_get_vartype(solver)
  if (length(vtypes) == 0L) {
    rep("C", solver_get_num_col(solver))
  } else {
    highs_variable_types()[as.integer(vtypes + 1L)]
  }
}


solver_get_amatrix <- function(solver) {
  mat <- solver_get_constraint_matrix(solver)
  mat_formats <- c("colwise", "rowwise", "rowwise_partitioned")
  mat[["format"]] <- mat_formats[mat[["format"]]]
  structure(mat, class = "highs_sparse_matrix")
}


#' @noRd
#' @export
print.highs_control <- function(x, ...) {
  writeLines(sprintf(
    "A %ix%i highs sparse matrix.",
    x[["nrow"]],
    x[["ncol"]]
  ))
}


#' Highs Solver
#'
#' @description
#' Create a wrapper around the \code{HiGHS} solver. Manly usefull if one wants
#' a low level wrapper around highs with hot-start capabilities. \cr
#'
#' @param model an object of class \code{"highs_model"} created with \code{highs_model()}.
#' @param control an object of class \code{"highs_control"} created with \code{highs_control()}.
#'
#' @details
#' \strong{Methods} \cr
#' The following methods are provided by the \code{"highs_solver"} class.
#' \itemize{
#'    \item \code{solve(...)} method to be called to solve the optimization problem.
#'          Returns an integer giving the status code returned by \strong{HiGHS}.
#'    \item \code{status()} method to obtain the status from the solver.
#'    \item \code{status_message()} method to obtain the status message from the solver.
#'    \item \code{solution()} method to obtain the solution from the solver.
#'    \item \code{info()} info to obtain addtional information from the solver.
#'    \item \code{L(i, v)} method to get and set the linear part of the objective.
#'    \item \code{A(i, j, v)} method to get and set the constraint matrix coefficients.
#'    \item \code{cbounds(i, lhs, rhs)} method to get and set the constraint bounds
#'                                      (left hand-side and right hand-side).
#'    \item \code{types(i, v)} method to get and set the variable types.
#'    \item \code{vbounds(i, lower, upper)} method to get and set the variable bounds.
#'    \item \code{maximum(maximize)} method to get and set the sense of the problem.
#' }
#' \strong{Method arguments} \cr
#' \itemize{
#'    \item \code{...} optional control arguments, which can be used to alter the options
#'                     set via the \code{control} argument when initializing the solver.
#'    \item \code{i} a vector of integers giving the index (vector index or row index)
#'                   of the coeficcients to be altered.
#'    \item \code{j} a vector of integers giving the index (column index) of the coeficcients to be altered.
#'    \item \code{v} a vector of doubles giving the values of the coeficcients to be altered.
#'    \item \code{lhs} a vector of doubles giving left hand-side.
#'    \item \code{rhs} a vector of doubles giving right hand-side.
#'    \item \code{lower} a vector of doubles giving the lower bounds to be altered.
#'    \item \code{upper} a vector of doubles giving the upper bounds to be altered.
#' }
#'
#' @return an object of class \code{"highs_solver"}.
#'
#' @examples
#' A <- rbind(c(0, 1), c(1, 2), c(3, 2))
#' m <- highs_model(L = c(1.0, 1), lower = c(0, 1), upper = c(4, Inf),
#'                  A = A, lhs = c(-Inf, 5, 6), rhs = c(7, 15, Inf),
#'                  offset = 3)
#' solver <- highs_solver(m)
#'
#' @export
highs_solver <- function(model, control = highs_control()) {
  checkmate::assert_class(model, classes = "highs_model")
  checkmate::assert_class(control, classes = "highs_control")
  set_number_of_threads(control$threads)
  solver <- hi_new_solver(model)
  hi_solver_set_options(solver, control)
  solve <- function(...) {
    cntrl <- list(...)
    if (length(cntrl) == 0L) {
      hi_solver_get_options(solver)
    } else {
      hi_solver_set_options(solver, cntrl)
    }
    hi_solver_run(solver)
  }
  status <- function() {
    solver_status(solver)
  }
  status_message <- function() {
    solver_status_message(solver)
  }
  solution <- function() {
    solver_get_solution(solver)
  }
  info <- function() {
    solver_info(solver)
  }
  L <- function(i, v) {
    if (missing(i) && missing(v)) {
      solver_get_lp_costs(solver)
    } else {
      assert_integerish(i, lower = 1, any.missing = FALSE)
      assert_numeric(v, any.missing = FALSE, len = length(i))
      solver_set_objective(solver, as.integer(i) - 1L, v)
    }
  }
  A <- function(i, j, v) {
    if (missing(i) && missing(j) && missing(v)) {
      solver_get_amatrix(solver)
    } else {
      assert_integerish(i, lower = 1, any.missing = FALSE)
      assert_integerish(
        j,
        lower = 1,
        any.missing = FALSE,
        len = length(i)
      )
      assert_numeric(v, any.missing = FALSE, len = length(i))
      solver_set_coeff(solver, as.integer(i) - 1L, as.integer(j) - 1L, v)
    }
  }
  cbounds <- function(i, lhs, rhs) {
    if (missing(i) && missing(lhs) && missing(rhs)) {
      cbnd <- solver_get_constraint_bounds(solver)
      matrix(cbnd, ncol = 2L, dimnames = list(NULL, c("lhs", "rhs")))
    } else {
      assert_integerish(i, lower = 1, any.missing = FALSE)
      assert_numeric(lhs, any.missing = FALSE, len = length(i))
      assert_numeric(rhs, any.missing = FALSE, len = length(i))
      solver_set_constraint_bounds(solver, as.integer(i) - 1L, lhs, rhs)
    }
  }
  types <- function(i, v) {
    if (missing(i) && missing(v)) {
      solver_get_types(solver)
    } else {
      assert_integerish(i, lower = 1, any.missing = FALSE)
      if (is.character(v)) {
        v <- match(v, highs_variable_types()) - 1L
      } else {
        v <- v - 1L
      }
      assert_integerish(
        v,
        lower = 0,
        upper = 4L,
        any.missing = FALSE,
        len = length(i)
      )
      solver_set_integrality(solver, as.integer(i) - 1L, as.integer(v))
    }
  }
  vbounds <- function(i, lower, upper) {
    if (missing(i) && missing(lower) && missing(upper)) {
      vbnd <- solver_get_variable_bounds(solver)
      matrix(vbnd, ncol = 2L, dimnames = list(NULL, c("lower", "upper")))
    } else {
      assert_integerish(i, lower = 1, any.missing = FALSE)
      assert_numeric(lower, any.missing = FALSE, len = length(i))
      assert_numeric(upper, any.missing = FALSE, len = length(i))
      solver_set_variable_bounds(solver, as.integer(i) - 1L, lower, upper)
    }
  }
  get_basis <- function() {
    solver_get_basis(solver)
  }
  set_basis <- function(col_status, row_status) {
    solver_set_basis(solver, as.integer(col_status), as.integer(row_status))
  }
  clear_basis <- function() {
    solver_clear_basis(solver)
  }
  set_solution <- function(col_value, row_value, col_dual, row_dual,
                           value_valid = TRUE, dual_valid = TRUE) {
    solver_set_solution_obj(solver, value_valid, dual_valid,
                            col_value, col_dual, row_value, row_dual)
  }
  get_dual_ray <- function() {
    solver_get_dual_ray(solver)
  }
  get_primal_ray <- function() {
    solver_get_primal_ray(solver)
  }
  maximum <- function(maximize) {
    if (missing(maximize)) {
      as.logical(solver_get_sense(solver))
    } else {
      checkmate::assert_logical(maximize, len = 1, any.missing = FALSE)
      solver_set_sense(solve, maximize)
    }
  }
  structure(environment(), class = "highs_solver_wrapper")
}


#' @noRd
#' @export
print.highs_solver_wrapper <- function(x, ...) {
  writeLines("High Solver Wrapper")
  writeLines("  - model")
  writeLines("  - solver")
  writeLines("  - control")
  exclude <- c("model", "solver", "control")
  for (name in setdiff(ls(x), exclude)) {
    xargs <- trimws(gsub(
      "function",
      "",
      capture.output(str(get(name, envir = x)))
    ))
    msg <- sprintf("  - %s%s", name, xargs)
    writeLines(msg)
  }
}
