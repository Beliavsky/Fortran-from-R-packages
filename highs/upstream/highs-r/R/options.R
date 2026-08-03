#' Available Solver Options
#'
#' Reference for the available solver options.
#'
#' @returns A \code{data.frame} containing the available solver options.
#' @examples
#' highs_available_solver_options()
#' @export
highs_available_solver_options <- function() {
  HIGHS_OPTIONS
}


solver_get_option_typed <- function(solver, key, type) {
  checkmate::assert_choice(type, c("bool", "integer", "double", "string"))
  if (type == "bool") {
    solver_get_bool_option(solver, key)
  } else if (type == "integer") {
    solver_get_int_option(solver, key)
  } else if (type == "double") {
    solver_get_dbl_option(solver, key)
  } else if (type == "string") {
    solver_get_str_option(solver, key)
  }
}


#' Get a HiGHS Solver Option
#'
#' Retrieves the value of a specific option from a HiGHS solver instance.
#'
#' @param solver A HiGHS solver object of class \code{"highs_solver"}.
#' @param key A character string specifying the option name to retrieve.
#' @param type Type of the option. Can be one of "auto", "bool", "integer", "double", or "string".
#'             When set to "auto" (default), the function will attempt to determine the type from the available options list.
#'             Specify a type directly if the option is valid but not listed in the available options.
#'
#' @return The value of the specified option with the appropriate type.
#' @examples
#' solver <- example_solver()
#' hi_solver_get_option(solver, "output_flag")
#' hi_solver_get_option(solver, "solver", type = "string")
#' @export
hi_solver_get_option <- function(solver, key, type = c("auto", "bool", "integer", "double", "string")) {
  checkmate::assert_class(solver, classes = "highs_solver")
  checkmate::assert_character(key, len = 1, any.missing = FALSE)
  type <- match.arg(type)
  if (type == "auto") {
    meta_data <- highs_available_solver_options()
    checkmate::assert_choice(key, meta_data[["option"]])
    type <- meta_data[["type"]][which(meta_data[["option"]] == key)]
  }
  solver_get_option_typed(solver, key, type)
}


#' Get multiple HiGHS Solver Options
#'
#' Retrieves the values of multiple options from a HiGHS solver instance.
#'
#' @param solver A HiGHS solver object of class \code{"highs_solver"}.
#' @param keys A character vector of option names to retrieve.
#'
#' @return A named list of option values with the appropriate types.
#'
#' @examples
#' solver <- example_solver()
#' hi_solver_get_options(solver, c("output_flag", "solver"))
#' @export
hi_solver_get_options <- function(solver, keys = NULL) {
  checkmate::assert_class(solver, classes = "highs_solver")
  checkmate::assert_character(keys, null.ok = TRUE)
  meta_data <- highs_available_solver_options()
  getters <- list(
    bool = solver_get_bool_option, integer = solver_get_int_option,
    double = solver_get_dbl_option, string = solver_get_str_option
  )
  if (length(keys) == 0L) {
    idx <- seq_len(nrow(meta_data))
  } else {
    idx <- match(meta_data[["option"]], keys)
    idx <- idx[!is.na(idx)]
    if (length(idx) != length(keys)) {
      emsg <- paste("The following options were not found in the available solver options:",
        paste(setdiff(keys, meta_data[["option"]]), collapse = ", "),
        "Please check the option names.",
        sep = "\n"
      )
      stop(emsg)
    }
  }

  opts <- setNames(vector("list", length(idx)), meta_data[["option"]][idx])
  for (i in idx) {
    opts[[i]] <- solver_get_option_typed(solver, meta_data[["option"]][i], meta_data[["type"]][i])
  }
  return(opts)
}


as_bool <- function(x) {
  if (is.character(x)) {
    if (tolower(x) %in% c("1", "on", "t")) {
      return(TRUE)
    }
    if (tolower(x) %in% c("0", "off", "f")) {
      return(FALSE)
    }
  }
  as.logical(x)
}


force_type <- function(obj, type) {
  if (is.na(type)) {
    return(obj)
  }
  force <- list(
    bool = as_bool,
    integer = as.integer,
    double = as.double,
    string = as.character
  )[[type]]
  force(obj)
}


test_type <- function(obj, type) {
  if (type == "bool") {
    checkmate::test_logical(obj, len = 1L, any.missing = FALSE)
  } else if (type == "integer") {
    checkmate::test_integer(obj, len = 1L, any.missing = FALSE)
  } else if (type == "double") {
    checkmate::test_double(obj, len = 1L, any.missing = FALSE)
  } else if (type == "string") {
    checkmate::test_string(obj)
  } else {
    stop(sprintf("type '%s' is not allowed", type))
  }
}



#' Set a HiGHS Solver Option
#'
#' Sets the value of a specific option for a HiGHS solver instance.
#'
#' @param solver A HiGHS solver object of class \code{"highs_solver"}.
#' @param key A character string specifying the option name to set.
#' @param value The value to set for the specified option. Will be coerced to the appropriate type.
#' @param type Type of the option. Can be one of "auto", "bool", "integer", "double", or "string".
#'             When set to "auto" (default), the function will attempt to determine the type from the available options list.
#'             Specify a type directly if the option is valid but not listed in the available options.
#'
#' @return Invisibly returns NULL.
#'
#' @examples
#' solver <- example_solver()
#' hi_solver_set_option(solver, "output_flag", FALSE)
#' hi_solver_set_option(solver, "solver", "simplex", type = "string")
#'
#' @export
hi_solver_set_option <- function(solver, key, value, type = c("auto", "bool", "integer", "double", "string")) {
  checkmate::assert_class(solver, classes = "highs_solver")
  checkmate::assert_character(key, len = 1, any.missing = FALSE)
  checkmate::assert_scalar(value)
  type <- match.arg(type)
  if (type == "auto") {
    meta_data <- highs_available_solver_options()
    checkmate::assert_choice(key, meta_data[["option"]])
    meta_data <- highs_available_solver_options()
    type <- meta_data[["type"]][which(meta_data[["option"]] == key)]
  }
  value <- force_type(value, type)
  solver_set_option(solver, key, value)
}


#' Set Multiple HiGHS Solver Options
#'
#' Sets multiple options for a HiGHS solver instance at once.
#'
#' @param solver A HiGHS solver object of class \code{"highs_solver"}.
#' @param control A named list of options to set. Names should be valid option names and values will be coerced to the appropriate types.
#'
#' @return Invisibly returns NULL.
#'
#' @examples
#' solver <- example_solver()
#' hi_solver_set_options(solver, list(output_flag = FALSE, solver = "simplex"))
#'
#' control <- list(
#'   presolve = "on",
#'   solver = "simplex",
#'   parallel = "on",
#'   ranging = "off",
#'   time_limit = 100.0,
#'   primal_feasibility_tolerance = 1e-7,
#'   dual_feasibility_tolerance = 1e-7,
#'   random_seed = 1234,
#'   threads = 4,
#'   output_flag = TRUE,
#'   log_to_console = TRUE,
#'   run_crossover = "on",
#'   allow_unbounded_or_infeasible = FALSE,
#'   mip_detect_symmetry = TRUE,
#'   mip_max_nodes = 10000,
#'   mip_max_leaves = 5000,
#'   mip_feasibility_tolerance = 1e-6
#' )
#' hi_solver_set_options(solver, control)
#' @export
hi_solver_set_options <- function(solver, control = list()) {
  if (length(control) == 0) {
    return(invisible(NULL))
  }
  checkmate::assert_class(solver, classes = "highs_solver")
  checkmate::assert_character(names(control), min.len = 1L, any.missing = FALSE)
  meta_data <- highs_available_solver_options()
  types <- meta_data[["type"]][match(names(control), meta_data[["option"]])]
  for (i in seq_along(control)) {
    key <- names(control)[i]
    value <- control[[i]]
    type <- types[i]
    solver_set_option(solver, key, force_type(value, type))
  }
  return(invisible(NULL))
}
