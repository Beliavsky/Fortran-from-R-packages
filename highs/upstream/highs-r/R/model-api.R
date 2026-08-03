#' Create new Highs Model
#'
#' @description
#' Create a new highs model object.
#' @return an object of class \code{"highs_model"}.
#' @examples
#' model <- hi_new_model() 
#' @export
hi_new_model <- function() {
    model <- new_model()
    class(model) <- c("highs_model", class(model))
    model
}


#' Sets the number of columns in the model
#'
#' This function sets the number of columns in the given model.
#'
#' @param model an object of class \code{"highs_model"}.
#' @param ncol an integer giving the number of columns (variables) to set in the model
#' @return \code{NULL}
#'
#' @examples
#' model <- hi_new_model()
#' hi_model_set_ncol(model, 10L) # Sets the model to have 10 columns
#' @export
hi_model_set_ncol <- function(model, ncol) {
    checkmate::assert_class(model, classes = "highs_model")
    checkmate::assert_integer(ncol, lower = 0, any.missing = FALSE, len = 1)
    model_set_ncol(model, ncol)
}


#' Set the number of rows in the model
#'
#' This function sets the number of rows in the given model.
#'
#' @param model an object of class \code{"highs_model"}.
#' @param nrow an integer giving the number of rows (variables) to set in the model
#' @return \code{NULL}
#'
#' @examples
#' model <- hi_new_model()
#' hi_model_set_nrow(model, 5L) # Sets the model to have 5 rows
#' @export
hi_model_set_nrow <- function(model, nrow) {
    checkmate::assert_class(model, classes = "highs_model")
    checkmate::assert_integer(nrow, lower = 0, any.missing = FALSE, len = 1)
    model_set_nrow(model, nrow)
}


#' Set the sense of the optimization model
#'
#' This function sets the sense of the optimization model to either maximization or minimization.
#'
#' @param model an object of class \code{"highs_model"}.
#' @param maximum a boolean value indicating whether the model should be set to maximization (`TRUE`) or minimization (`FALSE`).
#' @return \code{NULL}
#'
#' @examples
#' model <- hi_new_model()
#' hi_model_set_sense(model, TRUE) # Set the model to maximization
#' hi_model_set_sense(model, FALSE) # Set the model to minimization
#'
#' @export
hi_model_set_sense <- function(model, maximum) {
    checkmate::assert_class(model, classes = "highs_model")
    checkmate::assert_logical(maximum, len = 1, any.missing = FALSE)
    model_set_sense(model, maximum)
}



#' Set Offset for Highs Model
#'
#' This function sets the offset for a given Highs model.
#'
#' @param model an object of class \code{"highs_model"}.
#' @param offset a numeric value of length 1. The offset value to be set for the model.
#'
#' @return \code{NULL}
#'
#' @examples
#' model <- hi_new_model()
#' hi_model_set_offset(model, 10)
#'
#' @export
hi_model_set_offset <- function(model, offset) {
    checkmate::assert_class(model, classes = "highs_model")
    checkmate::assert_numeric(offset, len = 1, any.missing = FALSE)
    model_set_offset(model, offset)
}



#' Set Objective for Highs Model
#'
#' This function sets the objective for a given Highs model.
#'
#' @param model an object of class \code{"highs_model"}.
#' @param objective a numeric vector giving the objective values to be set for the model.
#'
#' @return \code{NULL}
#'
#' @examples
#' model <- hi_new_model()
#' objective <- c(1, 2, 3)
#' hi_model_set_objective(model, objective)
#'
#' @export
hi_model_set_objective <- function(model, objective) {
    checkmate::assert_class(model, classes = "highs_model")
    checkmate::assert_numeric(objective, any.missing = FALSE)
    model_set_objective(model, objective)
}


#' Set Lower Bounds for Highs Model
#'
#' This function sets the lower bounds for a given Highs model.
#'
#' @param model an object of class \code{"highs_model"}.
#' @param lower a numeric vector giving the lower bounds.
#'
#' @return \code{NULL}
#'
#' @examples
#' model <- hi_new_model()
#' lower_bounds <- c(0, 1, 2)
#' hi_model_set_lower(model, lower_bounds)
#'
#' @export
hi_model_set_lower <- function(model, lower) {
    checkmate::assert_class(model, classes = "highs_model")
    checkmate::assert_numeric(lower, any.missing = FALSE)
    model_set_lower(model, lower)
}


#' Set Upper Bounds for a Highs Model
#'
#' This function sets the upper bounds for a given Highs model.
#'
#' @param model an object of class \code{"highs_model"}.
#' @param upper a numeric vector giving the upper bounds.
#'
#' @return \code{NULL}
#'
#' @examples
#' model <- hi_new_model()
#' upper_bounds <- c(10, 20, 30)
#' hi_model_set_upper(model, upper_bounds)
#'
#' @export
hi_model_set_upper <- function(model, upper) {
    checkmate::assert_class(model, classes = "highs_model")
    checkmate::assert_numeric(upper, any.missing = FALSE)
    model_set_upper(model, upper)
}


#' Set Constraint Matrix for Highs Model
#'
#' This function sets the constraint matrix for a given Highs model.
#'
#' @param model an object of class \code{"highs_model"}.
#' @param matrix a matrix giving the Hessian matrix.
#'  Allowed matrix classes are \code{"matrix"}, \code{"dgCMatrix"}, \code{"matrix.csc"},
#'  and \code{"simple_triplet_matrix"}.
#'
#' @return \code{NULL}
#'
#' @examples
#' model <- hi_new_model()
#' matrix <- matrix(c(1, 0, 0, 1), nrow = 2)
#' hi_model_set_constraint_matrix(model, matrix)
#'
#' @export
hi_model_set_constraint_matrix <- function(model, matrix) {
    checkmate::assert_class(model, classes = "highs_model")
    checkmate::assert_numeric(matrix, any.missing = FALSE)
    model_set_constraint_matrix(model, matrix)
}


#' Set Left Hand Side for a Highs Model
#'
#' This function sets the left hand side for a given Highs model.
#'
#' @param model an object of class \code{"highs_model"}.
#' @param lhs a numeric vector giving the left hand side values.
#'
#' @return \code{NULL}
#'
#' @examples
#' model <- hi_new_model()
#' model <- hi_model_set_lhs(model, c(0, 1, 2))
#'
#' @export
hi_model_set_lhs <- function(model, lhs) {
    checkmate::assert_class(model, classes = "highs_model")
    checkmate::assert_numeric(lhs, any.missing = FALSE)
    model_set_lhs(model, lhs)
}


#' Set Right Hand Side for a Highs Model
#'
#' This function sets the left hand side for a given Highs model.
#'
#' @param model an object of class \code{"highs_model"}.
#' @param rhs a numeric vector giving the left hand side values.
#'
#' @return \code{NULL}
#'
#' @examples
#' model <- hi_new_model()
#' model <- hi_model_set_rhs(model, c(0, 1, 2))
#'
#' @export
hi_model_set_rhs <- function(model, rhs) {
    checkmate::assert_class(model, classes = "highs_model")
    checkmate::assert_numeric(rhs, any.missing = FALSE)
    model_set_rhs(model, rhs)
}


#' Set Hessian Matrix for Highs Model
#'
#' This function sets the Hessian matrix for a given Highs model.
#'
#' @param model an object of class \code{"highs_model"}.
#' @param matrix a matrix giving the Hessian matrix.
#'  Allowed matrix classes are \code{"matrix"}, \code{"dgCMatrix"}, \code{"matrix.csc"},
#'  and \code{"simple_triplet_matrix"}.
#'
#' @return \code{NULL}
#'
#' @examples
#' model <- hi_new_model()
#' hessian_matrix <- matrix(c(1, 0, 0, 1), nrow = 2)
#' hi_model_set_hessian(model, hessian_matrix)
#'
#' @export
hi_model_set_hessian <- function(model, matrix) {
    checkmate::assert_class(model, classes = "highs_model")
    checkmate::assert_numeric(matrix, any.missing = FALSE)
    model_set_hessian(model, matrix)
}


#' Set Variable Types in a Highs Model
#'
#' This function sets the variable types in a given Highs model.
#'
#' @param model an object of class \code{"highs_model"}.
#' @param types an integer vector specifying the types of the variables.
#'
#' @return The function does not return a value. It modifies the `model` object in place.
#'
#' @examples
#' model <- hi_new_model()
#' types <- c(1, 2, 1, 0)
#' hi_model_set_vartype(model, types)
#'
#' @export
hi_model_set_vartype <- function(model, types) {
    checkmate::assert_class(model, classes = "highs_model")
    checkmate::assert_integerish(types, any.missing = FALSE)
    model_set_vartype(model, as.integer(types))
}


#' Get Number of Variables in a Highs Model
#'
#' This function retrieves the number of variables in a given Highs model.
#'
#' @param model A `highs_model` object. The model from which to get the number of variables.
#'
#' @return An integer representing the number of variables in the model.
#'
#' @examples
#' model <- hi_new_model()
#' hi_model_get_nvars(model)
#'
#' @export
hi_model_get_nvars <- function(model) {
    checkmate::assert_class(model, classes = "highs_model")
    model_get_nvars(model)
}


#' Get Number of Constraints in a Model
#'
#' This function retrieves the number of constraints in a given `highs_model` object.
#'
#' @param model A `highs_model` object. The model from which to get the number of variables.
#' 
#' @return An integer representing the number of constraints in the model.
#' @examples
#' model <- hi_new_model()
#' hi_model_get_ncons(model)
#' 
#' @export
hi_model_get_ncons <- function(model) {
    checkmate::assert_class(model, classes = "highs_model")
    model_get_ncons(model)
}
