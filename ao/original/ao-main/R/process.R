#' Process Object
#'
#' @description
#' This object specifies an AO process.
#'
#' @param target \[`character()` | `NULL`\]\cr
#' The name(s) of the argument(s) over which \code{f} is optimized.
#'
#' This can only be \code{numeric} arguments.
#'
#' If `NULL` (default), the first argument of `f` is optimized.
#'
#' @param npar \[`integer()`\]\cr
#' The length(s) of the target argument(s).
#'
#' @param partition \[`character(1)` | `list()`\]\cr
#' Defines the parameter partition. It can be
#'
#' * `"sequential"` for treating each parameter separately,
#' * `"random"` for a random partition in each iteration,
#' * `"none"` for no partition (which is equivalent to joint optimization),
#' * or a `list` of vectors of parameter indices, specifying a custom
#'   partition for the AO process.
#'
#' @param new_block_probability \[`numeric(1)`\]\cr
#' Only relevant if `partition = "random"`.
#'
#' The probability of creating a new parameter block when creating a random
#' partition.
#'
#' Values close to 0 result in larger parameter blocks, values close to 1
#' result in smaller parameter blocks.
#'
#' @param minimum_block_number \[`integer(1)`\]\cr
#' Only relevant if `partition = "random"`.
#'
#' The minimum number of blocks in random partitions.
#'
#' @param verbose \[`logical(1)`\]\cr
#' Print tracing details during the AO process?
#'
#' @param minimize \[`logical(1)`\]\cr
#' Minimize during the AO process?
#'
#' If \code{FALSE}, maximization is performed.
#'
#' @param iteration_limit \[`integer(1)` | `Inf`\]\cr
#' The maximum number of iterations through the parameter partition before
#' the AO process is terminated.
#'
#' Can also be `Inf` for no iteration limit.
#'
#' @param seconds_limit \[`numeric(1)`\]\cr
#' The time limit in seconds before the AO process is terminated.
#'
#' Can also be `Inf` for no time limit.
#'
#' Note that this stopping criterion is only checked *after* a sub-problem is
#' solved and not *within* solving a sub-problem, so the actual process
#' time can exceed this limit.
#'
#' @param tolerance_value \[`numeric(1)`\]\cr
#' A non-negative tolerance value. The AO process terminates
#' if the absolute difference between the current function value and the value
#' from \code{tolerance_history} iterations earlier is smaller than
#' \code{tolerance_value}.
#'
#' Can be `0` for no value threshold.
#'
#' @param tolerance_parameter \[`numeric(1)`\]\cr
#' A non-negative tolerance value. The AO process terminates if
#' the distance between the current estimate and the estimate from
#' \code{tolerance_history} iterations earlier is smaller than
#' \code{tolerance_parameter}.
#'
#' Can be `0` for no parameter threshold.
#'
#' By default, the distance is measured using the Euclidean norm, but another
#' norm can be specified via the \code{tolerance_parameter_norm} field.
#'
#' @param tolerance_parameter_norm \[`function`\]\cr
#' The norm that measures the distance between two estimates. If the distance
#' is smaller than
#' \code{tolerance_parameter}, the AO process is terminated.
#'
#' It must be of the form \code{function(x, y)} for two vector inputs
#' \code{x} and \code{y}, and return a single \code{numeric} value.
#' By default, the Euclidean norm \code{function(x, y) sqrt(sum((x - y)^2))}
#' is used.
#'
#' @param tolerance_history \[`integer(1)`\]\cr
#' The number of iterations to look back to determine whether
#' \code{tolerance_value} or \code{tolerance_parameter} has been reached.
#'
#' @param which_iteration \[`integer()`\]\cr
#' Selects the iteration(s).
#'
#' Can also be \code{NULL} to select all iterations.
#'
#' @param which_block \[`character(1)` | `integer()`\]\cr
#' Selects the parameter block in the partition and can be one of
#'
#' - \code{"first"} for the first parameter block,
#' - \code{"last"} for the last parameter block,
#' - an `integer` vector of parameter indices,
#' - or \code{NULL} for all parameter blocks.
#'
#' @param which_column \[`character()`\]\cr
#' Selects the columns in the `details` part of the output and can be one or
#' more of
#'
#' - `"iteration"`,
#' - `"value"`,
#' - `"parameter"`,
#' - `"block"`,
#' - and `"seconds"`.
#'
#' @param keep_iteration_column \[`logical(1)`\]\cr
#' Keep the column containing the information about the iteration in the output?
#'
#' @param keep_block_columns \[`logical(1)`\]\cr
#' Keep the column containing the information about the active parameter block
#' in the output?
#'
#' @param parameter_type \[`character(1)`\]\cr
#' Selects the parameter type and can be one of
#'
#' - \code{"full"} (default) to get the full parameter vector,
#' - \code{"block"} to get the parameter values for the current block,
#'   i.e., the parameters with the indices `self$block`
#' - \code{"fixed"} to get the parameter values which are currently fixed,
#'   i.e., all except for those with the indices `self$block`
#'
#' @param add_details \[`logical(1)`\]\cr
#' Add details about the AO process to the output?

Process <- R6::R6Class("Process",
   cloneable = FALSE,
   public = list(

   #' @description
   #' Creates a new object of this [R6][R6::R6Class] class.

   initialize = function(
     target = NULL,
     npar = integer(),
     partition = "sequential",
     new_block_probability = 0.3,
     minimum_block_number = 1,
     verbose = FALSE,
     minimize = TRUE,
     iteration_limit = Inf,
     seconds_limit = Inf,
     tolerance_value = 1e-6,
     tolerance_parameter = 1e-6,
     tolerance_parameter_norm = function(x, y) sqrt(sum((x - y)^2)),
     tolerance_history = 1,
     add_details = TRUE
    ) {
     oeli::input_check_response(
       check = checkmate::check_integerish(npar, lower = 0),
       var_name = "npar"
     )
     oeli::input_check_response(
       check = checkmate::check_character(
         target, len = length(npar), unique = TRUE, null.ok = TRUE
        ),
       var_name = "target"
     )
     private$.npar <- npar
     private$.target <- target
     self$partition <- partition
     self$new_block_probability <- new_block_probability
     self$minimum_block_number <- minimum_block_number
     self$verbose <- verbose
     self$minimize <- minimize
     self$iteration_limit <- iteration_limit
     self$seconds_limit <- seconds_limit
     self$tolerance_value <- tolerance_value
     self$tolerance_parameter <- tolerance_parameter
     self$tolerance_parameter_norm <- tolerance_parameter_norm
     self$tolerance_history <- tolerance_history
     self$add_details <- add_details
     invisible(self)
   },

   #' @description
   #' Prints a status message.
   #'
   #' @param message \[`character(1)`\]\cr
   #' A status message.
   #'
   #' @param message_type \[`integer(1)`\]\cr
   #' The message type, one of the following:
   #'
   #' * `1` for `cli::cli_h1()`
   #' * `2` for `cli::cli_h2()`
   #' * `3` for `cli::cli_h3()`
   #' * `4` for `cli::cli_alert_success()`
   #' * `5` for `cli::cli_alert_info()`
   #' * `6` for `cli::cli_alert_warning()`
   #' * `7` for `cli::cli_alert_danger()`
   #' * `8` for `cli::cat_line()`

   print_status = function(message, message_type = 8, verbose = self$verbose) {
     checkmate::assert_string(message, min.chars = 1)
     checkmate::assert_int(message_type, lower = 1, upper = 8)
     checkmate::assert_flag(verbose)
     if (isTRUE(verbose)) {
       switch(message_type,
              cli::cli_h1(message),
              cli::cli_h2(message),
              cli::cli_h3(message),
              cli::cli_alert_success(message),
              cli::cli_alert_info(message),
              cli::cli_alert_warning(message),
              cli::cli_alert_danger(message),
              cli::cat_line(message)
       )
     }
     invisible(self)
   },

   #' @description
   #' Initializes the `details` part of the output.
   #'
   #' @param initial_parameter \[`numeric()`\]\cr
   #' The starting parameter values for the AO process.
   #'
   #' @param initial_value \[`numeric(1)`\]\cr
   #' The function value at the initial parameters.

   initialize_details = function(initial_parameter, initial_value) {
     private$.details <- structure(
       data.frame(
         t(c(0L, initial_value, initial_parameter, rep(0, sum(self$npar)), 0))
       ),
       names = c(
         "iteration", "value", paste0("p", seq_len(sum(self$npar))),
         paste0("b", seq_len(sum(self$npar))), "seconds"
       )
     )
     invisible(self)
   },

   #' @description
   #' Updates the `details` part of the output.
   #'
   #' @param value \[`numeric(1)`\]\cr
   #' The updated function value.
   #'
   #' @param parameter_block \[`numeric()`\]\cr
   #' The updated parameter values for the active parameter block.
   #'
   #' @param seconds \[`numeric(1)`\]\cr
   #' The time in seconds for solving the sub-problem.
   #'
   #' @param error \[`logical(1)`\]\cr
   #' Did solving the sub-problem result in an error?
   #'
   #' @param error_message \[`character(1)`\]\cr
   #' An error message if `error = TRUE`.
   #'
   #' @param block \[`integer()`\]\cr
   #' The currently active parameter block, represented as parameter indices.

   update_details = function(
     value, parameter_block, seconds, error, error_message, block = self$block
   ) {

     ### check inputs
     check_block <- checkmate::check_integerish(
       block,
       unique = TRUE, lower = 1, upper = sum(self$npar)
     )
     check_value <- checkmate::check_number(
       value,
       na.ok = FALSE, null.ok = FALSE, finite = TRUE
     )
     check_parameter_block <- oeli::check_numeric_vector(
       parameter_block,
       any.missing = FALSE, len = length(block)
     )
     check_seconds <- checkmate::check_number(
       seconds,
       na.ok = FALSE, lower = 0, finite = TRUE, null.ok = FALSE
     )
     check_error <- checkmate::check_flag(error)
     check_results <- all(sapply(
       c(
         check_block, check_value, check_parameter_block, check_seconds,
         check_error, isFALSE(error)
       ),
       isTRUE
     ))

     ### error occurred?
     if (isTRUE(error)) {
       private$.error <- TRUE
       private$.error_message <- paste(
         "solving for block [", paste(block, collapse = ","), "] failed:",
         if (checkmate::test_string(error_message, min.chars = 1)) {
           error_message
         } else {
           "unknown reason"
         }
       )
     }

     ### update details
     rows <- nrow(private$.details)
     if (isTRUE(check_results)) {
       parameter <- self$get_parameter_latest("full")
       parameter[block] <- parameter_block
       private$.details[rows + 1, "value"] <- value
       parameter_columns <- which(startsWith(colnames(private$.details), "p"))
       private$.details[rows + 1, parameter_columns] <- parameter
     } else {
       private$.details[rows + 1, ] <- private$.details[rows, ]
     }
     private$.details[rows + 1, "iteration"] <- self$iteration
     private$.details[rows + 1, "seconds"] <- seconds
     block_columns <- which(startsWith(colnames(private$.details), "b"))
     private$.details[rows + 1, block_columns[block]] <- 1
     private$.details[rows + 1, block_columns[-block]] <- 0

     ### return information
     self$print_status("function value:", 8)
     if (self$verbose) cli::cat_print(value)
     self$print_status("parameter:", 8)
     if (self$verbose) cli::cat_print(self$get_parameter_latest("full"))
     invisible(self)
   },

   #' @description
   #' Get a parameter partition.

   get_partition = function() {
     if (checkmate::test_string(self$partition)) {
       switch(
         self$partition,
         "sequential" = as.list(seq_len(sum(self$npar))),
         "random"     = generate_random_partition(
            x = seq_len(sum(self$npar)),
            p = self$new_block_probability,
            min = self$minimum_block_number
          ),
          "none"       = list(seq_len(sum(self$npar)))
       )
     } else {
       self$partition
     }
   },

   #' @description
   #' Get the `details` part of the output.

   get_details = function(
      which_iteration = NULL,
      which_block = NULL,
      which_column = c("iteration", "value", "parameter", "block", "seconds")
    ) {

     ### input checks
     oeli::input_check_response(
       check = checkmate::check_integerish(
         which_iteration,
         lower = 0, null.ok = TRUE, min.len = 1, any.missing = FALSE
       ),
       var_name = "which_iteration"
     )
     oeli::input_check_response(
       check = list(
         checkmate::test_choice(
           which_block, c("first", "last"), null.ok = TRUE
         ),
         checkmate::test_integerish(
           which_block,
           lower = 1, upper = sum(self$npar), unique = TRUE,
           min.len = 1, max.len = sum(self$npar), any.missing = FALSE
         )
       ),
       var_name = "which_block"
     )
     oeli::input_check_response(
       check = checkmate::check_subset(
         which_column, c("iteration", "value", "parameter", "block", "seconds"),
         empty.ok = TRUE
       ),
       var_name = "which_column"
     )

     ### filter details
     details <- private$.details

     # define block columns early (needed for filtering)
     block_columns <- which(startsWith(colnames(details), "b"))

     if (!is.null(which_iteration)) {
       ids <- details$iteration %in% which_iteration
       details <- details[ids, , drop = FALSE]
     }

     if (!is.null(which_block)) {
       if (is.character(which_block)) {
         rows <- oeli::vector_occurrence(details$iteration, which_block)
       } else {
         target <- sort(as.integer(which_block))
         block_details <- details[, block_columns, drop = FALSE]
         rows <- which(apply(block_details, 1, function(x) {
           active <- which(as.numeric(x) == 1)
           identical(active, target)
         }))
       }
     } else {
       rows <- seq_len(nrow(details))
     }

     ### select columns
     iteration_column <- which(colnames(details) == "iteration")
     value_column <- which(colnames(details) == "value")
     parameter_columns <- which(startsWith(colnames(details), "p"))
     block_columns <- which(startsWith(colnames(details), "b"))
     seconds_column <- which(colnames(details) == "seconds")
     columns <- c(
       if ("iteration" %in% which_column) iteration_column,
       if ("value" %in% which_column) value_column,
       if ("parameter" %in% which_column) parameter_columns,
       if ("block" %in% which_column) block_columns,
       if ("seconds" %in% which_column) seconds_column
     )

     ### return details
     details[rows, columns]
   },

   #' @description
   #' Get the function value in different steps of the AO process.

   get_value = function(which_iteration = NULL,
                        which_block = NULL,
                        keep_iteration_column = FALSE,
                        keep_block_columns = FALSE) {
     ### input checks
     oeli::input_check_response(
       check = checkmate::test_flag(keep_iteration_column),
       var_name = "keep_iteration_column"
     )
     oeli::input_check_response(
       check = checkmate::test_flag(keep_block_columns),
       var_name = "keep_block_columns"
     )

     ### return values
     self$get_details(
       which_iteration = which_iteration, which_block = which_block,
       which_column = c(
         "value",
         if (keep_iteration_column) "iteration",
         if (keep_block_columns) "block"
       )
     )
   },

   #' @description
   #' Get the function value in the latest step of the AO process.

   get_value_latest = function() {
     private$.details[nrow(private$.details), "value"]
   },

   #' @description
   #' Get the optimum function value in the AO process.

   get_value_best = function() {
     values <- private$.details[, "value"]
     id_best <- if (self$minimize) {
       which.min(values)
     } else {
       which.max(values)
     }
     private$.details[id_best, "value"]
   },

   #' @description
   #' Get the parameter values in different steps of the AO process.

   get_parameter = function(which_iteration = self$iteration,
                            which_block = NULL,
                            keep_iteration_column = FALSE,
                            keep_block_columns = FALSE) {
     ### input checks
     oeli::input_check_response(
       check = checkmate::test_flag(keep_iteration_column),
       var_name = "keep_iteration_column"
     )
     oeli::input_check_response(
       check = checkmate::test_flag(keep_block_columns),
       var_name = "keep_block_columns"
     )

     ### get parameters
     self$get_details(
       which_iteration = which_iteration, which_block = which_block,
       which_column = c(
         "parameter",
         if (keep_iteration_column) "iteration",
         if (keep_block_columns) "block"
       )
     )
   },

   #' @description
   #' Get the parameter value in the latest step of the AO process.

   get_parameter_latest = function(parameter_type = "full") {

     ### input checks
     oeli::input_check_response(
       check = checkmate::check_choice(
         parameter_type,
         c("full", "block", "fixed")
       ),
       var_name = "parameter_type"
     )

     ### get latest parameter
     details <- private$.details
     parameter_columns <- which(startsWith(colnames(details), "p"))
     latest_parameter <- as.numeric(details[nrow(details), parameter_columns])

     ### filter for 'parameter_type'
     if (parameter_type == "full") {
       return(latest_parameter)
     } else {
       if (parameter_type == "block") {
         return(latest_parameter[self$block])
       } else {
         return(latest_parameter[-self$block])
       }
     }
   },

   #' @description
   #' Get the optimum parameter value in the AO process.

   get_parameter_best = function(parameter_type = "full") {

     ### input checks
     oeli::input_check_response(
       check = checkmate::check_choice(
         parameter_type,
         c("full", "block", "fixed")
       ),
       var_name = "parameter_type"
     )

     ### find best parameter
     details <- private$.details
     id_best <- if (self$minimize) {
       which.min(details[, "value"])
     } else {
       which.max(details[, "value"])
     }
     parameter_columns <- which(startsWith(colnames(details), "p"))
     best_parameter <- as.numeric(details[id_best, parameter_columns])

     ### filter for 'parameter_type'
     if (parameter_type == "full") {
       return(best_parameter)
     } else {
       if (parameter_type == "block") {
         return(best_parameter[self$block])
       } else {
         return(best_parameter[-self$block])
       }
     }
   },

   #' @description
   #' Get the optimization time in seconds in different steps of the AO process.

   get_seconds = function(which_iteration = NULL,
                          which_block = NULL,
                          keep_iteration_column = FALSE,
                          keep_block_columns = FALSE) {

     ### input checks
     oeli::input_check_response(
       check = checkmate::test_flag(keep_iteration_column),
       var_name = "keep_iteration_column"
     )
     oeli::input_check_response(
       check = checkmate::test_flag(keep_block_columns),
       var_name = "keep_block_columns"
     )

     ### get parameters
     self$get_details(
       which_iteration = which_iteration, which_block = which_block,
       which_column = c(
         "seconds",
         if (keep_iteration_column) "iteration",
         if (keep_block_columns) "block"
       )
     )
   },

   #' @description
   #' Get the total optimization time in seconds of the AO process.

   get_seconds_total = function() {
     sum(self$get_seconds(
       which_iteration = NULL, which_block = NULL,
       keep_iteration_column = FALSE, keep_block_columns = FALSE
     ), na.rm = TRUE)
   },

   #' @description
   #' Checks if the AO process can be terminated.

   check_stopping = function() {

     ### check stopping criteria
     stopping <- FALSE
     while (TRUE) {

       ### check error occurred
       if (isTRUE(private$.error)) {
         message <- private$.error_message
         stopping <- TRUE
         break
       }

       ### check iteration limit
       if (self$iteration >= self$iteration_limit) {
         message <- paste("iteration limit of", self$iteration_limit, "reached")
         stopping <- TRUE
         break
       }

       ### check time limit
       if (self$get_seconds_total() >= self$seconds_limit) {
         message <- paste(
           "time limit of", self$seconds_limit, "seconds reached"
         )
         stopping <- TRUE
         break
       }

       ### only check tolerance if at least 'tolerance_history' iterations have
       ### been performed
       if (self$iteration >= self$tolerance_history) {
         ### check value tolerance
         abs_value_change <- abs(
           self$get_value_latest() - self$get_value(
             which_iteration = self$iteration - self$tolerance_history,
             which_block = "first"
           )
         )
         if (abs_value_change < self$tolerance_value) {
           message <- paste(
             "change in function value between", self$tolerance_history,
             ifelse(self$tolerance_history > 1, "iterations", "iteration"),
             "is <", self$tolerance_value
           )
           stopping <- TRUE
           break
         }

         ### check parameter tolerance
         parameter_change <- self$tolerance_parameter_norm(
           self$get_parameter_latest(),
           self$get_parameter(
             which_iteration = self$iteration - self$tolerance_history,
             which_block = "first"
           )
         )
         if (parameter_change < self$tolerance_parameter) {
           message <- paste(
             "change in function parameters between", self$tolerance_history,
             ifelse(self$tolerance_history > 1, "iterations", "iteration"),
             "is <", self$tolerance_parameter
           )
           stopping <- TRUE
           break
         }
       }
       break
     }

     ### decide for stopping
     if (isTRUE(stopping)) {
       private$.stopping_reason <- message
       self$print_status("\n", 8)
       self$print_status(paste("process is terminated:", message), 4)
       self$print_status("\n", 8)
     }
     return(stopping)
   }

 ),

 active = list(

   #' @field npar \[`integer()`\]\cr
   #' The length(s) of the target argument(s).

   npar = function(value) {
     if (missing(value)) {
       private$.npar
     } else {
       cli::cli_abort(
         "Field {.var npar} is read-only",
         call = NULL
       )
     }
   },

   #' @field partition \[`character(1)` | `list()`\]\cr
   #' Defines the parameter partition. It can be
   #'
   #' * `"sequential"` for treating each parameter separately,
   #' * `"random"` for a random partition in each iteration,
   #' * `"none"` for no partition (which is equivalent to joint optimization),
   #' * or a `list` of vectors of parameter indices, specifying a custom
   #'   partition for the AO process.

   partition = function(value) {
     if (missing(value)) {
       private$.partition
     } else {
       if (checkmate::test_string(value)) {
         oeli::input_check_response(
           check = checkmate::check_choice(
             value,
             c("sequential", "random", "none")
           ),
           var_name = "partition"
         )
       } else if (checkmate::test_list(value)) {
         oeli::input_check_response(
           check = checkmate::check_list(value),
           var_name = "partition"
         )

         ### only parameter indices 1,...,'sum(self$npar)' allowed
         oeli::input_check_response(
           check = checkmate::check_integerish(
             unlist(value),
             lower = 1, upper = sum(self$npar), any.missing = FALSE
           ),
           var_name = "partition",
           prefix = "Elements in {.cls list} {.var partition} are bad:"
         )
       } else {
         oeli::input_check_response(
           check = "Bad specification",
           var_name = "partition"
         )
       }
       private$.partition <- value
     }
   },

   #' @field new_block_probability \[`numeric(1)`\]\cr
   #' Only relevant if `partition = "random"`.
   #'
   #' The probability of creating a new parameter block when creating a random
   #' partition.
   #'
   #' Values close to 0 result in larger parameter blocks, values close to 1
   #' result in smaller parameter blocks.

   new_block_probability = function(value) {
     if (missing(value)) {
       private$.new_block_probability
     } else {
       oeli::input_check_response(
         check = checkmate::check_number(value, lower = 0, upper = 1),
         var_name = "new_block_probability"
       )
       private$.new_block_probability <- value
     }
   },

   #' @field minimum_block_number \[`integer(1)`\]\cr
   #' Only relevant if `partition = "random"`.
   #'
   #' The minimum number of blocks in random partitions.

   minimum_block_number = function(value) {
     if (missing(value)) {
       private$.minimum_block_number
     } else {
       oeli::input_check_response(
         check = checkmate::check_int(
           value, lower = 1, upper = sum(self$npar)
         ),
         var_name = "minimum_block_number"
       )
       private$.minimum_block_number <- value
     }
   },

   #' @field verbose \[`logical(1)`\]\cr
   #' Print tracing details during the AO process?

   verbose = function(value) {
     if (missing(value)) {
       private$.verbose
     } else {
       oeli::input_check_response(
         check = checkmate::check_flag(value),
         var_name = "verbose"
       )
       private$.verbose <- value
     }
   },

   #' @field minimize \[`logical(1)`\]\cr
   #' Minimize during the AO process?
   #'
   #' If \code{FALSE}, maximization is performed.

   minimize = function(value) {
     if (missing(value)) {
       private$.minimize
     } else {
       oeli::input_check_response(
         check = checkmate::check_flag(value),
         var_name = "minimize"
       )
       private$.minimize <- value
     }
   },

   #' @field iteration_limit \[`integer(1)` | `Inf`\]\cr
   #' The maximum number of iterations through the parameter partition before
   #' the AO process is terminated.
   #'
   #' Can also be `Inf` for no iteration limit.

   iteration_limit = function(value) {
     if (missing(value)) {
       private$.iteration_limit
     } else {
       oeli::input_check_response(
         check = checkmate::check_number(value, lower = 1, finite = FALSE),
         var_name = "iteration_limit"
       )
       if (is.finite(value)) {
         value <- as.integer(value)
       }
       private$.iteration_limit <- value
     }
   },

   #' @field seconds_limit \[`numeric(1)`\]\cr
   #' The time limit in seconds before the AO process is terminated.
   #'
   #' Can also be `Inf` for no time limit.
   #'
   #' Note that this stopping criterion is only checked *after* a sub-problem is
   #' solved and not *within* solving a sub-problem, so the actual process
   #' time can exceed this limit.

   seconds_limit = function(value) {
     if (missing(value)) {
       private$.seconds_limit
     } else {
       oeli::input_check_response(
         check = checkmate::check_number(
           value,
           lower = 0, finite = FALSE, null.ok = FALSE, na.ok = FALSE
         ),
         var_name = "seconds_limit"
       )
       private$.seconds_limit <- value
     }
   },

   #' @field tolerance_value \[`numeric(1)`\]\cr
   #' A non-negative tolerance value. The AO process terminates
   #' if the absolute difference between the current function value and the
   #' value from \code{tolerance_history} iterations earlier is smaller than
   #' \code{tolerance_value}.
   #'
   #' Can be `0` for no value threshold.

   tolerance_value = function(value) {
     if (missing(value)) {
       private$.tolerance_value
     } else {
       oeli::input_check_response(
         check = checkmate::check_number(value, lower = 0, finite = TRUE),
         var_name = "tolerance_value"
       )
       private$.tolerance_value <- value
     }
   },

   #' @field tolerance_parameter \[`numeric(1)`\]\cr
   #' A non-negative tolerance value. The AO process terminates if
   #' the distance between the current estimate and the estimate from
   #' \code{tolerance_history} iterations earlier is smaller than
   #' \code{tolerance_parameter}.
   #'
   #' Can be `0` for no parameter threshold.
   #'
   #' By default, the distance is measured using the Euclidean norm, but another
   #' norm can be specified via the \code{tolerance_parameter_norm} field.

   tolerance_parameter = function(value) {
     if (missing(value)) {
       private$.tolerance_parameter
     } else {
       oeli::input_check_response(
         check = checkmate::check_number(value, lower = 0, finite = FALSE),
         var_name = "tolerance_parameter"
       )
       private$.tolerance_parameter <- value
     }
   },

   #' @field tolerance_parameter_norm \[`function`\]\cr
   #' The norm that measures the distance between two estimates. If the distance
   #' is smaller than
   #' \code{tolerance_parameter}, the AO process is terminated.
   #'
   #' It must be of the form \code{function(x, y)} for two vector inputs
   #' \code{x} and \code{y}, and return a single \code{numeric} value.
   #' By default, the Euclidean norm \code{function(x, y) sqrt(sum((x - y)^2))}
   #' is used.

   tolerance_parameter_norm = function(value) {
     if (missing(value)) {
       private$.tolerance_parameter_norm[[1]]
     } else {
       oeli::input_check_response(
         check = checkmate::check_function(
           value,
           args = c("x", "y"),
           nargs = 2
         ),
         var_name = "tolerance_parameter_norm"
       )
       private$.tolerance_parameter_norm <- list(value)
     }
   },

   #' @field tolerance_history \[`integer(1)`\]\cr
   #' The number of iterations to look back to determine whether
   #' \code{tolerance_value} or \code{tolerance_parameter} has been reached.

   tolerance_history = function(value) {
     if (missing(value)) {
       private$.tolerance_history
     } else {
       oeli::input_check_response(
         check = checkmate::check_count(value, positive = TRUE),
         var_name = "tolerance_history"
       )
       private$.tolerance_history <- value
     }
   },

   #' @field add_details \[`logical(1)`\]\cr
   #' Add details about the AO process to the output?

   add_details = function(value) {
     if (missing(value)) {
       private$.add_details
     } else {
       oeli::input_check_response(
         check = checkmate::check_flag(value),
         var_name = "add_details"
       )
       private$.add_details <- value
     }
   },

   #' @field iteration \[`integer(1)`\]\cr
   #' The current iteration number.

   iteration = function(value) {
     if (missing(value)) {
       private$.iteration
     } else {
       oeli::input_check_response(
         check = checkmate::check_int(value, lower = 0, na.ok = FALSE),
         var_name = "iteration"
       )
       private$.iteration <- as.integer(value)
       self$print_status(
         paste("iteration", private$.iteration, "of", self$iteration_limit), 2
       )
     }
   },

   #' @field block \[`integer()`\]\cr
   #' The currently active parameter block, represented as parameter indices.

   block = function(value) {
     if (missing(value)) {
       private$.block
     } else {
       oeli::input_check_response(
         check = checkmate::check_integerish(
           value,
           unique = TRUE, lower = 1, upper = sum(self$npar)
         ),
         var_name = "block"
       )
       private$.block <- value
       self$print_status(
         paste("block [", paste(value, collapse = ","), "] : "), 3
       )
     }
   },

   #' @field output \[`list()`, read-only\]\cr
   #' The output of the AO process, which is a \code{list} with the following
   #' elements:
   #'
   #' * \code{estimate} is the parameter vector at termination.
   #' * \code{value} is the function value at termination.
   #' * \code{details} is a `data.frame` with full information about the
   #'   AO process.
   #'   For each iteration (column `iteration`) it contains the function value
   #'   (column `value`), parameter values (columns starting with `p`
   #'   followed by the parameter index), the active parameter block
   #'   (columns starting with `b` followed by the parameter index, where
   #'   `1` stands for a parameter contained in the active parameter block
   #'   and `0` if not), and computation times in seconds (column `seconds`).
   #'   Only available if `add_details = TRUE`.
   #' * \code{seconds} is the overall computation time in seconds.
   #' * \code{stopping_reason} is a message explaining why the AO process
   #'   terminated.

   output = function(value) {
     if (missing(value)) {
       estimate <- self$get_parameter_best("full")
       c(
         list("estimate" = estimate),
         if (!is.null(private$.target)) {
           list("estimate_split" = split_by_target(
             estimate, private$.target, self$npar
           ))
         },
         list("value" = self$get_value_best()),
         if (self$add_details) list("details" = self$get_details()),
         list(
           "seconds" = self$get_seconds_total(),
           "stopping_reason" = private$.stopping_reason
         )
       )
     } else {
       cli::cli_abort(
         "Field {.var output} is read-only",
         call = NULL
       )
     }
   }
 ),

 private = list(

   ### inputs
   .target = character(),
   .npar = integer(),
   .partition = NULL,
   .new_block_probability = numeric(),
   .minimum_block_number = integer(),
   .verbose = logical(),
   .minimize = logical(),
   .iteration_limit = integer(),
   .seconds_limit = numeric(),
   .tolerance_value = numeric(),
   .tolerance_parameter = numeric(),
   .tolerance_history = integer(),

   ### must store 'tolerance_parameter_norm' inside list
   .tolerance_parameter_norm = list(function(x, y) sqrt(sum((x - y)^2))),

   ### results
   .iteration = 0L,
   .block = integer(),
   .details = data.frame(),
   .add_details = TRUE,
   .stopping_reason = "not terminated yet",

   ### error status
   .error = FALSE,
   .error_message = character()
 )
)
