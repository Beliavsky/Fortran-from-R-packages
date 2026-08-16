#' Retrieve Implemented Test Problems for the SOLNP Suite
#'
#' Returns a list (or a single object) of implemented test problems corresponding to a selected suite.
#' Problem functions must follow the naming convention \sQuote{problem_name_problem} and return a list
#' describing the optimization problem (e.g., objective, constraints, bounds).
#'
#' @param suite Character. The test suite to draw from. Must be one of \dQuote{Hock-Schittkowski} or \dQuote{Other}.
#'              Default is \dQuote{Hock-Schittkowski}.
#' @param number Integer or vector of integers. One or more problem numbers to retrieve. Ignored if return_all = TRUE.
#' @param return_all Logical. If TRUE, returns all implemented problems in the specified suite. Default is FALSE.
#'
#' @return If one problem is requested and implemented, the evaluated problem object is returned directly.
#'         Otherwise, an unnamed list of evaluated problem objects is returned.
#'
#' @details
#' - Problems are matched by number within the selected suite, using the table from [solnp_problems_table()].
#' - If a requested problem is valid but not yet implemented (i.e., the corresponding function does not exist),
#'   a message will inform the user.
#' - If a problem number exceeds the allowable range (e.g., > 306 for Hock-Schittkowski), an error is raised.
#'
#' @examples
#' \dontrun{
#' # Retrieve a single HS problem
#' prob <- solnp_problem_suite(number = 1)
#'
#' # Retrieve multiple HS problems
#' probs <- solnp_problem_suite(number = c(1, 2, 3))
#'
#' # Retrieve problem in "Other" suite
#' other_prob <- solnp_problem_suite(suite = "Other", number = 1)
#' }
#'
#' @seealso [solnp_problems_table()]
#' @export
solnp_problem_suite <- function(suite = "Hock-Schittkowski", number = 1, return_all = FALSE)
{
    suite <- match.arg(suite[1], c("Hock-Schittkowski","Other"))
    v_problems <- solnp_problems_table()
    v_problems <- v_problems[v_problems$Suite == suite,]
    if (!return_all) {
        number <- unique(as.integer(number))
        if (any(number < 1)) stop("Problem number must be strictly positive.")
        if (suite == "Hock-Schittkowski" & any(number > 306)) stop("Some requested problem numbers do not exist (maximum is 306 for Hock-Schittkowski).")
        known <- v_problems$Number
        not_implemented <- setdiff(number, known)
        implemented <- intersect(number, known)

        if (length(not_implemented) > 0) {
            message("The following problems are not yet implemented: ",
                    paste(not_implemented, collapse = ", "))
        }
        funs <- lapply(implemented, function(n) {
            fname <- paste0(v_problems$Problem[v_problems$Number == n], "_problem")
            if (!exists(fname, mode = "function", inherits = TRUE)) {
                warning(sprintf("Function '%s' not found in the environment.", fname))
                return(NULL)
            }
            get(fname, mode = "function", inherits = TRUE)()
        })
        funs <- funs[!sapply(funs, is.null)]
        if (length(funs) == 1) {
            return(funs[[1]])
        } else {
            return(unname(funs))
        }
    } else {
        # return_all = TRUE: Return all implemented hsXX_problem functions
        all_funs <- lapply(seq_len(nrow(v_problems)), function(i) {
            fname <- paste0(v_problems$Problem[i], "_problem")
            if (exists(fname, mode = "function", inherits = TRUE)) {
                get(fname, mode = "function", inherits = TRUE)()
            } else {
                NULL
            }
        })
        all_funs <- all_funs[!sapply(all_funs, is.null)]
        return(unname(all_funs))
    }
}

#' List of Valid Test Problems for the SOLNP Suite
#'
#' Returns a data.frame of known and registered test problems used with the SOLNP solver.
#' The list includes problems from the Hock-Schittkowski suite as well as a selection
#' of other classic optimization problems.
#'
#' @return A data.frame with the following columns:
#' \describe{
#'   \item{Suite}{A character string indicating the suite the problem belongs to.
#'                One of \dQuote{Hock-Schittkowski} or \dQuote{Other}.}
#'   \item{Problem}{The base name of the problem function (without the \sQuote{_problem} suffix).}
#'   \item{Number}{An integer identifier used to index or request problems programmatically.}
#' }
#'
#' @details
#' - All problem functions are expected to follow the naming convention
#'   \sQuote{Problem_problem} (e.g., \sQuote{hs01_problem}).
#' - For Hock-Schittkowski problems, numbers range from 1 to 64, with a few
#'   selected extras (e.g., 110, 118, 119).
#' - The \dQuote{Other} suite includes named problems like \sQuote{box},
#'   \sQuote{alkylation}, \sQuote{entropy}, \sQuote{garch}, etc.,
#' and are numbered sequentially.
#'
#' @examples
#' # View all known problems
#' tail(solnp_problems_table())
#'
#' # Filter only HS problems
#' head(subset(solnp_problems_table(), Suite == "Hock-Schittkowski"))
#'
#' @seealso [solnp_problem_suite()]
#' @export
solnp_problems_table <- function()
{
    hs_problems <- c("hs01","hs02","hs03","hs04","hs05","hs06","hs07",
                     "hs08","hs09", paste0("hs",10:64),"hs110","hs118","hs119")
    hs_number <- c(1:64,110,118,119)
    other_problems <- sort(c("alkylation","wright4","wright9","box","entropy","rosen_suzuki","powell",
                             "himmelblau5","garch"))
    other_number <- c(1:length(other_problems))
    rbind(data.frame("Suite" = "Hock-Schittkowski", "Problem" = hs_problems, "Number" = hs_number),
          data.frame("Suite" = "Other", "Problem" = other_problems, "Number" = other_number))
}

#' Standardize an Optimization Problem to NLP Standard Form
#'
#' Converts a problem specified with two-sided inequalities and nonzero equality right-hand sides
#' to the standard nonlinear programming (NLP) form.
#'
#'
#' @param prob A list specifying the problem in SOLNP-compatible format, with components
#'   \code{fn}, \code{eq_fn}, \code{eq_jac}, \code{eq_b}, \code{ineq_fn}, \code{ineq_jac},
#'   \code{ineq_lower}, \code{ineq_upper}, and others.
#'
#' @return A list with the same structure as the input, but with \code{eq_fn} and \code{ineq_fn}
#' standardized to the forms \eqn{e(x) = 0} and \eqn{g(x) \leq 0}, and with
#' \code{eq_b}, \code{ineq_lower}, and \code{ineq_upper} removed.
#'
#' @details
#' The standard form given by the following set of equations:
#'
#' \deqn{ \min_x\ f(x) }
#' \deqn{ \textrm{subject to}\quad e(x) = 0 }
#' \deqn{ \qquad\qquad\qquad g(x) \leq 0 }
#'
#' Specifically:
#' \itemize{
#'   \item All equality constraints are standardized to \eqn{e(x) = e(x) - b = 0}
#'   \item Each two-sided inequality \eqn{l \leq g(x) \leq u} is converted to one or two
#'   one-sided constraints: \eqn{l - g(x) \leq 0}, \eqn{g(x) - u \leq 0}
#' }
#'
#' The returned problem object has all equalities as \eqn{e(x) = 0}, all inequalities as \eqn{g(x) \leq 0},
#' and any right-hand side or bounds are absorbed into the standardized constraint functions.
#'
#' @examples
#' # Alkylation problem
#' p <- solnp_problem_suite(suite = "Other", number = 1)
#' ps <- solnp_standardize_problem(p)
#' ps$eq_fn(ps$start)    # standardized equalities: e(x) = 0
#' ps$ineq_fn(ps$start)  # standardized inequalities: g(x) <= 0
#'
#' @seealso \code{\link{solnp_problem_suite}}
#'
#' @export
solnp_standardize_problem <- function(prob) {
    # Objective and gradient
    fn <- prob$fn
    gr <- prob$gr

    # Equalities
    eq_fn <- prob$eq_fn
    eq_jac <- prob$eq_jac
    eq_b <- if (!is.null(prob$eq_b)) prob$eq_b else NULL

    # Standardize equalities: e(x) = eq_fn(x) - eq_b = 0
    new_eq_fn <- if (is.null(eq_fn)) {
        NULL
    } else if (is.null(eq_b) || all(eq_b == 0)) {
        eq_fn
    } else {
        function(x) eq_fn(x) - eq_b
    }

    new_eq_jac <- eq_jac  # No transformation needed if analytic jacobian is for eq_fn

    # Inequalities
    ineq_fn <- prob$ineq_fn
    ineq_jac <- prob$ineq_jac
    ineq_lower <- prob$ineq_lower
    ineq_upper <- prob$ineq_upper

    # Standardize inequalities to g(x) <= 0 form
    new_ineq_fn <- if (is.null(ineq_fn)) {
        NULL
    } else {
        function(x) {
            gx <- ineq_fn(x)
            g <- c()
            if (!is.null(ineq_lower)) {
                g <- c(g, ineq_lower - gx)
            }
            if (!is.null(ineq_upper)) {
                g <- c(g, gx - ineq_upper)
            }
            if (length(g) == 0) NULL else g
        }
    }

    # Standardize inequality jacobian if supplied
    new_ineq_jac <- if (is.null(ineq_jac)) {
        NULL
    } else {
        function(x) {
            J <- ineq_jac(x)
            Jlist <- list()
            if (!is.null(ineq_lower)) {
                for (i in seq_along(ineq_lower)) {
                    Jlist[[length(Jlist) + 1]] <- -J[i, , drop = FALSE]
                }
            }
            if (!is.null(ineq_upper)) {
                for (i in seq_along(ineq_upper)) {
                    Jlist[[length(Jlist) + 1]] <- J[i, , drop = FALSE]
                }
            }
            if (length(Jlist) == 0) NULL else do.call(rbind, Jlist)
        }
    }

    # Return standardized problem
    prob_std <- prob
    prob_std$eq_fn <- new_eq_fn
    prob_std$eq_jac <- new_eq_jac
    prob_std$ineq_fn <- new_ineq_fn
    prob_std$ineq_jac <- new_ineq_jac
    prob_std$eq_b <- NULL      # standardized out
    prob_std$ineq_lower <- NULL
    prob_std$ineq_upper <- NULL
    prob_std
}
