#' SCIP solver control parameters
#'
#' Create a list of control parameters for the SCIP solver. Parameters are
#' organized into logical groups: output, limits, tolerances, presolving,
#' LP, branching, and heuristics. Any SCIP parameter can also be set directly
#' using its native path via \code{...}.
#'
#' @section Output:
#' \describe{
#'   \item{\code{verbose}}{Logical; print solver output. Default \code{TRUE}.}
#'   \item{\code{verbosity_level}}{Integer 0--5; verbosity detail (0 = none,
#'     3 = normal, 5 = full). Default \code{3}. Ignored if \code{verbose = FALSE}.}
#'   \item{\code{display_freq}}{Integer; display a status line every this many nodes
#'     (\code{-1} = never). Default \code{100}.}
#' }
#'
#' @section Termination Limits:
#' \describe{
#'   \item{\code{time_limit}}{Numeric; maximum solving time in seconds.
#'     Default \code{Inf} (no limit).}
#'   \item{\code{node_limit}}{Integer; maximum number of branch-and-bound nodes.
#'     Default \code{-1L} (no limit).}
#'   \item{\code{stall_node_limit}}{Integer; nodes without improvement before stopping.
#'     Default \code{-1L} (no limit).}
#'   \item{\code{sol_limit}}{Integer; stop after finding this many feasible solutions.
#'     Default \code{-1L} (no limit).}
#'   \item{\code{best_sol_limit}}{Integer; stop after this many improving solutions.
#'     Default \code{-1L} (no limit).}
#'   \item{\code{mem_limit}}{Numeric; memory limit in MB. Default \code{Inf}
#'     (no limit).}
#'   \item{\code{restart_limit}}{Integer; maximum restarts. Default \code{-1L}
#'     (no limit).}
#' }
#'
#' @section Tolerances:
#' \describe{
#'   \item{\code{gap_limit}}{Numeric; relative MIP gap tolerance; solver stops
#'     when the gap falls below this value. Default \code{0} (prove optimality).}
#'   \item{\code{abs_gap_limit}}{Numeric; absolute gap between primal and dual
#'     bound. Default \code{0}.}
#'   \item{\code{feastol}}{Numeric; feasibility tolerance for LP constraints.
#'     Default \code{1e-6}.}
#'   \item{\code{dualfeastol}}{Numeric; dual feasibility tolerance.
#'     Default \code{1e-7}.}
#'   \item{\code{epsilon}}{Numeric; absolute values below this are treated as zero.
#'     Default \code{1e-9}.}
#' }
#'
#' @section Presolving:
#' \describe{
#'   \item{\code{presolving}}{Logical; enable presolving. Default \code{TRUE}.}
#'   \item{\code{presolve_rounds}}{Integer; maximum presolving rounds
#'     (\code{-1} = unlimited). Default \code{-1L}.}
#' }
#'
#' @section LP:
#' \describe{
#'   \item{\code{lp_threads}}{Integer; number of threads for LP solver.
#'     Default \code{1L}.}
#'   \item{\code{lp_iteration_limit}}{Integer; LP iteration limit per solve
#'     (\code{-1} = no limit). Default \code{-1L}.}
#'   \item{\code{lp_scaling}}{Logical; enable LP scaling. Default \code{TRUE}.}
#' }
#'
#' @section Branching:
#' \describe{
#'   \item{\code{branching_score}}{Character; branching score function:
#'     \code{"s"} (sum), \code{"p"} (product), \code{"q"} (quotient).
#'     Default \code{"p"}.}
#' }
#'
#' @section Heuristics:
#' \describe{
#'   \item{\code{heuristics_emphasis}}{Character; heuristic emphasis setting:
#'     \code{"default"}, \code{"aggressive"}, \code{"fast"}, or \code{"off"}.
#'     Default \code{"default"}.}
#' }
#'
#' @section Parallel:
#' \describe{
#'   \item{\code{threads}}{Integer; number of threads for concurrent solving.
#'     Default \code{1L}. Note: concurrent solving may require a SCIP build
#'     compiled with parallel support (e.g., \code{PARASCIP=true}); not all
#'     installations provide this.}
#' }
#'
#' @param verbose Logical; print solver output. Default \code{TRUE}.
#' @param verbosity_level Integer 0--5; verbosity detail. Default \code{3}.
#' @param display_freq Integer; node display frequency. Default \code{100}.
#' @param time_limit Numeric; time limit in seconds. Default \code{Inf}.
#' @param node_limit Integer; max nodes. Default \code{-1L}.
#' @param stall_node_limit Integer; stall nodes. Default \code{-1L}.
#' @param sol_limit Integer; solution limit. Default \code{-1L}.
#' @param best_sol_limit Integer; improving solution limit. Default \code{-1L}.
#' @param mem_limit Numeric; memory limit in MB. Default \code{Inf}.
#' @param restart_limit Integer; restart limit. Default \code{-1L}.
#' @param gap_limit Numeric; relative MIP gap. Default \code{0}.
#' @param abs_gap_limit Numeric; absolute gap. Default \code{0}.
#' @param feastol Numeric; feasibility tolerance. Default \code{1e-6}.
#' @param dualfeastol Numeric; dual feasibility tolerance. Default \code{1e-7}.
#' @param epsilon Numeric; zero tolerance. Default \code{1e-9}.
#' @param presolving Logical; enable presolving. Default \code{TRUE}.
#' @param presolve_rounds Integer; presolve rounds. Default \code{-1L}.
#' @param lp_threads Integer; LP solver threads. Default \code{1L}.
#' @param lp_iteration_limit Integer; LP iteration limit. Default \code{-1L}.
#' @param lp_scaling Logical; LP scaling. Default \code{TRUE}.
#' @param branching_score Character; score function. Default \code{"p"}.
#' @param heuristics_emphasis Character; heuristic emphasis. Default \code{"default"}.
#' @param threads Integer; parallel solving threads. Default \code{1L}. See
#'   Parallel section for caveats.
#' @param ... Additional SCIP parameters as name-value pairs, using SCIP's
#'   native hierarchical parameter paths (e.g., \code{"lp/fastmip" = 1},
#'   \code{"conflict/enable" = FALSE}). See the
#'   \href{https://www.scipopt.org/doc/html/PARAMETERS.php}{SCIP parameter documentation}
#'   for the full list.
#' @return A named list of class \code{"scip_control"} with components:
#' \describe{
#'   \item{\code{verbose}}{Logical; verbosity flag.}
#'   \item{\code{scip_params}}{Named list; all parameters as SCIP native paths.}
#' }
#' @examples
#' ## Quick solve with 60-second time limit
#' ctrl <- scip_control(time_limit = 60)
#'
#' ## Quiet solve with 1% gap tolerance
#' ctrl <- scip_control(verbose = FALSE, gap_limit = 0.01)
#'
#' ## Aggressive heuristics, no presolving
#' ctrl <- scip_control(heuristics_emphasis = "aggressive", presolving = FALSE)
#'
#' ## Pass a native SCIP parameter directly
#' ctrl <- scip_control("conflict/enable" = FALSE, "separating/maxrounds" = 5L)
#' @seealso \code{\link{scip_solve}}, \code{\link{scip_set_param}}
#' @export
scip_control <- function(verbose = TRUE,
                         verbosity_level = 3L,
                         display_freq = 100L,
                         ## Limits
                         time_limit = Inf,
                         node_limit = -1L,
                         stall_node_limit = -1L,
                         sol_limit = -1L,
                         best_sol_limit = -1L,
                         mem_limit = Inf,
                         restart_limit = -1L,
                         ## Tolerances
                         gap_limit = 0,
                         abs_gap_limit = 0,
                         feastol = 1e-6,
                         dualfeastol = 1e-7,
                         epsilon = 1e-9,
                         ## Presolving
                         presolving = TRUE,
                         presolve_rounds = -1L,
                         ## LP
                         lp_threads = 1L,
                         lp_iteration_limit = -1L,
                         lp_scaling = TRUE,
                         ## Branching
                         branching_score = "p",
                         ## Heuristics
                         heuristics_emphasis = "default",
                         ## Parallel
                         threads = 1L,
                         ...) {

    ## -- Validate --
    verbose <- as.logical(verbose)[1L]
    verbosity_level <- as.integer(verbosity_level)[1L]
    display_freq <- as.integer(display_freq)[1L]
    time_limit <- as.double(time_limit)[1L]
    node_limit <- as.integer(node_limit)[1L]
    stall_node_limit <- as.integer(stall_node_limit)[1L]
    sol_limit <- as.integer(sol_limit)[1L]
    best_sol_limit <- as.integer(best_sol_limit)[1L]
    mem_limit <- as.double(mem_limit)[1L]
    restart_limit <- as.integer(restart_limit)[1L]
    gap_limit <- as.double(gap_limit)[1L]
    abs_gap_limit <- as.double(abs_gap_limit)[1L]
    feastol <- as.double(feastol)[1L]
    dualfeastol <- as.double(dualfeastol)[1L]
    epsilon <- as.double(epsilon)[1L]
    presolving <- as.logical(presolving)[1L]
    presolve_rounds <- as.integer(presolve_rounds)[1L]
    lp_threads <- as.integer(lp_threads)[1L]
    lp_iteration_limit <- as.integer(lp_iteration_limit)[1L]
    lp_scaling <- as.logical(lp_scaling)[1L]
    branching_score <- match.arg(branching_score, c("s", "p", "q"))
    heuristics_emphasis <- match.arg(heuristics_emphasis,
                                     c("default", "aggressive", "fast", "off"))
    threads <- as.integer(threads)[1L]

    ## -- Build the SCIP parameter list --
    ## Only include parameters that differ from SCIP defaults or were
    ## explicitly set by the user, to keep things clean.
    params <- list()

    ## Display / verbosity
    if (!verbose) {
        params[["display/verblevel"]] <- 0L
    } else if (verbosity_level != 3L) {
        params[["display/verblevel"]] <- verbosity_level
    }
    if (display_freq != 100L) {
        params[["display/freq"]] <- display_freq
    }

    ## Limits
    if (is.finite(time_limit)) {
        params[["limits/time"]] <- time_limit
    }
    if (node_limit > 0L) {
        params[["limits/nodes"]] <- node_limit
    }
    if (stall_node_limit > 0L) {
        params[["limits/stallnodes"]] <- stall_node_limit
    }
    if (sol_limit > 0L) {
        params[["limits/solutions"]] <- sol_limit
    }
    if (best_sol_limit > 0L) {
        params[["limits/bestsol"]] <- best_sol_limit
    }
    if (is.finite(mem_limit)) {
        params[["limits/memory"]] <- mem_limit
    }
    if (restart_limit > 0L) {
        params[["limits/restarts"]] <- restart_limit
    }

    ## Tolerances
    if (gap_limit > 0) {
        params[["limits/gap"]] <- gap_limit
    }
    if (abs_gap_limit > 0) {
        params[["limits/absgap"]] <- abs_gap_limit
    }
    if (feastol != 1e-6) {
        params[["numerics/feastol"]] <- feastol
    }
    if (dualfeastol != 1e-7) {
        params[["numerics/dualfeastol"]] <- dualfeastol
    }
    if (epsilon != 1e-9) {
        params[["numerics/epsilon"]] <- epsilon
    }

    ## Presolving
    if (!presolving) {
        params[["presolving/maxrounds"]] <- 0L
    } else if (presolve_rounds != -1L) {
        params[["presolving/maxrounds"]] <- presolve_rounds
    }

    ## LP
    if (lp_threads != 1L) {
        params[["lp/threads"]] <- lp_threads
    }
    if (lp_iteration_limit > 0L) {
        params[["lp/iterlim"]] <- lp_iteration_limit
    }
    if (!lp_scaling) {
        params[["lp/scaling"]] <- FALSE
    }

    ## Branching
    if (branching_score != "p") {
        params[["branching/scorefunc"]] <- branching_score
    }

    ## Parallel
    if (threads != 1L) {
        params[["parallel/maxnthreads"]] <- threads
    }

    ## Extra SCIP parameters from ...
    extra <- list(...)
    if (length(extra) > 0L) {
        nms <- names(extra)
        if (is.null(nms) || any(nms == "")) {
            stop("All extra parameters in '...' must be named with SCIP parameter paths")
        }
        for (nm in nms) {
            params[[nm]] <- extra[[nm]]
        }
    }

    ## Return structure
    ctrl <- list(
        verbose = verbose,
        scip_params = params
    )

    ## Store heuristics emphasis separately for the C layer
    if (heuristics_emphasis != "default") {
        ctrl$heuristics_emphasis <- heuristics_emphasis
    }

    structure(ctrl, class = "scip_control")
}

#' Print method for scip_control objects
#'
#' @param x A \code{scip_control} object.
#' @param ... Ignored.
#' @return Invisible \code{x}.
#' @export
print.scip_control <- function(x, ...) {
    cat("SCIP Control Parameters\n")
    cat("-----------------------\n")
    cat("verbose:", x$verbose, "\n")
    if (!is.null(x$heuristics_emphasis)) {
        cat("heuristics emphasis:", x$heuristics_emphasis, "\n")
    }
    params <- x$scip_params
    if (length(params) > 0L) {
        cat("\nSCIP parameters:\n")
        nms <- names(params)
        width <- max(nchar(nms))
        for (i in seq_along(params)) {
            cat(sprintf("  %-*s = %s\n", width, nms[i], format(params[[i]])))
        }
    } else {
        cat("\nAll SCIP parameters at defaults.\n")
    }
    invisible(x)
}
