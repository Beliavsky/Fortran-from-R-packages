# This file is part of PIQP-R. It is partly based on osqp-r
# (https://github.com/osqp/osqp-r) which is licensed
# under Apache License 2.0
#
# Copyright (c) 2023 piqp-r authors
# Copyright (c) 2019 Paul Goulart, Bartolomeo Stellato
#
# This source code is licensed under the BSD 2-Clause License found in the
# LICENSE file in the root directory of this source tree.

#' @importFrom Rcpp evalCpp
#' @importFrom S7 new_class new_property class_any class_list class_logical new_generic new_external_generic method

#' @title The PIQP Solver Model Class
#'
#' @description An S7 class wrapping the PIQP C++ Solver. Users will never
#'   need to directly create instances of this class and should use the
#'   more user-friendly functions [piqp()] and [solve_piqp()].
#' @param solver_ptr external pointer to the C++ solver object
#' @param dims a named list with elements `n`, `p`, and `m`
#' @param dense_backend logical flag indicating if the dense solver is used
#' @export
piqp_model <- S7::new_class("piqp_model",
  properties = list(
    solver_ptr = S7::new_property(S7::class_any),
    dims = S7::new_property(S7::class_list),
    dense_backend = S7::new_property(S7::class_logical)
  )
)

# External generics for base::solve and stats::update
base_solve <- S7::new_external_generic("base", "solve", "a")
stats_update <- S7::new_external_generic("stats", "update", "object")

S7::method(base_solve, piqp_model) <- function(a, b, ...) {
  .Call('_piqp_solve_model', PACKAGE = 'piqp', a@solver_ptr, a@dense_backend)
}

S7::method(stats_update, piqp_model) <- function(object, P = NULL, c = NULL, A = NULL, b = NULL, G = NULL,
                                                   h_l = NULL, h_u = NULL, x_l = NULL, x_u = NULL, ...) {
  dims <- object@dims
  n <- dims$n; p <- dims$p; m <- dims$m
  if (! (length(c) %in% c(0, n) &&
           length(b) %in% c(0, p) &&
           length(h_l) %in% c(0, m) &&
           length(h_u) %in% c(0, m) &&
           length(x_l) %in% c(0, n) &&
           length(x_u) %in% c(0, n)) ) {
    stop("Update parameters do not match original problem dimensions")
  }

  sparse_backend <- !object@dense_backend

  if (!is.null(P)) {
    if (NCOL(P) != n || NROW(P) != n) {
      stop("P dimension not matching original problem")
    }
    if (sparse_backend) P <- ensure_dgc_matrix(P)
  }

  if (!is.null(A)) {
    if (NCOL(A) != n || NROW(A) != p) {
      stop("A dimension not matching original problem")
    }
    if (sparse_backend) A <- ensure_dgc_matrix(A)
  }

  if (!is.null(G)) {
    if (NCOL(G) != n || NROW(G) != m) {
      stop("G dimension not matching original problem")
    }
    if (sparse_backend) G <- ensure_dgc_matrix(G)
  }

  if (sparse_backend) {
    .Call('_piqp_piqp_update_sparse', PACKAGE = 'piqp', object@solver_ptr, P, c, A, b, G, h_l, h_u, x_l, x_u)
  } else {
    .Call('_piqp_piqp_update_dense', PACKAGE = 'piqp', object@solver_ptr, P, c, A, b, G, h_l, h_u, x_l, x_u)
  }
  invisible(object)
}

#' Get settings of a PIQP model
#'
#' @param model A `piqp_model` object
#' @param ... not used
#' @return a list of current settings
#' @export
get_settings <- S7::new_generic("get_settings", "model")

S7::method(get_settings, piqp_model) <- function(model) {
  .Call('_piqp_get_settings', PACKAGE = 'piqp', model@solver_ptr, model@dense_backend)
}

#' Get dimensions of a PIQP model
#'
#' @param model A `piqp_model` object
#' @param ... not used
#' @return a list with named elements `n`, `p`, and `m`
#' @export
get_dims <- S7::new_generic("get_dims", "model")

S7::method(get_dims, piqp_model) <- function(model) {
  model@dims
}

#' Update settings of a PIQP model
#'
#' @param model A `piqp_model` object
#' @param ... not used
#' @return invisible NULL
#' @export
update_settings <- S7::new_generic("update_settings", "model")

S7::method(update_settings, piqp_model) <- function(model, new_settings = list()) {
  invisible(.Call('_piqp_update_settings', PACKAGE = 'piqp', model@solver_ptr, model@dense_backend, new_settings))
}

#' PIQP Solver object
#'
#' @importFrom Matrix sparseMatrix
#' @inheritParams solve_piqp
#' @return An S7 object of class "piqp_model" with methods
#'   `solve()`, `update()`, [get_settings()], [get_dims()], [update_settings()]
#'   which can be used to solve the problem with updated settings / parameters.
#' @seealso [solve_piqp()],  [piqp_settings()]
#' @section Usage:
#' \preformatted{model = piqp(P = NULL, c = NULL, A = NULL, b = NULL, G = NULL,
#'              h_l = NULL, h_u = NULL, x_l = NULL, x_u = NULL,
#'              settings = piqp_settings(),
#'              backend = c("auto", "sparse", "dense"))
#'
#' solve(model)
#' update(model, P = NULL, c = NULL, A = NULL, b = NULL, G = NULL,
#'        h_l = NULL, h_u = NULL, x_l = NULL, x_u = NULL)
#' get_settings(model)
#' get_dims(model)
#' update_settings(model, new_settings = piqp_settings())
#'
#' }
#' @details
#' Allows one to solve a parametric
#' problem with for example warm starts between updates of the parameter, c.f. the examples.
#' The object returned by \code{piqp} contains several methods which can be used to either update/get details of the
#' problem, modify the optimization settings or attempt to solve the problem.
#' @examples
#' ## example, adapted from PIQP documentation
#' library(piqp)
#' library(Matrix)
#'
#' P <- Matrix(c(6., 0.,
#'               0., 4.), 2, 2, sparse = TRUE)
#' c <- c(-1., -4.)
#' A <- Matrix(c(1., -2.), 1, 2, sparse = TRUE)
#' b <- c(1.)
#' G <- Matrix(c(1., 2., -1., 0.), 2, 2, sparse = TRUE)
#' h_u <- c(0.2, -1.)
#' x_l <- c(-1., -Inf)
#' x_u <- c(1., Inf)
#'
#' settings <- list(verbose = TRUE)
#'
#' model <- piqp(P, c, A, b, G, h_u = h_u, x_l = x_l, x_u = x_u, settings = settings)
#'
#' # Solve
#' res <- solve(model)
#' res$x
#'
#' # Define new data
#' A_new <- Matrix(c(1., -3.), 1, 2, sparse = TRUE)
#' h_u_new <- c(2., 1.)
#'
#' # Update model and solve again
#' update(model, A = A_new, h_u = h_u_new)
#' res <- solve(model)
#' res$x
#'
#' @export piqp
piqp <- function(P = NULL, c = NULL, A = NULL, b = NULL, G = NULL,
                 h_l = NULL, h_u = NULL, x_l = NULL, x_u = NULL,
                 settings = list(), backend = c("auto", "sparse", "dense")) {

  # match possible options
  backend <- match.arg(backend)

  sparse_backend <- (backend == "sparse") || inherits(P, "simple_triplet_matrix") ||
    inherits(A, "simple_triplet_matrix") || inherits(G, "simple_triplet_matrix") ||
    inherits(P, "sparseMatrix") || inherits(A, "sparseMatrix") || inherits(G, "sparseMatrix")

  if (is.null(P)) {
    n <- length(c)
  } else {
    n <- NCOL(P)
  }

  if (n == 0) {
    stop("At least one of P and c must be supplied")
  }

  if (!sparse_backend) { ## dense

    if (is.null(P)) {
      P <- matrix(0, n, n)
    }

    if (is.null(A)) {
      p <- 0
      A <- matrix(0, 0, n)
      b <- numeric(0)
    } else {
      p <- nrow(A)
      if (length(b) != p)
        stop(sprintf("b length %d must match A number of rows %d", length(b), p))
      if (NCOL(A) != n)
        stop(sprintf("A should have %d columns", n))
    }

    if (is.null(G)) {
      m <- 0
      G <- matrix(0, 0, n)
      h_l <- numeric(0)
      h_u <- numeric(0)
    } else {
      m <- nrow(G)
      if (!is.null(h_u) && length(h_u) != m)
        stop(sprintf("h_u length %d must match G number of rows %d", length(h_u), m))
      if (!is.null(h_l) && length(h_l) != m)
        stop(sprintf("h_l length %d must match G number of rows %d", length(h_l), m))
      if (NCOL(G) != n)
        stop(sprintf("G should have %d columns", n))
    }

  } else { ## sparse

    if (is.null(P)) {
      P <- Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), dims = c(n, n))
    } else {
      P <- ensure_dgc_matrix(P, n, n)
    }

    if (is.null(A)) {
      p <- 0
      A <- Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0), dims = c(0, n))
      b <- numeric(0)
    } else {
      p <- NROW(A)
      if (length(b) != p) stop(sprintf("b length %d must match A number of rows %d", length(b), p))
      if (NCOL(A) != n) stop(sprintf("A should have %d columns", n))
      A <- ensure_dgc_matrix(A, p, n)
    }

    if (is.null(G)) {
      m <- 0
      G <- Matrix::sparseMatrix(i = integer(), j = integer(0), x = numeric(0), dims = c(0, n))
      h_l <- numeric(0)
      h_u <- numeric(0)
    } else {
      m <- NROW(G)
      if (!is.null(h_u) && length(h_u) != m) stop(sprintf("h_u length %d must match G number of rows %d", length(h_u), m))
      if (!is.null(h_l) && length(h_l) != m) stop(sprintf("h_l length %d must match G number of rows %d", length(h_l), m))
      if (NCOL(G) != n) stop(sprintf("G should have %d columns", n))
      G <- ensure_dgc_matrix(G, m, n)
    }
  }
  if (is.null(c)) {
    c <- numeric(n)
  } else {
    if (length(c) != n) {
      stop(sprintf("P and c must have the same dimension %d", n))
    }
  }

  if (is.null(h_l)) {
    h_l <- rep(-Inf, m)
  }

  if (is.null(h_u)) {
    h_u <- rep(Inf, m)
  }

  if (is.null(x_l)) {
    x_l <- rep(-Inf, n)
  } else {
    if (length(x_l) != n) stop(sprintf("x_l length should be %d", n))
  }

  if (is.null(x_u)) {
    x_u <- rep(Inf, n)
  } else {
    if (length(x_u) != n) stop(sprintf("x_u length should be %d", n))
  }

  dense_backend <- !sparse_backend
  if (dense_backend) {
    solver_ptr <- .Call('_piqp_piqp_dense_setup', PACKAGE = 'piqp', P, c, A, b, G, h_l, h_u, x_l, x_u, settings)
  } else {
    solver_ptr <- .Call('_piqp_piqp_sparse_setup', PACKAGE = 'piqp', P, c, A, b, G, h_l, h_u, x_l, x_u, settings)
  }

  piqp_model(solver_ptr = solver_ptr,
             dims = list(n = n, p = p, m = m),
             dense_backend = dense_backend)
}


#' Settings parameters with default values and types in parenthesis
#'
#' @param rho_init Initial value for the primal proximal penalty parameter rho (default = 1e-6)
#' @param delta_init Initial value for the augmented lagrangian penalty parameter delta (default = 1e-4)
#' @param eps_abs Absolute tolerance (default = 1e-8)
#' @param eps_rel Relative tolerance (default = 1e-9)
#' @param check_duality_gap Check terminal criterion on duality gap (default = TRUE)
#' @param eps_duality_gap_abs Absolute tolerance on duality gap (default = 1e-8)
#' @param eps_duality_gap_rel Relative tolerance on duality gap (default = 1e-9)
#' @param infeasibility_threshold Threshold for infeasibility detection (default = 0.9)
#' @param reg_lower_limit Lower limit for regularization (default = 1e-10)
#' @param reg_finetune_lower_limit Fine tune lower limit regularization (default = 1e-13)
#' @param reg_finetune_primal_update_threshold Threshold of number of no primal updates to transition to fine tune mode (default = 7)
#' @param reg_finetune_dual_update_threshold Threshold of number of no dual updates to transition to fine tune mode (default = 7)
#' @param max_iter Maximum number of iterations (default = 250)
#' @param max_factor_retires Maximum number of factorization retires before failure (default = 10)
#' @param preconditioner_scale_cost Scale cost in Ruiz preconditioner (default = FALSE)
#' @param preconditioner_reuse_on_update Reuse preconditioner on problem update (default = FALSE)
#' @param preconditioner_iter Maximum of preconditioner iterations (default = 10)
#' @param tau Maximum interior point step length (default = 0.99)
#' @param iterative_refinement_always_enabled Always run iterative refinement and not only on factorization failure (default = FALSE)
#' @param iterative_refinement_eps_abs Iterative refinement absolute tolerance (default = 1e-12)
#' @param iterative_refinement_eps_rel Iterative refinement relative tolerance (default = 1e-12)
#' @param iterative_refinement_max_iter Maximum number of iterations for iterative refinement (default = 10)
#' @param iterative_refinement_min_improvement_rate Minimum improvement rate for iterative refinement (default = 5.0)
#' @param iterative_refinement_static_regularization_eps Static regularization for KKT system for iterative refinement (default = 1e-8)
#' @param iterative_refinement_static_regularization_rel Static regularization w.r.t. the maximum abs diagonal term of KKT system. (default = .Machine$double.eps^2)
#' @param verbose Verbose printing (default = FALSE)
#' @param compute_timings Measure timing information internally (default = FALSE)
#' @return a list containing the settings parameters.
#' @export piqp_settings
piqp_settings <- function(
                          ## Main algorithm settings
                          rho_init = 1e-6,
                          delta_init = 1e-4,
                          eps_abs = 1e-8,
                          eps_rel = 1e-9,
                          check_duality_gap = TRUE,
                          eps_duality_gap_abs = 1e-8,
                          eps_duality_gap_rel = 1e-9,
                          infeasibility_threshold = 0.9,
                          reg_lower_limit = 1e-10,
                          reg_finetune_lower_limit = 1e-13,
                          reg_finetune_primal_update_threshold = 7L,
                          reg_finetune_dual_update_threshold = 7L,
                          max_iter = 250L,
                          max_factor_retires = 10L,
                          preconditioner_scale_cost = FALSE,
                          preconditioner_reuse_on_update = FALSE,
                          preconditioner_iter = 10L,
                          tau = 0.99,
                          iterative_refinement_always_enabled = FALSE,
                          iterative_refinement_eps_abs = 1e-12,
                          iterative_refinement_eps_rel = 1e-12,
                          iterative_refinement_max_iter = 10L,
                          iterative_refinement_min_improvement_rate = 5.0,
                          iterative_refinement_static_regularization_eps = 1e-8,
                          iterative_refinement_static_regularization_rel = .Machine$double.eps^2,
                          verbose = FALSE,
                          compute_timings = FALSE) {

  params <- as.list(environment())

  bool_params <- c("check_duality_gap",
                   "preconditioner_scale_cost",
                   "preconditioner_reuse_on_update",
                   "iterative_refinement_always_enabled",
                   "verbose",
                   "compute_timings")

  int_params <- c("reg_finetune_primal_update_threshold",
                  "reg_finetune_dual_update_threshold",
                  "max_iter",
                  "max_factor_retires",
                  "preconditioner_iter",
                  "iterative_refinement_max_iter")

  if (any(sapply(params, length) != 1L)) stop("piqp_settings: arguments should be scalars!")
  if (any(unlist(params[int_params]) < 0)) stop("piqp_settings: integer arguments should be >= 0!")

  ## The rest
  float_params <- setdiff(names(params), c(bool_params, int_params))

  for (x in bool_params) {
    params[[x]] <- as.logical(params[[x]])
  }
  for (x in int_params) {
    params[[x]] <- as.integer(params[[x]])
  }

  for (x in float_params) {
    params[[x]] <- as.numeric(params[[x]])
  }
  params
}

#' Return the solver status description string
#' @param code a valid solver return code
#' @return a status description string
#' @examples
#' status_description(1) ## for solved problem
#' status_description(-1) ## for max iterations limit reached
#' @export status_description
status_description <- function(code) {
  ## Solver descriptions (using return codes in reverse order and offset by 11)
  ## https://predict-epfl.github.io/piqp/interfaces/c_cpp/getting_started_cpp
  desc <- c("Invalid settings were provided to the solver.", #(-10)
            "The problem is unsolved, i.e., `solve` was never called.", #(-9)
            "Numerical error occurred during solving.", #(-8),
            NA,
            NA,
            NA,
            NA,
            "The problem is dual infeasible.", #(-3)
            "The problem is primal infeasible.", #(-2)
            "Iteration limit was reached.",#(-1),
            NA,
            "Solver solved problem up to given tolerance." #(1)
            )
  desc[code + 11]
}
