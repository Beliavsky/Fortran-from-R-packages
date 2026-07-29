# =============================================================================
# Geometry-adaptive sampling engine for bayesianOU (net-new, self-contained).
#
# A small R-native Hamiltonian Monte Carlo engine with a PLUGGABLE metric, used
# to probe and remedy hard posterior geometry in the nested OU model. The
# Level-2 latent block of ou_nested.stan was diagnosed (Sesion 9) as a
# multiplicative non-identification ridge (NCP funnel under an informative
# anchor); the centered reparametrization is the primary cure, and this engine
# is the complementary lever: a dense (constant) or Riemannian (position-
# dependent, SoftAbs) metric that navigates curvature Stan's diagonal NUTS
# cannot, plus a diagnostic of the realised energy mixing (E-BFMI).
#
# PROVENANCE: ported and adapted from the gdpar geometry engine (a sister
# research package, NOT production-ready and NOT a dependency). Only the lean
# subset bayesianOU needs is copied here -- the cmdstan target wrapper, the
# Euclidean (diagonal/dense) and Riemannian (SoftAbs / supplied-Fisher) metrics,
# the static HMC loop with the explicit and generalised-implicit leapfrog, and
# the E-BFMI. The advanced levels (relativistic, sub-Riemannian, GP-amortised
# Fisher, the adaptive controller and the synthetic suite) are deliberately
# NOT ported; they can be added if a future block needs them.
#
# Correctness vs efficiency: every metric is a PRECONDITIONER, not part of the
# target; the Metropolis correction uses the exact log-density, so the choice of
# geometry governs efficiency only, never the validity of the draws.
# =============================================================================

# ---- lean internal validators ---------------------------------------------
.ou_geom_pos_int <- function(x, nm) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 1 ||
      x != as.integer(x)) {
    stop(sprintf("'%s' must be a positive integer scalar.", nm), call. = FALSE)
  }
  as.integer(x)
}
.ou_geom_pos_num <- function(x, nm) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 0) {
    stop(sprintf("'%s' must be a non-negative numeric scalar.", nm),
         call. = FALSE)
  }
  as.numeric(x)
}

# ---------------------------------------------------------------------------
# Sampling target
# ---------------------------------------------------------------------------

#' Build a geometry-engine sampling target
#'
#' Wrap an unconstrained log-density (and its gradient, optionally its Hessian)
#' into the target object \code{\link{ou_geom_hmc}} integrates. Accepts either a
#' cmdstan fit/model compiled with \code{compile_model_methods = TRUE} (whose
#' \code{$log_prob} / \code{$grad_log_prob} / \code{$hessian} act on the
#' unconstrained scale) or explicit R closures.
#'
#' @param object Optional cmdstan fit/model (with methods) or a list carrying
#'   \code{log_prob}/\code{grad_log_prob}/\code{dim}.
#' @param log_prob,grad_log_prob R closures of the unconstrained parameter vector
#'   (used when \code{object} is \code{NULL}).
#' @param dim Integer unconstrained dimension (required for closures or a
#'   cmdstan model).
#' @param hessian Optional Hessian closure (used by the SoftAbs Riemannian
#'   metric); finite-differenced from the gradient when absent.
#' @param param_names Optional parameter names.
#' @param data Stan data list (used when \code{object} is an uncompiled-to-fit
#'   \code{CmdStanModel} that must be sampled once to expose its methods).
#'
#' @return An object of class \code{ou_geom_target}.
#' @export
ou_geom_target <- function(object = NULL, log_prob = NULL, grad_log_prob = NULL,
                           dim = NULL, hessian = NULL, param_names = NULL,
                           data = NULL) {
  if (!is.null(object)) {
    if (is.list(object) && is.function(object$log_prob) &&
        is.function(object$grad_log_prob)) {
      d <- dim %||% object$dim
      return(.ou_geom_target_obj(object$log_prob, object$grad_log_prob,
                                 hessian %||% object$hessian, d,
                                 param_names %||% object$param_names, "closure"))
    }
    if (is.function(object$grad_log_prob) ||
        inherits(object, c("CmdStanFit", "CmdStanMCMC", "CmdStanModel"))) {
      return(.ou_geom_target_cmdstan(object, dim, data, param_names))
    }
    stop(paste("Unrecognised 'object': supply a cmdstanr fit/model compiled",
               "with compile_model_methods = TRUE, or use the log_prob /",
               "grad_log_prob / dim arguments."), call. = FALSE)
  }
  if (!is.function(log_prob) || !is.function(grad_log_prob)) {
    stop("Both 'log_prob' and 'grad_log_prob' must be functions.", call. = FALSE)
  }
  if (is.null(dim)) stop("Argument 'dim' is required for a closure target.",
                         call. = FALSE)
  .ou_geom_target_obj(log_prob, grad_log_prob, hessian,
                      .ou_geom_pos_int(dim, "dim"), param_names, "closure")
}

.ou_geom_target_obj <- function(log_prob, grad_log_prob, hessian, dim,
                                param_names, backend) {
  obj <- list(log_prob = log_prob, grad_log_prob = grad_log_prob,
              hessian = hessian, dim = as.integer(dim),
              param_names = param_names %||% paste0("theta[", seq_len(dim), "]"),
              backend = backend)
  class(obj) <- c("ou_geom_target", "list")
  obj
}

.ou_geom_target_cmdstan <- function(object, dim, data, param_names) {
  fit <- object
  if (inherits(object, "CmdStanModel")) {
    fit <- object$sample(data = data, chains = 1, iter_warmup = 1,
                         iter_sampling = 1, refresh = 0, show_messages = FALSE,
                         show_exceptions = FALSE)
  }
  if (!is.function(fit$grad_log_prob)) {
    stop(paste("The cmdstan object does not expose grad_log_prob; compile the",
               "model with compile_model_methods = TRUE."), call. = FALSE)
  }
  if (is.null(dim)) stop("Argument 'dim' is required for a cmdstan target.",
                         call. = FALSE)
  lp <- function(theta) fit$log_prob(unconstrained_variables = theta)
  gl <- function(theta) as.numeric(fit$grad_log_prob(
    unconstrained_variables = theta))
  he <- if (is.function(fit$hessian)) {
    function(theta) fit$hessian(unconstrained_variables = theta)$hessian
  } else NULL
  .ou_geom_target_obj(lp, gl, he, .ou_geom_pos_int(dim, "dim"), param_names,
                      "cmdstan")
}

#' @export
print.ou_geom_target <- function(x, ...) {
  cat("<ou_geom_target> backend: ", x$backend, " | dim: ", x$dim, "\n", sep = "")
  cat("  hessian available: ", !is.null(x$hessian), "\n", sep = "")
  invisible(x)
}

# ---------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------

#' Euclidean (constant) metric for the geometry engine
#'
#' A position-independent mass matrix: the identity (diagonal) by default, or a
#' supplied symmetric positive-definite matrix / positive variance vector (dense
#' preconditioner; the remedy for a straight anisotropic canyon).
#'
#' @param dim Integer dimension (required when \code{M} is \code{NULL}).
#' @param M Optional \code{dim x dim} SPD matrix or length-\code{dim} positive
#'   diagonal vector.
#' @return A \code{ou_geom_metric} (position-independent).
#' @export
ou_geom_metric_euclidean <- function(dim = NULL, M = NULL) {
  if (is.null(M)) {
    if (is.null(dim)) stop("Supply either 'dim' or a mass matrix 'M'.",
                           call. = FALSE)
    M <- diag(.ou_geom_pos_int(dim, "dim"))
  } else if (is.vector(M) && !is.matrix(M)) {
    if (any(!is.finite(M)) || any(M <= 0)) {
      stop("Diagonal mass vector 'M' must be finite and positive.",
           call. = FALSE)
    }
    M <- diag(M, nrow = length(M))
  }
  M <- as.matrix(M)
  if (nrow(M) != ncol(M)) stop("Mass matrix 'M' must be square.", call. = FALSE)
  ch <- tryCatch(chol(M), error = function(e) NULL)
  if (is.null(ch)) stop("Mass matrix 'M' must be symmetric positive-definite.",
                        call. = FALSE)
  d <- nrow(M); Minv <- chol2inv(ch); L <- t(ch); ld <- 2 * sum(log(diag(ch)))
  obj <- list(position_dependent = FALSE, dim = d,
              mass = function(theta) M, inv_mass = function(theta) Minv,
              chol_mass = function(theta) L, logdet = function(theta) ld)
  class(obj) <- c("ou_geom_metric", "list")
  obj
}

#' Riemannian (position-dependent) metric for the geometry engine
#'
#' A position-dependent mass \eqn{M(\theta)} adapting the sampler's local
#' distance to the curvature of the log-posterior -- the remedy for a funnel
#' (variable curvature). \code{curvature = "softabs"} maps the eigenvalues of
#' the Hessian of \eqn{-\log\pi} through \eqn{\lambda\coth(\alpha\lambda)}
#' (Betancourt 2013), flooring nearly flat directions; it needs only the Hessian
#' (from the target or finite-differenced). \code{curvature = "fisher"} uses a
#' supplied expected-Fisher function (positive-definite by construction).
#'
#' @param target A \code{\link{ou_geom_target}}.
#' @param curvature \code{"softabs"} (default; observed Hessian) or
#'   \code{"fisher"} (supplied expected Fisher).
#' @param fisher,dfisher For \code{"fisher"}: \code{fisher(theta)} returning the
#'   expected Fisher matrix and optionally its derivative list.
#' @param alpha SoftAbs softening (larger tracks curvature more faithfully).
#' @param floor Minimum eigenvalue keeping \eqn{M(\theta)} positive-definite.
#' @param fd_step Finite-difference step for the Hessian / metric derivative.
#' @return A \code{ou_geom_metric} (position-dependent).
#' @export
ou_geom_metric_riemannian <- function(target, curvature = c("softabs", "fisher"),
                                       fisher = NULL, dfisher = NULL,
                                       alpha = 1e6, floor = 1e-8,
                                       fd_step = 1e-4) {
  curvature <- match.arg(curvature)
  if (!inherits(target, "ou_geom_target")) target <- ou_geom_target(target)
  d <- target$dim
  alpha <- .ou_geom_pos_num(alpha, "alpha"); floor <- .ou_geom_pos_num(floor, "floor")
  fd_step <- .ou_geom_pos_num(fd_step, "fd_step")
  if (identical(curvature, "fisher")) {
    if (!is.function(fisher)) {
      stop("curvature = 'fisher' requires a 'fisher' function.", call. = FALSE)
    }
    mass_fn <- function(theta) .ou_geom_floor_spd(as.matrix(fisher(theta)), floor)
    dmass_fn <- if (is.function(dfisher)) {
      function(theta) lapply(dfisher(theta), as.matrix)
    } else function(theta) .ou_geom_fd_dmass(mass_fn, theta, fd_step)
  } else {
    hess_lp <- .ou_geom_hessian_fn(target, fd_step)
    mass_fn <- function(theta) .ou_geom_softabs_mass(-hess_lp(theta), alpha, floor)
    dmass_fn <- function(theta)
      .ou_geom_softabs_dmass(theta, hess_lp, alpha, fd_step)
  }
  .ou_geom_metric_from_mass(mass_fn, dmass_fn, d, curvature, alpha)
}

.ou_geom_metric_from_mass <- function(mass_fn, dmass_fn, d, kind, alpha) {
  cache <- new.env(parent = emptyenv()); cache$theta <- NULL
  ensure <- function(theta) {
    if (!is.null(cache$theta) && identical(cache$theta, theta)) return(invisible())
    M <- mass_fn(theta); ch <- .ou_geom_chol_spd(M)
    cache$theta <- theta; cache$M <- M; cache$ch <- ch
    cache$Minv <- chol2inv(ch); cache$ld <- 2 * sum(log(diag(ch)))
    invisible()
  }
  obj <- list(
    position_dependent = TRUE, dim = d, metric_kind = kind, alpha = alpha,
    mass = function(theta) { ensure(theta); cache$M },
    inv_mass = function(theta) { ensure(theta); cache$Minv },
    chol_mass = function(theta) { ensure(theta); t(cache$ch) },
    logdet = function(theta) { ensure(theta); cache$ld },
    dmass = function(theta) dmass_fn(theta))
  class(obj) <- c("ou_geom_metric", "list")
  obj
}

.ou_geom_chol_spd <- function(M) {
  M <- (M + t(M)) / 2
  ch <- tryCatch(chol(M), error = function(e) NULL)
  if (!is.null(ch)) return(ch)
  jit <- 1e-8 * max(abs(diag(M)), 1); chol(M + diag(jit, nrow(M)))
}
.ou_geom_floor_spd <- function(M, floor) {
  M <- (M + t(M)) / 2; ev <- eigen(M, symmetric = TRUE)
  lam <- pmax(ev$values, floor); ev$vectors %*% (lam * t(ev$vectors))
}
.ou_geom_softabs_vals <- function(lambda, alpha) {
  u <- alpha * lambda
  ucoth <- ifelse(abs(u) < 1e-4, 1 + u^2 / 3, u / tanh(u)); ucoth / alpha
}
.ou_geom_softabs_deriv <- function(lambda, alpha) {
  u <- alpha * lambda; small <- abs(u) < 1e-3; out <- numeric(length(u))
  un <- u[!small]; out[!small] <- 1 / tanh(un) - un / sinh(un)^2
  out[small] <- (2 / 3) * u[small]; out
}
.ou_geom_softabs_mass <- function(HU, alpha, floor) {
  HU <- (HU + t(HU)) / 2; ev <- eigen(HU, symmetric = TRUE)
  sav <- pmax(.ou_geom_softabs_vals(ev$values, alpha), floor)
  ev$vectors %*% (sav * t(ev$vectors))
}
.ou_geom_softabs_dmass <- function(theta, hess_lp_fn, alpha, h) {
  HU <- -hess_lp_fn(theta); HU <- (HU + t(HU)) / 2
  ev <- eigen(HU, symmetric = TRUE)
  Q <- ev$vectors; lam <- ev$values; sav <- .ou_geom_softabs_vals(lam, alpha)
  d <- length(lam); R <- matrix(0, d, d)
  for (i in seq_len(d)) for (j in seq_len(d)) {
    R[i, j] <- if (abs(lam[i] - lam[j]) > 1e-8) {
      (sav[i] - sav[j]) / (lam[i] - lam[j])
    } else .ou_geom_softabs_deriv((lam[i] + lam[j]) / 2, alpha)
  }
  lapply(seq_len(d), function(k) {
    e <- numeric(d); e[k] <- h
    HUp <- -hess_lp_fn(theta + e); HUm <- -hess_lp_fn(theta - e)
    dHU <- ((HUp + t(HUp)) - (HUm + t(HUm))) / (4 * h)
    A <- crossprod(Q, dHU %*% Q); Q %*% (R * A) %*% t(Q)
  })
}
.ou_geom_hessian_fn <- function(target, h) {
  if (is.function(target$hessian)) return(target$hessian)
  function(theta) {
    d <- length(theta); H <- matrix(0, d, d)
    for (k in seq_len(d)) {
      e <- numeric(d); e[k] <- h
      H[, k] <- (target$grad_log_prob(theta + e) -
                   target$grad_log_prob(theta - e)) / (2 * h)
    }
    (H + t(H)) / 2
  }
}
.ou_geom_fd_dmass <- function(mass_fn, theta, h) {
  d <- length(theta)
  lapply(seq_len(d), function(k) {
    e <- numeric(d); e[k] <- h
    (mass_fn(theta + e) - mass_fn(theta - e)) / (2 * h)
  })
}

#' @export
print.ou_geom_metric <- function(x, ...) {
  cat("<ou_geom_metric> ",
      if (x$position_dependent) paste0("position-dependent (", x$metric_kind, ")")
      else "euclidean (constant)", " | dim: ", x$dim, "\n", sep = "")
  invisible(x)
}

# ---------------------------------------------------------------------------
# HMC engine
# ---------------------------------------------------------------------------

.ou_geom_kinetic_gaussian <- function(metric) {
  list(
    value = function(theta, p) {
      Minv <- metric$inv_mass(theta)
      0.5 * as.numeric(crossprod(p, Minv %*% p)) + 0.5 * metric$logdet(theta)
    },
    grad_p = function(theta, p) as.numeric(metric$inv_mass(theta) %*% p),
    grad_theta = function(theta, p) {
      if (!metric$position_dependent) return(rep(0, length(theta)))
      Minv <- metric$inv_mass(theta); dM <- metric$dmass(theta)
      Minv_p <- as.numeric(Minv %*% p)
      vapply(seq_along(theta), function(i) {
        dMi <- dM[[i]]
        0.5 * sum(Minv * dMi) - 0.5 * as.numeric(crossprod(Minv_p, dMi %*% Minv_p))
      }, numeric(1))
    },
    draw_momentum = function(theta) {
      L <- metric$chol_mass(theta); as.numeric(L %*% stats::rnorm(nrow(L)))
    })
}
.ou_geom_dH_dtheta <- function(target, kinetic, theta, p) {
  -target$grad_log_prob(theta) + kinetic$grad_theta(theta, p)
}
.ou_geom_leapfrog_step <- function(theta, p, target, metric, kinetic, eps,
                                   fp_tol = 1e-9, fp_max = 100L) {
  if (!metric$position_dependent) {
    p <- p + (eps / 2) * target$grad_log_prob(theta)
    theta <- theta + eps * kinetic$grad_p(theta, p)
    p <- p + (eps / 2) * target$grad_log_prob(theta)
    return(list(theta = theta, p = p, converged = TRUE))
  }
  tryCatch({
    p_half <- p; ok1 <- FALSE
    for (it in seq_len(fp_max)) {
      dH <- .ou_geom_dH_dtheta(target, kinetic, theta, p_half)
      p_new <- p - (eps / 2) * dH; delta <- max(abs(p_new - p_half))
      if (!is.finite(delta)) break
      p_half <- p_new; if (delta < fp_tol) { ok1 <- TRUE; break }
    }
    drift0 <- as.numeric(metric$inv_mass(theta) %*% p_half); ok2 <- FALSE
    theta_new <- theta
    for (it in seq_len(fp_max)) {
      th_new <- theta + (eps / 2) *
        (drift0 + as.numeric(metric$inv_mass(theta_new) %*% p_half))
      delta <- max(abs(th_new - theta_new))
      if (!is.finite(delta)) break
      theta_new <- th_new; if (delta < fp_tol) { ok2 <- TRUE; break }
    }
    p_final <- p_half - (eps / 2) *
      .ou_geom_dH_dtheta(target, kinetic, theta_new, p_half)
    list(theta = theta_new, p = p_final, converged = ok1 && ok2)
  }, error = function(e) list(theta = theta, p = p, converged = FALSE))
}
.ou_geom_leapfrog_traj <- function(theta, p, target, metric, kinetic, eps, L,
                                   fp_tol = 1e-9, fp_max = 100L) {
  conv <- TRUE
  for (i in seq_len(L)) {
    st <- .ou_geom_leapfrog_step(theta, p, target, metric, kinetic, eps,
                                 fp_tol, fp_max)
    theta <- st$theta; p <- st$p
    if (isFALSE(st$converged)) { conv <- FALSE; break }
  }
  list(theta = theta, p = p, converged = conv)
}
.ou_geom_hamiltonian <- function(target, kinetic, theta, p) {
  -target$log_prob(theta) + kinetic$value(theta, p)
}
.ou_geom_ebfmi <- function(energy) {
  energy <- energy[is.finite(energy)]
  if (length(energy) < 2L) return(NA_real_)
  den <- sum((energy - mean(energy))^2)
  if (!is.finite(den) || den <= 0) return(NA_real_)
  sum(diff(energy)^2) / den
}

#' Static Hamiltonian Monte Carlo with a pluggable geometry
#'
#' Fixed step-size / fixed trajectory-length HMC with a Metropolis correction
#' over a pluggable metric. The default Euclidean metric is textbook HMC; a
#' Riemannian metric reuses the same loop with the generalised implicit leapfrog
#' of Girolami & Calderhead (2011). The metric is a preconditioner only, so the
#' returned draws are exact regardless of which metric is chosen.
#'
#' @param target A \code{\link{ou_geom_target}} (or object accepted by it).
#' @param metric A \code{\link{ou_geom_metric_euclidean}} /
#'   \code{\link{ou_geom_metric_riemannian}}. Defaults to the identity metric.
#' @param epsilon Leapfrog step size.
#' @param L Leapfrog steps per proposal.
#' @param n_iter Retained iterations.
#' @param n_warmup Burn-in iterations discarded (no adaptation).
#' @param init Optional initial position (defaults to zeros).
#' @param seed Optional integer seed (RNG state set and restored).
#' @param fp_tol,fp_max Fixed-point tolerance / max iterations for the implicit
#'   leapfrog (position-dependent metric only).
#'
#' @return A list of class \code{ou_geom_hmc} with \code{draws}, \code{accept_rate},
#'   \code{n_divergent}, \code{energy}, \code{ebfmi}, \code{epsilon}, \code{L} and
#'   \code{metric_type}.
#' @export
ou_geom_hmc <- function(target, metric = NULL, epsilon = 0.1, L = 20L,
                        n_iter = 1000L, n_warmup = 500L, init = NULL,
                        seed = NULL, fp_tol = 1e-9, fp_max = 100L) {
  if (!inherits(target, "ou_geom_target")) target <- ou_geom_target(target)
  epsilon <- .ou_geom_pos_num(epsilon, "epsilon"); L <- .ou_geom_pos_int(L, "L")
  n_iter <- .ou_geom_pos_int(n_iter, "n_iter")
  if (!is.numeric(n_warmup) || length(n_warmup) != 1L || n_warmup < 0) {
    stop("Argument 'n_warmup' must be a non-negative integer scalar.",
         call. = FALSE)
  }
  d <- target$dim
  if (is.null(metric)) metric <- ou_geom_metric_euclidean(dim = d)
  if (!is.null(seed)) {
    seed <- .ou_geom_pos_int(seed, "seed")
    if (exists(".Random.seed", envir = .GlobalEnv)) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv)
      on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
    }
    set.seed(seed)
  }
  kinetic <- .ou_geom_kinetic_gaussian(metric)
  theta <- if (is.null(init)) rep(0, d) else as.numeric(init)
  if (length(theta) != d) stop(sprintf("'init' must have length %d.", d),
                               call. = FALSE)
  total <- as.integer(n_warmup) + as.integer(n_iter)
  draws <- matrix(NA_real_, nrow = n_iter, ncol = d)
  colnames(draws) <- target$param_names
  energy <- numeric(n_iter); n_accept <- 0L; n_div <- 0L; kept <- 0L
  for (it in seq_len(total)) {
    p0 <- kinetic$draw_momentum(theta)
    H0 <- .ou_geom_hamiltonian(target, kinetic, theta, p0)
    prop <- tryCatch(
      .ou_geom_leapfrog_traj(theta, p0, target, metric, kinetic, epsilon, L,
                             fp_tol, fp_max),
      error = function(e) list(theta = theta, p = p0, converged = FALSE))
    H1 <- tryCatch(.ou_geom_hamiltonian(target, kinetic, prop$theta, prop$p),
                   error = function(e) Inf)
    dH <- H1 - H0
    divergent <- isFALSE(prop$converged) || !is.finite(dH) || abs(dH) > 1000
    if (divergent) n_div <- n_div + 1L
    accept <- !divergent && (log(stats::runif(1)) < -dH)
    if (accept) { theta <- prop$theta; n_accept <- n_accept + 1L }
    if (it > n_warmup) { kept <- kept + 1L; draws[kept, ] <- theta; energy[kept] <- H0 }
  }
  obj <- list(draws = draws, accept_rate = n_accept / total, n_divergent = n_div,
              energy = energy, ebfmi = .ou_geom_ebfmi(energy),
              epsilon = epsilon, L = as.integer(L),
              metric_type = if (metric$position_dependent) "position_dependent"
                            else "euclidean_constant")
  class(obj) <- c("ou_geom_hmc", "list")
  obj
}

#' @export
print.ou_geom_hmc <- function(x, ...) {
  cat("<ou_geom_hmc> ", nrow(x$draws), " draws x ", ncol(x$draws), " dims\n",
      sep = "")
  cat("  metric: ", x$metric_type, " | epsilon: ", x$epsilon, " | L: ", x$L,
      "\n", sep = "")
  cat("  accept: ", format(x$accept_rate, digits = 3), " | divergent: ",
      x$n_divergent, " | E-BFMI: ", format(x$ebfmi, digits = 3), "\n", sep = "")
  invisible(x)
}

# ---------------------------------------------------------------------------
# Bridge: compiled-with-methods cmdstan model + data -> engine target
# ---------------------------------------------------------------------------

#' Bridge a methods-enabled cmdstan model to the geometry engine
#'
#' Decoupled core (the improvement over the gdpar bridge, which required a fitted
#' object of a specific class and recompiled from its source): takes a
#' \code{CmdStanModel} ALREADY compiled with \code{compile_model_methods = TRUE}
#' together with its Stan data, samples one iteration to expose the standalone
#' \code{log_prob}/\code{grad_log_prob}/\code{hessian} methods, reads the
#' unconstrained dimension and a posterior-mean warm-start, and returns an
#' \code{\link{ou_geom_target}} plus a reference position. Nothing in the regular
#' fit path is touched.
#'
#' @param model A \code{CmdStanModel} compiled with
#'   \code{compile_model_methods = TRUE}.
#' @param stan_data The Stan data list.
#' @param dim Optional unconstrained dimension (derived from the one-iteration
#'   fit when \code{NULL}).
#' @param reference Optional unconstrained warm-start (posterior mean of the
#'   one-iteration fit when \code{NULL}).
#' @param hessian Logical; request the standalone Hessian (best-effort: falls
#'   back to gradient-only when higher-order autodiff does not compile).
#' @param methods_seed Integer seed forwarded to \code{init_model_methods()}.
#'
#' @return A list with \code{target} (an \code{ou_geom_target}), \code{dim},
#'   \code{reference}, \code{fit} (the one-iteration cmdstan fit with methods) and
#'   \code{has_hessian}.
#' @export
ou_geom_bridge <- function(model, stan_data, dim = NULL, reference = NULL,
                           hessian = TRUE, methods_seed = 1L) {
  if (!inherits(model, "CmdStanModel")) {
    stop("'model' must be a CmdStanModel compiled with compile_model_methods = TRUE.",
         call. = FALSE)
  }
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("Package 'posterior' is required to unconstrain the draws for the bridge.",
         call. = FALSE)
  }
  fit <- model$sample(data = stan_data, chains = 1L, iter_warmup = 1L,
                      iter_sampling = 1L, refresh = 0L, show_messages = FALSE,
                      show_exceptions = FALSE)
  has_hessian <- tryCatch({
    fit$init_model_methods(seed = as.integer(methods_seed), verbose = FALSE,
                           hessian = isTRUE(hessian))
    isTRUE(hessian)
  }, error = function(e) {
    tryCatch({
      fit$init_model_methods(seed = as.integer(methods_seed), verbose = FALSE,
                             hessian = FALSE)
      FALSE
    }, error = function(e2)
      stop(sprintf("init_model_methods failed: %s", conditionMessage(e2)),
           call. = FALSE))
  })
  um <- as.matrix(posterior::as_draws_matrix(fit$unconstrain_draws()))
  d <- dim %||% ncol(um)
  ref <- reference %||% unname(colMeans(um))
  tgt <- .ou_geom_target_cmdstan(fit, d, stan_data, NULL)
  if (!has_hessian) tgt$hessian <- NULL
  list(target = tgt, dim = as.integer(d), reference = ref, fit = fit,
       has_hessian = has_hessian)
}
