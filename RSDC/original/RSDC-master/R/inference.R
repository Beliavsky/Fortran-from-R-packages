# Inference helpers added in 1.4-0:
#  - per-observation scores (for OPG / sandwich standard errors, item 8)
#  - delta-method standard errors on the natural scale of transition
#    probabilities and regime correlations (item 7)
#  - Markov-chain regime diagnostics (item 5)

#' Element-wise sqrt that maps negative variances to NA (without NA-propagating
#' the whole vector, as \code{pmax(x, NA)} would). Used for standard errors.
#' @noRd
.rsdc_safe_sqrt <- function(v) {
  v[!is.finite(v) | v < 0] <- NA_real_
  sqrt(v)
}

#' Per-observation log-likelihood contributions at a parameter vector
#'
#' Mirrors the parameter packing of \code{rsdc_likelihood} but returns the vector
#' of per-time log-likelihood contributions (length \eqn{T}) used to build
#' outer-product-of-gradients and sandwich covariances.
#' @noRd
.rsdc_ll_t <- function(par, method, residuals, X, K, N) {
  if (method == "const") {
    R <- diag(K); R[lower.tri(R)] <- par; R[upper.tri(R)] <- t(R)[upper.tri(R)]
    R <- R + diag(1e-8, K)
    inv <- solve(R)
    ld  <- determinant(R, logarithm = TRUE)$modulus[1]
    quad <- rowSums((residuals %*% inv) * residuals)
    return(-0.5 * (ld + quad + K * log(2 * pi)))
  }
  exog <- if (method == "tvtp") X else NULL
  n_rho <- N * K * (K - 1) / 2
  if (is.null(exog)) {
    n_p <- N * (N - 1)
    trans <- par[seq_len(n_p)]
    rho_m <- matrix(par[(n_p + 1):(n_p + n_rho)], nrow = N, byrow = TRUE)
    P <- matrix(0, N, N)
    if (N == 2) {
      P <- matrix(c(trans[1], 1 - trans[1], 1 - trans[2], trans[2]), 2, byrow = TRUE)
    } else {
      for (i in seq_len(N)) {
        fr <- trans[((i - 1L) * (N - 1L) + 1L):(i * (N - 1L))]
        P[i, ] <- c(fr, 1 - sum(fr))
      }
    }
    beta <- NULL
  } else {
    n_beta <- if (N == 2) N * ncol(exog) else N * (N - 1L) * ncol(exog)
    beta  <- matrix(par[seq_len(n_beta)], nrow = N, byrow = TRUE)
    rho_m <- matrix(par[(n_beta + 1):(n_beta + n_rho)], nrow = N, byrow = TRUE)
    P <- NULL
  }
  res <- rsdc_hamilton(residuals, exog, beta, rho_m, K, N, P)
  if (is.null(res$loglik_t)) rep(NA_real_, nrow(residuals)) else res$loglik_t
}

#' T x npar matrix of per-observation scores (numerical gradient of \code{.rsdc_ll_t})
#' @importFrom numDeriv jacobian
#' @noRd
.rsdc_scores <- function(par, method, residuals, X, K, N) {
  tryCatch(
    numDeriv::jacobian(function(p) .rsdc_ll_t(p, method, residuals, X, K, N), par),
    error = function(e) NULL)
}

#' Robust variance-covariance from stored bread (Hessian) and scores
#' @noRd
.rsdc_robust_vcov <- function(object, type = c("hessian", "opg", "sandwich")) {
  type <- match.arg(type)
  if (type == "hessian") return(object$vcov)
  S <- object$scores
  if (is.null(S)) stop("Robust covariance needs per-observation scores; refit with ",
                       "control = list(compute_se = TRUE).")
  meat <- crossprod(S)                              # sum_t s_t s_t'
  if (type == "opg") {
    V <- tryCatch(solve(meat), error = function(e) NULL)
  } else {                                          # sandwich = bread %*% meat %*% bread
    bread <- object$vcov
    if (is.null(bread)) stop("Sandwich covariance needs the Hessian-based vcov.")
    V <- bread %*% meat %*% bread
  }
  if (!is.null(V)) {
    # object$vcov is NULL for OPG on a singular Hessian; fall back to coef names.
    dn <- if (!is.null(object$vcov)) dimnames(object$vcov) else
      list(names(object$coefficients), names(object$coefficients))
    dimnames(V) <- dn
  }
  V
}

#' Stationary (ergodic) distribution of a transition matrix P
#' @noRd
.rsdc_ergodic <- function(P) {
  N <- nrow(P)
  if (N == 1L) return(1)
  ev <- eigen(t(P))
  i  <- which.min(abs(ev$values - 1))
  pi <- Re(ev$vectors[, i]); pi <- pi / sum(pi)
  if (any(pi < -1e-8)) return(rep(NA_real_, N))
  pmax(pi, 0) / sum(pmax(pi, 0))
}

#' Regime diagnostics: stay probabilities, expected durations, ergodic distribution
#' @noRd
.rsdc_diagnostics <- function(object) {
  P <- object$transition_matrix
  N <- nrow(P)
  stay <- diag(P)
  data.frame(
    regime          = seq_len(N),
    stay_prob       = stay,
    exp_duration    = ifelse(stay < 1, 1 / (1 - stay), Inf),
    ergodic_prob    = .rsdc_ergodic(P),
    row.names = NULL)
}

#' Delta-method standard errors on the natural scale (transition probs + correlations)
#'
#' Correlations are parameters, so their SEs come directly from \code{vcov}. The
#' representative transition matrix (at the mean covariate for \code{"tvtp"}) is a
#' smooth function of the parameters; its SEs use a numerical Jacobian.
#' @importFrom numDeriv jacobian
#' @noRd
.rsdc_natural_se <- function(object) {
  V <- object$vcov
  if (is.null(V)) return(NULL)
  method <- object$method; N <- object$N; K <- object$K
  C <- K * (K - 1) / 2
  par <- object$par
  np <- length(par)

  # transition matrix as a function of the parameter vector
  P_of_par <- function(p) {
    if (method == "const") return(numeric(0))
    if (method == "noX") {
      trans <- p[seq_len(N * (N - 1))]
      P <- matrix(0, N, N)
      if (N == 2) P <- matrix(c(trans[1], 1 - trans[1], 1 - trans[2], trans[2]), 2, byrow = TRUE)
      else for (i in seq_len(N)) { fr <- trans[((i-1)*(N-1)+1):(i*(N-1))]; P[i, ] <- c(fr, 1 - sum(fr)) }
    } else {
      p_cov <- object$p; ax <- object$avg_X
      n_beta <- if (N == 2) N * p_cov else N * (N - 1L) * p_cov
      beta <- matrix(p[seq_len(n_beta)], nrow = N, byrow = TRUE)
      P <- matrix(0, N, N)
      for (i in seq_len(N)) {
        if (N == 2) { pii <- plogis(sum(ax * beta[i, ])); P[i, ] <- if (i == 1) c(pii, 1 - pii) else c(1 - pii, pii) }
        else {
          logits <- c(vapply(seq_len(N - 1L), function(j) sum(ax * beta[i, ((j-1)*p_cov+1):(j*p_cov)]), 0), 0)
          el <- exp(logits - max(logits)); P[i, ] <- el / sum(el)
        }
      }
    }
    as.vector(t(P))      # row-major vec of P
  }

  out <- list()
  # regime correlations (direct from vcov: rho are parameters, last N*C entries)
  rho_idx <- (np - N * C + 1):np
  rho_se  <- .rsdc_safe_sqrt(diag(V)[rho_idx])
  out$correlations <- data.frame(
    regime = rep(seq_len(N), each = C),
    estimate = par[rho_idx], se = rho_se, row.names = NULL)

  if (method != "const") {
    J <- tryCatch(numDeriv::jacobian(P_of_par, par), error = function(e) NULL)
    if (!is.null(J)) {
      Pvar <- diag(J %*% V %*% t(J))
      Pse  <- .rsdc_safe_sqrt(Pvar)
      Pest <- P_of_par(par)
      idx  <- as.vector(t(matrix(seq_len(N * N), N, N)))  # labels row-major
      lab  <- outer(seq_len(N), seq_len(N), function(i, j) paste0("p[", i, ",", j, "]"))
      out$transition <- data.frame(
        entry = as.vector(t(lab)), estimate = Pest, se = Pse, row.names = NULL)
    }
  }
  out
}
