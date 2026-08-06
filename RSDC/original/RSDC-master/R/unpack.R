# Internal helper: unpack a packed parameter vector into (beta, rho_matrix, P),
# mirroring the packing used by the estimators and rsdc_likelihood(). Shared by
# the band, multi-step-forecast and Viterbi routines added in 1.5-0.
#' @noRd
.rsdc_unpack_par <- function(par, method, X, K, N) {
  n_rho <- N * K * (K - 1) / 2
  if (method == "const") {
    return(list(beta = NULL,
                rho_matrix = matrix(par, nrow = 1),
                P = matrix(1, 1, 1)))
  }
  if (method == "noX") {
    n_p   <- N * (N - 1)
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
    return(list(beta = NULL, rho_matrix = rho_m, P = P))
  }
  # tvtp
  p_cov  <- ncol(X)
  n_beta <- if (N == 2) N * p_cov else N * (N - 1L) * p_cov
  list(beta       = matrix(par[seq_len(n_beta)], nrow = N, byrow = TRUE),
       rho_matrix = matrix(par[(n_beta + 1):(n_beta + n_rho)], nrow = N, byrow = TRUE),
       P          = NULL)
}

# Build the N x N transition matrix at a single covariate row x (length p) for a
# TVTP beta matrix, using the same logistic (N=2) / softmax (N>=3) link as the filter.
#' @noRd
.rsdc_P_at <- function(beta, x, N) {
  p <- length(x)
  P <- matrix(0, N, N)
  if (N == 2) {
    for (i in seq_len(N)) {
      pii <- stats::plogis(sum(x * beta[i, ]))
      P[i, i] <- pii
      P[i, -i] <- 1 - pii
    }
  } else {
    for (i in seq_len(N)) {
      logits <- c(vapply(seq_len(N - 1L), function(j)
        sum(x * beta[i, ((j - 1L) * p + 1L):(j * p)]), numeric(1)), 0)
      e <- exp(logits - max(logits))
      P[i, ] <- e / sum(e)
    }
  }
  P
}
