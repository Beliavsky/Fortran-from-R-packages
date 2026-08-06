# Reparameterized global search (single search; the multi-start loop lives in
# rsdc_estimate() via control$n_starts).
#
# DEoptim explores a space in which EVERY point of the box is a valid model:
#  - correlations: canonical partial correlations per regime (Joe 2006, see
#    partial_corr.R) -- absorbs the positive-definiteness constraint;
#  - noX head, N >= 3: per-row softmax logits bounded at +/-4 -- absorbs the
#    row-simplex constraint AND floors every transition probability (an
#    unbounded softmax lets the search kill a regime, leaving the likelihood
#    flat in the dead regime's correlations and freezing the refinement);
#    +/-4 is the historical [0.01, 0.99] box mapped to the log-odds scale
#    (logit(0.99) = 4.6), but note the softmax image is not that box: with
#    N = 3 an individual entry can reach ~3e-4 while a stay probability
#    reaches ~0.98. The guarantee is positivity, not a uniform floor.
#    noX with N = 2 keeps its natural head: the two
#    stay probabilities are unconstrained jointly, the box is fully feasible;
#  - tvtp head: the usual beta box (the logistic/softmax link is always valid).
#
# The DE stage only selects a basin: the best few population members are
# refined by L-BFGS-B and the best refined point is returned. The top
# candidates are required to be mutually distant (the best DE point is not
# always in the best basin; near-ranked distant members often are).

# Internal constants, validated on the Fama-French K = 5 and K = 10 benchmarks
# (see the 1.7-0 NEWS entry). Deliberately not user-settable.
.RSDC_TOP_K     <- 3L     # refined candidates per search
.RSDC_DIV_DIST  <- 0.10   # min normalized RMS distance between candidates
.RSDC_LOGIT_BOX <- 4      # |softmax logit| bound for the noX (N >= 3) head
.RSDC_PENALTY_MIN <- 1e9  # objective values at/above this are the infeasibility
                          # penalty (1e10), not a likelihood

# L-BFGS-B returns convergence 52 ("ABNORMAL_TERMINATION_IN_LNSRCH") when its
# line search cannot improve on the point it started from. Since 1.7-0 the
# global stage already returns an L-BFGS-B-refined point, so the polish step
# below routinely starts at a converged optimum and 52 is the expected, benign
# outcome there. Relabel it as converged when the objective did not
# deteriorate, and leave every other code untouched.
#' @noRd
.rsdc_polish_code <- function(convergence, value, start_value) {
  conv <- as.integer(convergence)
  if (identical(conv, 52L) && is.finite(value) && is.finite(start_value) &&
      value <= start_value + 1e-8) 0L else conv
}
.RSDC_AGREE_TOL <- 2.5    # log-lik tolerance of the replication (agreement)
                          # check across independent searches; refits of the
                          # same optimum from different seeds scatter by
                          # ~1-2 log-points, so 1.0 would false-alarm

#' Head-block specification for the theta parameterization
#'
#' Returns the box and the map from the theta head block to the packed head
#' parameters expected by the model objective (transition probabilities for
#' noX, beta coefficients for tvtp, nothing for const).
#' @noRd
.rsdc_head_spec <- function(method, N, p = 0L) {
  if (method == "const")
    return(list(lo = numeric(0), hi = numeric(0),
                to_par = function(h) numeric(0)))
  if (method == "noX" && N == 2L)
    return(list(lo = rep(0.01, 2L), hi = rep(0.99, 2L), to_par = identity))
  if (method == "noX") {
    nb <- N * (N - 1L)
    return(list(
      lo = rep(-.RSDC_LOGIT_BOX, nb), hi = rep(.RSDC_LOGIT_BOX, nb),
      to_par = function(h) unlist(lapply(seq_len(N), function(i) {
        lg <- c(h[((i - 1L) * (N - 1L) + 1L):(i * (N - 1L))], 0)
        e <- exp(lg - max(lg))
        (e / sum(e))[seq_len(N - 1L)]
      }))))
  }
  # tvtp: same beta box as the natural parameterization
  nb <- if (N == 2L) N * p else N * (N - 1L) * p
  list(lo = rep(-10, nb), hi = rep(10, nb), to_par = identity)
}

#' Greedy max-diversity selection of refinement candidates
#'
#' Walks the population in objective order and keeps a member only if its
#' normalized RMS distance to every kept member exceeds `min_dist`;
#' backfills by rank if fewer than `k` qualify.
#' @noRd
.rsdc_select_diverse <- function(pop, vals, lower, upper, k, min_dist) {
  ord <- order(vals)
  width <- upper - lower
  sel <- ord[1L]
  for (i in ord[-1L]) {
    if (length(sel) >= k) break
    far <- all(vapply(sel, function(s)
      sqrt(mean(((pop[i, ] - pop[s, ]) / width)^2)) > min_dist, logical(1)))
    if (far) sel <- c(sel, i)
  }
  if (length(sel) < k) sel <- unique(c(sel, ord))[seq_len(min(k, nrow(pop)))]
  sel
}

#' One reparameterized global search: DEoptim in theta space + top-k refits
#'
#' @param objective Function of one packed natural parameter vector returning
#'   the value to minimize (the same closure the caller's own L-BFGS-B
#'   refinement uses, so search and refit optimize the identical criterion).
#' @param method "const", "noX" or "tvtp".
#' @param K,N,p Series count, regime count, covariate count (tvtp only).
#' @param lower,upper Natural-parameterization bounds (used for the refits;
#'   widened per candidate exactly like the warm-start path, since a softmax
#'   head can produce probabilities below the 0.01 box edge).
#' @param con The backend control list (itermax, NP, steptol, maxit, do_trace).
#' @returns list(par, value): the best refined natural parameter vector and
#'   its objective value.
#' @noRd
.rsdc_theta_search <- function(objective, method, K, N, p = 0L,
                               lower, upper, con) {
  C <- K * (K - 1L) / 2L
  n_reg <- max(N, 1L)
  lt <- lower.tri(diag(K))
  hs <- .rsdc_head_spec(method, N, p)
  n_head <- length(hs$lo)
  th_lo <- c(hs$lo, rep(-0.999, n_reg * C))
  th_hi <- c(hs$hi, rep(0.999, n_reg * C))
  d <- length(th_lo)

  theta_to_par <- function(theta) {
    zs <- if (n_head > 0L) theta[-seq_len(n_head)] else theta
    c(hs$to_par(theta[seq_len(n_head)]),
      unlist(lapply(seq_len(n_reg), function(r)
        .rsdc_pc_to_corr(zs[((r - 1L) * C + 1L):(r * C)], K)[lt])))
  }
  obj_theta <- function(theta) objective(theta_to_par(theta))

  de <- DEoptim::DEoptim(
    fn = obj_theta, lower = th_lo, upper = th_hi,
    control = DEoptim::DEoptim.control(
      itermax = con$itermax,
      NP = if (is.null(con$NP)) 10L * d else con$NP,
      strategy = 2, F = 0.8, CR = 0.9,
      trace = con$do_trace, reltol = 1e-8, steptol = con$steptol))

  pop <- de$member$pop
  vals <- apply(pop, 1L, obj_theta)
  sel <- .rsdc_select_diverse(pop, vals, th_lo, th_hi,
                              k = .RSDC_TOP_K, min_dist = .RSDC_DIV_DIST)

  refits <- lapply(sel, function(i) {
    par0 <- theta_to_par(pop[i, ])
    tryCatch(
      stats::optim(par0, objective, method = "L-BFGS-B",
                   lower = pmin(lower, par0), upper = pmax(upper, par0),
                   control = list(maxit = con$maxit)),
      error = function(e) NULL)
  })
  ok <- !vapply(refits, is.null, logical(1))
  if (!any(ok))   # degenerate corner: fall back to DE's raw best point
    return(list(par = theta_to_par(de$optim$bestmem), value = de$optim$bestval))
  refits <- refits[ok]
  best <- refits[[which.min(vapply(refits, `[[`, numeric(1), "value"))]]
  list(par = best$par, value = best$value)
}
