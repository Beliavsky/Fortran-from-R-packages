#' Subseries-Based Cauchy Combination Test (SCT) for Spanning Hypotheses
#'
#' Computes robust p-values for testing spanning-related linear restrictions using the residual-based
#' subseries method. This function supports high-dimensional inference on hypotheses such as
#' \eqn{H_0^\delta}, \eqn{H_0^\alpha}, and their joint null \eqn{H_0^{\alpha, \delta}} using
#' a simulation-based approximation and aggregation via the Cauchy Combination Test (CCT).
#'
#' @param u A numeric vector of outcomes (e.g., returns or residuals) of length \eqn{T}.
#' @param x A numeric matrix of regressors (T x K), where the first column is used as baseline and others are compared.
#' @param ks A numeric vector of subseries exponents (e.g., \code{1/3}); each value defines subseries size as \code{floor(T^k)}. Default is \code{c(1/3)}.
#' @param L A numeric vector controlling the strength of randomization applied to residual scores. Default is \code{c(0, 2)}.
#'
#' @return A named numeric vector of CCT p-values. Each name encodes the test type, L value, and subseries size:
#' \describe{
#'   \item{CCTd}{Test of \eqn{H_0^\delta}: no directional (slope) deviation.}
#'   \item{CCTad}{Joint test of \eqn{H_0^{\alpha,\delta}}: no intercept or slope deviation.}
#'   \item{CCTa}{Test of \eqn{H_0^\alpha}: no intercept deviation.}
#' }
#' Names are suffixed with the values of \code{L} and the subseries exponent index, e.g., \code{CCTa_L2_k1}.
#'
#' @details
#' This function builds score vectors from OLS residuals and applies randomized weightings
#' (via \code{f_prods}) to simulate perturbations. These perturbed scores are passed through a subseries
#' t-test pipeline (via \code{f_testbm}), and resulting p-values are aggregated using the
#' Cauchy Combination Test (Liu & Xie, 2020), which is valid under arbitrary dependence.
#'
#' The methodology is motivated by Carlstein's (1986) subseries inference and extends it using
#' residual-based oracle approximations that remain valid when K increases with T.
#'
#' @references
#' \insertRef{ArdiaSessinou2025}{spantest} \cr
#'
#' @importFrom stats cov na.omit pf pnorm pt qnorm rnorm runif
#'
#' @keywords internal
#'
#' @noRd
#'
f_getpv <- function(u, x, ks = c(1/3), L = c(0, 2), B = 1L, seed = 123L) {
  one <- rep(1, nrow(x))
  y <- u - x[, 1]

  if (ncol(x) == 1) {
    x_centered <- matrix(nrow = nrow(x), ncol = 0)  # empty
    x_augmented <- x  # only one column
  } else {
    x_centered <- sweep(x[, -1, drop = FALSE], 1, x[, 1], "-")
    x_augmented <- cbind(x[, 1], x_centered)
  }

  x <- x_augmented

  X_main <- cbind(1, x)
  qr_main <- qr(X_main)
  res <- qr.resid(qr_main, y)

  xy_centered <- cbind(y, x_centered)

  f_compute_resid <- function(X, y) {
    qr <- qr(X)
    qr.resid(qr, y)
  }

  ew <- f_compute_resid(cbind(1, xy_centered), x[, 1])
  ew1 <- f_compute_resid(cbind(y, x), one)

  res_ew <- res * ew
  score <- score2 <- matrix(res_ew, ncol = 1)

  # Score processing: randomized perturbation then subseries Cauchy p-value.
  # With B > 1 the perturbation is drawn B times and the resulting p-values are
  # Cauchy-merged, so the reported value does not hinge on one arbitrary draw
  # (see the note in f_getpv_batch). B = 1 is the single-draw behaviour.
  f_process_scores <- function(score_mat, sn) {
    nB <- if (sn > 0L) max(1L, as.integer(B)) else 1L
    M  <- matrix(NA_real_, length(ks), nB)
    for (b in seq_len(nB)) {
      scoreb <- score_mat * f_prods(score_mat, sn, cseed = seed + b - 1L)
      M[, b] <- vapply(ks, function(k) f_testbm(scoreb, k = k)[1], numeric(1))
    }
    if (nB == 1L) return(M[, 1L])
    apply(M, 1L, function(p) {
      p <- stats::na.omit(p)
      if (!length(p)) NA_real_ else f_cauchypv(p)
    })
  }

  # Delta = 0
  test1 <- lapply(L, function(sn_value) {
    f_process_scores(score, sn = sn_value)
  })

  # Delta = alpha = 0 (joint)
  score_combo <- cbind(score2, res * ew1)
  test2 <- lapply(L, function(sn_value) {
    f_process_scores(score_combo, sn = sn_value)
  })

  # Alpha = 0
  score_alpha <- matrix(res * ew1, ncol = 1)
  test4 <- lapply(L, function(sn_value) {
    f_process_scores(score_alpha, sn = sn_value)
  })

  result <- unlist(c(test1, test2, test4))

  # Dynamic naming: CCT{d,ad,a} x L x k
  test_names <- paste0(
    rep(c("CCTd", "CCTad", "CCTa"),
        each = length(L) * length(ks)),
    "_L", rep(L,
              each = length(ks),
              times = length(c("CCTd", "CCTad", "CCTa"))),
    "_k", rep(seq_along(ks),
              times = length(L) * length(c("CCTd", "CCTad", "CCTa")))
  )

  names(result) <- test_names
  return(result)
}

#' Vectorized per-asset CCT p-values (batch over the test cross-section)
#'
#' Computes the same per-asset subseries CCT p-values as \code{f_getpv()} but for
#' all \eqn{N} test assets at once. Because every test asset shares the same
#' benchmark design, the benchmark-only QR decompositions are formed once and the
#' three swap regressions are obtained by Frisch--Waugh partialling, turning the
#' original \eqn{O(N)} loop of full QR factorizations into a handful of shared
#' factorizations plus vectorized arithmetic. It is numerically equivalent to
#' looping \code{f_getpv()} over the columns of \code{test}.
#'
#' @param bench Numeric \eqn{T \times K} matrix of benchmark returns.
#' @param test  Numeric \eqn{T \times N} matrix of test-asset returns.
#' @param ks,L  As in \code{f_getpv()}.
#' @param B     Number of independent perturbation draws merged per asset when
#'   \code{L > 0} (Cauchy rule). \code{B = 1} is the single-draw behaviour.
#' @param seed  Seed of the first draw; draw \eqn{b} uses \code{seed + b - 1}.
#'
#' @return A named list; each element is a length-\eqn{N} vector of per-asset
#' p-values, named as \code{CCT{d,ad,a}_L{L}_k{i}}.
#'
#' @keywords internal
#'
#' @noRd
#'
f_getpv_batch <- function(bench, test, ks = c(1/3), L = c(0, 2),
                          B = 1L, seed = 123L) {
  x  <- bench
  Tn <- nrow(x)
  K  <- ncol(x)
  N  <- ncol(test)
  x1 <- x[, 1]
  one <- rep(1, Tn)
  Y  <- test - x1                                   # T x N: y_j = test_j - x1

  if (K == 1) {
    xc   <- matrix(nrow = Tn, ncol = 0)
    xaug <- x
  } else {
    xc   <- sweep(x[, -1, drop = FALSE], 1, x1, "-")
    xaug <- cbind(x1, xc)
  }

  # Benchmark-only QR factorizations (shared across all test assets)
  qr_main <- qr(cbind(1, xaug))   # residualize on [1, x]
  qr_bd   <- qr(cbind(1, xc))     # FWL base for ew  (add y as the extra column)
  qr_ba   <- qr(xaug)             # FWL base for ew1 (add y as the extra column)

  # res_j = residual of y_j on [1, x]  (batched over all assets)
  res <- qr.resid(qr_main, Y)

  # ew_j  = residual of x1 on [1, xc, y_j]  via Frisch-Waugh on base [1, xc]
  x1_p <- qr.resid(qr_bd, x1)
  yp1  <- qr.resid(qr_bd, Y)
  b1   <- colSums(yp1 * x1_p) / colSums(yp1 * yp1)
  ew   <- x1_p - sweep(yp1, 2, b1, "*")

  # ew1_j = residual of 1 on [x, y_j]  via Frisch-Waugh on base x
  one_p <- qr.resid(qr_ba, one)
  yp2   <- qr.resid(qr_ba, Y)
  b2    <- colSums(yp2 * one_p) / colSums(yp2 * yp2)
  ew1   <- one_p - sweep(yp2, 2, b2, "*")

  D <- res * ew    # delta scores, T x N
  A <- res * ew1   # alpha scores, T x N

  # For each (L, k): student p-values for every asset via one f_ttest per score
  # matrix (columns are independent), then Cauchy-combine the two for CCTad.
  #
  # DE-RANDOMISATION (B). For L > 0 the score is multiplied by a random weight
  # vector w. That vector is SHARED by every test asset, so its effect does not
  # average out across the cross-section: a single draw is a single realisation
  # of the test, and with few subseries (short T) the resulting p-value can swing
  # over orders of magnitude from one draw to the next. With B > 1 we draw B
  # independent weight vectors and merge the B p-values per asset with the Cauchy
  # rule -- valid under arbitrary dependence, and self-correcting in that one
  # lucky draw among B does not carry the merged value. B = 1 reproduces the
  # single-draw behaviour exactly.
  out <- list()
  for (sn in L) {
    nB <- if (sn > 0L) max(1L, as.integer(B)) else 1L   # L = 0: nothing to average
    accD <- accA <- accAD <- vector("list", length(ks))
    for (ki in seq_along(ks))
      accD[[ki]] <- accA[[ki]] <- accAD[[ki]] <- matrix(NA_real_, N, nB)

    for (b in seq_len(nB)) {
      w  <- f_prods(matrix(0, Tn, 1), sn, cseed = seed + b - 1L)
      Dw <- D * w
      Aw <- A * w
      for (ki in seq_along(ks)) {
        sD  <- f_ttest(Dw, k = ks[ki])$student
        sA  <- f_ttest(Aw, k = ks[ki])$student
        accD[[ki]][, b]  <- sD
        accA[[ki]][, b]  <- sA
        accAD[[ki]][, b] <- vapply(seq_len(N),
                                   function(j) f_cauchypv(c(sD[j], sA[j])), numeric(1))
      }
    }

    merge_draws <- function(M) {
      if (ncol(M) == 1L) return(M[, 1L])
      apply(M, 1L, function(p) {
        p <- stats::na.omit(p)
        if (!length(p)) NA_real_ else f_cauchypv(p)
      })
    }
    for (ki in seq_along(ks)) {
      out[[paste0("CCTd_L",  sn, "_k", ki)]] <- merge_draws(accD[[ki]])
      out[[paste0("CCTa_L",  sn, "_k", ki)]] <- merge_draws(accA[[ki]])
      out[[paste0("CCTad_L", sn, "_k", ki)]] <- merge_draws(accAD[[ki]])
    }
  }
  out
}

#' Ardia and Sessinou (2025) Subseries-Based Cauchy Combination Test (SCT) for Spanning
#'
#' Computes robust p-values for linear spanning restrictions using a residual-based
#' subseries procedure with Cauchy Combination Test (CCT) aggregation. Supports high-dimensional
#' inference for \eqn{H_0^\delta} (variance/spread slopes), \eqn{H_0^\alpha} (intercepts),
#' and the joint null \eqn{H_0^{\alpha,\delta}}.
#'
#' @param bench Numeric matrix of benchmark returns, dimension \eqn{T \times K}.
#' @param test  Numeric matrix of test-asset returns, dimension \eqn{T \times N}.
#' @param control Optional list passed to internal computation:
#' \describe{
#'   \item{\code{ks}}{Numeric vector of subseries exponents (block size approximately \code{floor(T^k)}); default \code{c(1/3)}.}
#'   \item{\code{L}}{Numeric vector of perturbation scales for randomized projections; default \code{c(0, 2)}.}
#'   \item{\code{B}}{Number of independent perturbation draws to merge when \code{L > 0}; default \code{1}. See \sQuote{Choosing B}.}
#'   \item{\code{seed}}{Seed of the first perturbation draw; default \code{123}. Draw \eqn{b} uses \code{seed + b - 1}.}
#' }
#'
#' @return A named list of global (combined) p-values. Names encode hypothesis and settings:
#' \itemize{
#'   \item \code{CCTd_L{L}_k{i}} — variance (slope) spanning, \eqn{\delta = 0};
#'   \item \code{CCTa_L{L}_k{i}} — alpha spanning, \eqn{\alpha = 0};
#'   \item \code{CCTad_L{L}_k{i}} — joint mean–variance spanning, \eqn{\alpha = 0,\ \delta = 0}.
#' }
#'
#' @details
#' For each \code{k} in \code{ks}, data are partitioned into overlapping subseries of length \code{floor(T^k)}.
#' Residual perturbations controlled by \code{L} generate test statistics robust to serial and
#' cross-sectional dependence and conditional heteroskedasticity. Resulting sub-p-values are aggregated
#' by the Cauchy Combination Test (CCT), which remains valid under dependence and retains power in high dimensions.
#'
#' @section Choosing B:
#' When \code{L > 0} the score is multiplied by a random weight vector. That
#' vector is shared by every test asset, so a single draw is a single realisation
#' of the test and its effect does \emph{not} average out across the
#' cross-section. The test is valid for any fixed draw --- its size is correct ---
#' but the p-value it returns on one data set can vary substantially from draw to
#' draw, especially when \code{T} is short (few subseries) or \code{N} is large.
#' \code{B > 1} draws \code{B} independent weight vectors and merges the
#' resulting p-values per asset with the Cauchy rule, which is valid under
#' arbitrary dependence and is not carried by a single extreme draw. The default
#' \code{B = 1} reproduces the historical single-draw behaviour exactly; for
#' applications, and for any result that will be reported, \code{B = 100} or more
#' is recommended. Cost is modest: the expensive residual construction is done
#' once and only the subseries statistic is recomputed per draw.
#'
#' Draw \eqn{b} is generated from \code{seed + b - 1}. Consecutive seeds are used
#' rather than an explicit substream generator so that \code{B = 1} reproduces
#' the historical single-draw result exactly. R scrambles the seed when
#' initialising the Mersenne-Twister, so the resulting weight vectors behave as
#' independent draws: across 100 consecutive seeds at \eqn{T = 250} the mean
#' absolute pairwise correlation is 0.051 and the maximum 0.23, against 0.050 and
#' about 0.22 expected under exact independence. Pass an explicit \code{seed} to
#' obtain a different, equally valid family of draws.
#'
#' @references
#' \insertRef{ArdiaSessinou2025}{spantest} \cr
#'
#' @examples
#' set.seed(123)
#' bench <- matrix(rnorm(300), 100, 3)
#' test  <- matrix(rnorm(200), 100, 2)
#' out <- span_as(bench, test)
#' out$CCTa_L0_k1; out$CCTd_L0_k1; out$CCTad_L0_k1
#'
#' @family Alpha Spanning Tests
#' @family Variance Spanning Tests
#' @family Joint Mean-Variance Spanning Tests
#'
#' @importFrom stats na.omit
#' @export
span_as <- function(bench, test, control = list()) {

  # Set control parameters
  con <- list(ks = c(1/3), L = c(0, 2), B = 1L, seed = 123L)
  con[names(control)] <- control
  k_values <- con$ks
  l_values <- con$L
  # B and seed must be whole numbers: as.integer() truncates, so a fractional
  # B = 1.9 would silently run a single draw rather than the two the caller meant.
  stopifnot(length(con$B) == 1L, is.finite(con$B), con$B >= 1,
            isTRUE(all.equal(con$B, round(con$B))),
            length(con$seed) == 1L, is.finite(con$seed),
            isTRUE(all.equal(con$seed, round(con$seed))))

  # Generate explicit template names
  test_types <- c("CCTd", "CCTad", "CCTa")
  template_names <- paste0(
    rep(test_types, each = length(l_values) * length(k_values)),
    "_L", rep(l_values, each = length(k_values), times = length(test_types)),
    "_k", rep(seq_along(k_values), times = length(l_values) * length(test_types))
  )

  # Per-asset p-values for the whole cross-section, computed in batch.
  pv <- f_getpv_batch(bench, test, ks = k_values, L = l_values,
                      B = as.integer(con$B), seed = as.integer(con$seed))

  # Combine per-asset p-values across the cross-section via the Cauchy method.
  combined_results <- vapply(template_names,
                             function(nm) f_cauchypv(na.omit(pv[[nm]])),
                             numeric(1))

  return(as.list(combined_results))
}
