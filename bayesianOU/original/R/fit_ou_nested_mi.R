#' Fit the nested OU model under multiple imputation of the market price
#'
#' Couples the disaggregation posterior to the OU model by multiple imputation
#' (Rubin's rule). The disaggregation step produces \eqn{M} posterior draws of
#' the sector-level market price \eqn{\varphi} (each a T x S matrix); this driver
#' refits the OU model once per imputation and combines the analyses, propagating
#' the disaggregation uncertainty into the OU posterior. Two complementary
#' summaries are returned: Rubin's moment-pooled estimates (within + between
#' imputation variance) and the mixture posterior over the imputations (the
#' "complete-draws" pooling used for joint probabilities such as the time-scale
#' separation).
#'
#' Unified external multiple imputation (D-IMPL-13.5). Beyond the market price,
#' \code{X} (the production-price anchor \eqn{\Phi}), \code{V_value} (the value
#' anchor), \code{CAPITAL_TOTAL} (\eqn{K}) and \code{COM} may also vary by
#' imputation. Each of these four accepts either a fixed matrix (the legacy
#' behaviour: used verbatim for every imputation, bit-identical results) or a
#' per-imputation set -- a list of \eqn{D} matrices or a 3-D array
#' \code{[T, S, D]} -- paired one-to-one with \code{phi_draws}. This lets a single
#' external imputation generator propagate ALL the construction uncertainty
#' (disaggregation of the market, the imputed split of \eqn{v} and \eqn{p} that
#' feeds \eqn{K} and \eqn{V}, and the reconstructed \eqn{\Phi}) through the
#' convergent K-DET / Variant-A engine, pooling by Rubin's rule. The aggregates
#' \code{TMG} and \code{Gprime} stay fixed (they are observed, preserved by the
#' aggregate-preserving renormalization of the generator). The engine is untouched
#' (R-side only): \eqn{n \in \{1, 2\}} and the all-fixed path remain bit-identical.
#'
#' @param phi_draws Imputations of the market price. Either a list of \eqn{D}
#'   numeric matrices (each T x S) or a 3-D array \code{[T, S, D]}. These are
#'   posterior draws of \eqn{\varphi} from the disaggregation model.
#' @param X Production-price index \eqn{\Phi} (the Level-2 anchor). Either a fixed
#'   numeric matrix (T x S) -- used verbatim for every imputation, the legacy
#'   behaviour -- or a per-imputation set (a list of \eqn{D} matrices or a 3-D
#'   array \code{[T, S, D]}) paired one-to-one with \code{phi_draws}: draw \eqn{m}
#'   of \code{X} is used with draw \eqn{m} of \code{phi_draws}.
#' @param TMG Numeric vector (length T). Aggregate profit-rate series. Fixed across
#'   imputations (observed aggregate).
#' @param COM Composition of capital by sector. A fixed matrix (T x S) or a
#'   per-imputation list/array \code{[T, S, D]} paired with \code{phi_draws}
#'   (see \code{X}).
#' @param CAPITAL_TOTAL Total capital advanced \eqn{K} by sector. A fixed matrix
#'   (T x S) or a per-imputation list/array \code{[T, S, D]} paired with
#'   \code{phi_draws} (see \code{X}).
#' @param Gprime Numeric vector (length T) or NULL. Aggregate profit rate driving
#'   the Level-2 mean; required when \code{n_levels >= 2}. Fixed across imputations
#'   (observed aggregate, preserved by construction).
#' @param V_value Value anchor (direct prices \eqn{c+v+p} = DirectPrices_Index) for
#'   the Level-3 model; required when \code{n_levels = 3}, NULL otherwise. A fixed
#'   matrix (T x S) or a per-imputation list/array \code{[T, S, D]} paired with
#'   \code{phi_draws} (see \code{X}); passed to \code{\link{fit_ou_nested}}, which
#'   standardizes it on the training window.
#' @param n_levels Integer, 1, 2 or 3. Passed to \code{\link{fit_ou_nested}}.
#' @param M Integer. Number of imputations to use. Default 25. If fewer draws are
#'   supplied, all are used (with a warning); if more, \code{M} are taken on an
#'   evenly spaced grid.
#' @param pool_params Character vector or NULL. Parameters to Rubin-pool. NULL
#'   selects a sensible set for the chosen \code{n_levels}.
#' @param keep_fits Logical. Keep the per-imputation fit objects (memory heavy).
#'   Default FALSE (only compact per-imputation summaries are retained).
#' @param keep_kappa_mixture Logical. When TRUE (and \code{n_levels >= 2}), retain
#'   in \code{$kappa_mixture} a thinned sample of the mixture posterior of the
#'   per-sector reversion speeds \code{kappa_m_base} and \code{kappa_p} (stacked
#'   across imputations), for plotting the sectoral distribution of half-lives and
#'   computing low-kappa-trap probabilities at a chosen horizon. Default FALSE.
#' @param kappa_mixture_draws Integer. Target number of thinned mixture draws to
#'   retain when \code{keep_kappa_mixture = TRUE}. Default 4000.
#' @param seed Integer. Base seed; imputation \code{m} uses \code{seed + m}.
#' @param verbose Logical. Print progress. Default FALSE.
#' @param checkpoint_dir Character path or NULL. When NULL (default) the function
#'   behaves exactly as before (no files written). When a directory is given,
#'   each imputation's compact pooling contribution is persisted to
#'   \code{checkpoint_dir/imp_<j>.rds} as soon as it is computed, and a later
#'   call with the same directory reloads the completed imputations and skips
#'   their (expensive) Stan refit -- resuming a run interrupted by a crash or
#'   shutdown. Only the pooling contributions are stored, never the full fits, so
#'   it is incompatible with \code{keep_fits = TRUE}. A stored checkpoint is
#'   honoured only when its draw index still matches the current imputation grid,
#'   so changing \code{M} (or the number of available draws) invalidates stale
#'   checkpoints automatically. Because the post-loop Rubin pooling is untouched
#'   and each per-imputation fit is deterministic given \code{seed + m}, a
#'   resumed run yields an object identical to an uninterrupted one.
#' @param ... Passed to \code{\link{fit_ou_nested}} (e.g. \code{chains},
#'   \code{iter}, \code{warmup}, \code{theta_separation}, \code{kappa_cap}).
#'
#' @return An object of class \code{"ou_nested_mi"}: a list with
#'   \describe{
#'     \item{rubin}{Data frame of Rubin-pooled estimates per parameter
#'       (estimate, total SD, 95\% CI, Barnard-Rubin df, fraction of missing
#'       information).}
#'     \item{phi_latent_pooled}{(\code{n_levels >= 2}) Rubin-pooled latent
#'       \eqn{\Phi} path: \code{mean}, \code{sd}, \code{q2.5}, \code{q97.5} as
#'       T x S matrices.}
#'     \item{separation_pooled}{(\code{n_levels >= 2}) mixture-posterior evidence
#'       for \eqn{\kappa^m > \kappa^p}: \code{prob_sep_by_sector} (per-sector
#'       \eqn{P(\kappa^m > \kappa^p)}, the sectoral distribution to report -- at
#'       \code{n_levels = 3} the joint probability typically collapses because
#'       \eqn{\kappa^p} is re-identified once the value term absorbs the production
#'       mean, see D-IMPL-10.1, which is NOT "no separation"), \code{prob_sep_joint},
#'       and the per-sector reversion summaries \code{kappa_m_median} /
#'       \code{kappa_p_median} and half-lives \eqn{(\ln 2)/\kappa}
#'       (\code{halflife_*_median}, \code{halflife_*_q2.5}, \code{halflife_*_q97.5})
#'       from the mixture posterior. Since \eqn{\kappa} is positive by construction,
#'       \eqn{P(\kappa>0)=1} is uninformative; the half-life (and whether it exceeds
#'       the observed span -- the low-kappa trap) is the substantive per-sector
#'       quantity.}
#'     \item{kappa_mixture}{(only when \code{keep_kappa_mixture = TRUE}) a thinned
#'       mixture posterior of \code{kappa_m_base} and \code{kappa_p} (draws x S),
#'       stacked across imputations, for sectoral half-life distributions and
#'       low-kappa-trap probabilities at a chosen horizon.}
#'     \item{per_imputation}{List of compact per-imputation summaries or full
#'       fits if \code{keep_fits = TRUE}. Each summary carries the chain-aware
#'       convergence over the pooled parameters (\code{rhat_max_pooled},
#'       \code{ess_bulk_min_pooled}, \code{ess_tail_min_pooled}) -- the binding
#'       diagnostic for Rubin's pooling -- plus the global \code{rhat_max},
#'       \code{rhat_share}, \code{ess_bulk_min} (over all parameters, latent
#'       states included; reference only) and \code{divergences}.}
#'     \item{config}{The imputation/pooling configuration.}
#'   }
#'
#' @references Rubin DB (1987) Multiple Imputation for Nonresponse in Surveys.
#'   Barnard J, Rubin DB (1999) Small-sample degrees of freedom with multiple
#'   imputation. Biometrika 86(4):948-955.
#'
#' @seealso \code{\link{fit_ou_nested}}.
#' @export
fit_ou_nested_mi <- function(phi_draws, X, TMG, COM, CAPITAL_TOTAL,
                             Gprime = NULL, V_value = NULL, n_levels = 2L, M = 25,
                             pool_params = NULL, keep_fits = FALSE,
                             keep_kappa_mixture = FALSE, kappa_mixture_draws = 4000L,
                             seed = 1234, verbose = FALSE, checkpoint_dir = NULL,
                             ...) {
  n_levels <- as.integer(n_levels)
  if (!n_levels %in% c(1L, 2L, 3L)) stop("`n_levels` must be 1, 2 or 3.", call. = FALSE)
  if (!is.null(checkpoint_dir) && isTRUE(keep_fits)) {
    stop("`checkpoint_dir` is incompatible with `keep_fits = TRUE`: the full ",
         "per-imputation fits are not persisted (only the compact pooling ",
         "contributions are), so a resumed run could not restore them. Use ",
         "`keep_fits = FALSE` (the default) with `checkpoint_dir`.", call. = FALSE)
  }
  if (n_levels == 3L && is.null(V_value)) {
    stop("`V_value` (indice de valor / precio directo c+v+p = DirectPrices_Index, ",
         "matriz T x S) es obligatorio para n_levels = 3 (Nivel 3, valores). ",
         "Es fijo entre imputaciones, como `X`.", call. = FALSE)
  }

  phi_list <- .as_phi_list(phi_draws)
  D <- length(phi_list)
  if (D < 1L) stop("`phi_draws` is empty.", call. = FALSE)

  # ---- Unified external MI (D-IMPL-13.5): the four model inputs that the
  #      imputation generator may reconstruct per draw. Each is either a fixed
  #      matrix (legacy: used verbatim every imputation, bit-identical) or a
  #      per-imputation set (list of D matrices / [T, S, D] array) paired with
  #      phi_draws by the SAME draw index. Validated against D up-front so a
  #      mismatched length fails before any Stan fit. TMG and Gprime stay fixed
  #      (observed aggregates). ----
  X_imp <- .as_imp_input(X, "X", D)
  V_imp <- .as_imp_input(V_value, "V_value", D)
  K_imp <- .as_imp_input(CAPITAL_TOTAL, "CAPITAL_TOTAL", D)
  C_imp <- .as_imp_input(COM, "COM", D)

  M <- as.integer(M)
  if (M < 2L) stop("`M` must be >= 2 for Rubin pooling.", call. = FALSE)
  if (D < M) {
    warning(sprintf("Only %d imputations supplied (< M = %d); using all %d.",
                    D, M, D), call. = FALSE)
    idx <- seq_len(D); M <- D
  } else if (D > M) {
    idx <- unique(round(seq(1, D, length.out = M)))
    M <- length(idx)
  } else {
    idx <- seq_len(D)
  }

  if (is.null(pool_params)) {
    pool_params <- if (n_levels >= 2L) {
      base2 <- c("kappa_m_base", "kappa_p", "sigma_p", "mu_const", "m1",
                 "sigma_phi_meas", "a3_s", "beta1", "gamma", "nu")
      # Level 3 (D-IMPL-10): add the value coupling m_v to the pooled set.
      if (n_levels == 3L) c(base2, "m_v") else base2
    } else {
      c("kappa_s", "a3_s", "theta_s", "beta0_s", "beta1", "gamma", "nu")
    }
  }

  # Grid dimensions and output column names from a representative anchor (the
  # first used draw): the [T, S] grid is invariant across imputations.
  X_rep <- .pick_imp(X_imp, idx[1L])
  Tn <- nrow(X_rep); S <- ncol(X_rep)
  Qm <- list(); Um <- list()                 # per-param: M x P matrices of means/vars
  sep_km <- list(); sep_kp <- list()          # stacked draws for the mixture
  Phi_Q <- if (n_levels >= 2L) matrix(NA_real_, M, Tn * S) else NULL
  Phi_U <- if (n_levels >= 2L) matrix(NA_real_, M, Tn * S) else NULL
  per_imp <- vector("list", M)
  fits    <- if (keep_fits) vector("list", M) else NULL

  # ---- Per-imputation checkpointing / resume (crash-resilience) ----
  # When `checkpoint_dir` is set, each imputation's pooling contribution is
  # persisted as `imp_<j>.rds` after it is computed; on a later call with the
  # same directory, completed imputations are reloaded and their Stan fit is
  # skipped. The post-loop Rubin pooling is untouched and each per-imputation
  # fit is deterministic given `seed + m`, so a resumed run yields an object
  # identical to an uninterrupted one. `checkpoint_dir = NULL` (default) is a
  # no-op: the legacy path is bit-identical. A checkpoint is only honoured when
  # its stored draw index matches `idx[j]` (so it is invalidated automatically
  # if `M` / the imputation grid changes).
  ckpt_loaded <- NULL
  if (!is.null(checkpoint_dir)) {
    dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
    ckpt_loaded <- vector("list", M)
    for (j in seq_len(M)) {
      f <- file.path(checkpoint_dir, sprintf("imp_%d.rds", j))
      if (file.exists(f)) {
        cj <- tryCatch(readRDS(f), error = function(e) NULL)
        if (!is.null(cj) && identical(cj$m, idx[j])) ckpt_loaded[[j]] <- cj
      }
    }
    vmsg(sprintf("Checkpoint: %d/%d imputations restored from %s",
                 sum(!vapply(ckpt_loaded, is.null, logical(1))), M, checkpoint_dir),
         verbose)
  }

  for (j in seq_len(M)) {
    m <- idx[j]
    # Resume: reload a completed imputation's contribution and skip its refit.
    if (!is.null(ckpt_loaded) && !is.null(ckpt_loaded[[j]])) {
      cj <- ckpt_loaded[[j]]
      for (p in names(cj$qm)) {
        if (is.null(Qm[[p]])) {
          Qm[[p]] <- matrix(NA_real_, M, length(cj$qm[[p]]),
                            dimnames = list(NULL, names(cj$qm[[p]])))
          Um[[p]] <- Qm[[p]]
        }
        Qm[[p]][j, ] <- cj$qm[[p]]
        Um[[p]][j, ] <- cj$um[[p]]
      }
      if (n_levels >= 2L) {
        if (!is.null(cj$phi_q)) { Phi_Q[j, ] <- cj$phi_q; Phi_U[j, ] <- cj$phi_u }
        if (!is.null(cj$sep_km)) { sep_km[[j]] <- cj$sep_km; sep_kp[[j]] <- cj$sep_kp }
      }
      per_imp[[j]] <- cj$per_imp
      vmsg(sprintf("Imputation %d/%d (draw %d) restored from checkpoint", j, M, m),
           verbose)
      next
    }
    vmsg(sprintf("Imputation %d/%d (draw %d)", j, M, m), verbose)
    Ym <- phi_list[[m]]
    fitm <- fit_ou_nested(Y = Ym, X = .pick_imp(X_imp, m), TMG = TMG,
                          COM = .pick_imp(C_imp, m),
                          CAPITAL_TOTAL = .pick_imp(K_imp, m), Gprime = Gprime,
                          V_value = .pick_imp(V_imp, m), n_levels = n_levels,
                          seed = seed + m, verbose = FALSE, ...)
    sf <- fitm$factor_ou$stan_fit

    for (p in pool_params) {
      Mp <- .draws_matrix(sf, p)
      if (is.null(Mp)) next
      if (is.null(Qm[[p]])) {
        Qm[[p]] <- matrix(NA_real_, M, ncol(Mp), dimnames = list(NULL, colnames(Mp)))
        Um[[p]] <- Qm[[p]]
      }
      Qm[[p]][j, ] <- colMeans(Mp)
      Um[[p]][j, ] <- apply(Mp, 2, stats::var)
    }

    if (n_levels >= 2L) {
      Mphi <- .draws_matrix(sf, "Phi")
      if (!is.null(Mphi)) {
        cn <- colnames(Mphi)
        ord <- match(.phi_colnames(Tn, S), cn)
        Phi_Q[j, ] <- colMeans(Mphi)[ord]
        Phi_U[j, ] <- apply(Mphi, 2, stats::var)[ord]
      }
      km <- .draws_matrix(sf, "kappa_m_base"); kp <- .draws_matrix(sf, "kappa_p")
      if (!is.null(km) && !is.null(kp)) { sep_km[[j]] <- km; sep_kp[[j]] <- kp }
    }

    # Chain-aware convergence over the POOLED parameters (those that govern the
    # reliability of Rubin's pooling): max R-hat and min ESS-bulk/tail. This is
    # the binding per-imputation diagnostic. The global rhat_max / ess_bulk_min
    # (over ALL parameters, latent states included) are kept only for reference:
    # a stray latent cell can drag the global extreme without bearing on the
    # pooled estimates.
    pr <- pe_b <- pe_t <- numeric(0)
    if (requireNamespace("posterior", quietly = TRUE)) {
      for (p in pool_params) {
        da <- tryCatch(sf$draws(p), error = function(e) NULL)
        if (is.null(da)) next
        ss <- suppressWarnings(
          posterior::summarise_draws(da, "rhat", "ess_bulk", "ess_tail"))
        pr <- c(pr, ss$rhat); pe_b <- c(pe_b, ss$ess_bulk); pe_t <- c(pe_t, ss$ess_tail)
      }
    }
    fmin <- function(x) if (length(x)) min(x, na.rm = TRUE) else NA_real_
    fmax <- function(x) if (length(x)) max(x, na.rm = TRUE) else NA_real_
    per_imp[[j]] <- list(
      imputation = m,
      rhat_max = fitm$diagnostics$rhat_max,
      rhat_share = fitm$diagnostics$rhat_share,
      rhat_max_pooled = fmax(pr),
      ess_bulk_min_pooled = fmin(pe_b),
      ess_tail_min_pooled = fmin(pe_t),
      ess_bulk_min = suppressWarnings(min(fitm$diagnostics$ess, na.rm = TRUE)),
      divergences = fitm$diagnostics$divergences
    )
    if (keep_fits) fits[[j]] <- fitm

    # Persist this imputation's pooling contribution (compact; never the full
    # fit). Written last, so a half-finished imputation leaves no checkpoint and
    # is simply refit on resume (deterministic given `seed + m`).
    if (!is.null(checkpoint_dir)) {
      qm_j <- um_j <- list()
      for (p in pool_params) {
        if (!is.null(Qm[[p]])) { qm_j[[p]] <- Qm[[p]][j, ]; um_j[[p]] <- Um[[p]][j, ] }
      }
      cj <- list(
        j = j, m = m, qm = qm_j, um = um_j,
        phi_q   = if (n_levels >= 2L) Phi_Q[j, ] else NULL,
        phi_u   = if (n_levels >= 2L) Phi_U[j, ] else NULL,
        sep_km  = if (n_levels >= 2L && length(sep_km) >= j) sep_km[[j]] else NULL,
        sep_kp  = if (n_levels >= 2L && length(sep_kp) >= j) sep_kp[[j]] else NULL,
        per_imp = per_imp[[j]])
      # Escritura atómica: tmp + rename en el mismo FS. Un kill a mitad de la
      # serialización deja a lo sumo un .tmp huérfano, nunca un imp_<j>.rds
      # truncado que el resume pudiera leer a medias.
      ckpt_f   <- file.path(checkpoint_dir, sprintf("imp_%d.rds", j))
      ckpt_tmp <- paste0(ckpt_f, ".tmp")
      saveRDS(cj, ckpt_tmp)
      file.rename(ckpt_tmp, ckpt_f)
    }
  }

  # ---- Rubin pooling of the scalar/vector parameters ----
  rubin_rows <- lapply(names(Qm), function(p) {
    r <- .rubin_combine(Qm[[p]], Um[[p]], M)
    data.frame(parameter = colnames(Qm[[p]]) %||% p,
               estimate = r$estimate, total_sd = r$total_sd,
               q2.5 = r$lo, q97.5 = r$hi, df = r$df, fmi = r$fmi,
               row.names = NULL, stringsAsFactors = FALSE)
  })
  rubin <- do.call(rbind, rubin_rows)

  out <- list(rubin = rubin,
              per_imputation = if (keep_fits) fits else per_imp,
              config = list(M = M, n_levels = n_levels, pool_params = pool_params,
                            n_available = D, keep_fits = keep_fits))

  if (n_levels >= 2L) {
    rphi <- .rubin_combine(Phi_Q, Phi_U, M)
    rs   <- function(v) matrix(v, Tn, S, dimnames = list(NULL, colnames(X_rep)))
    out$phi_latent_pooled <- list(mean = rs(rphi$estimate), sd = rs(rphi$total_sd),
                                  q2.5 = rs(rphi$lo), q97.5 = rs(rphi$hi))
    if (length(sep_km) > 0) {
      KM <- do.call(rbind, sep_km); KP <- do.call(rbind, sep_kp)  # mixture posterior
      cn_sec <- colnames(X_rep)
      # kappa is strictly positive by construction (kappa = kappa_cap *
      # inv_logit(.) in (0, kappa_cap)), so P(kappa > 0) = 1 is uninformative.
      # The substantive per-sector quantity (spec S2.3 / S2.6) is the half-life
      # (ln 2)/kappa and whether it exceeds the observed span (the "low-kappa
      # trap": reversion too slow to be told apart from a unit root). The mixture
      # median / 95% interval is reported rather than a moment-based delta on
      # kappa because kappa_p piles near the lower boundary exactly in the slow
      # regime, where a Gaussian/t moment approximation is least faithful. Note
      # the half-life is decreasing in kappa, so its 2.5% quantile maps to the
      # 97.5% quantile of kappa and vice versa.
      hl  <- function(K, s, q) log(2) / stats::quantile(K[, s], probs = q,
                                                        names = FALSE, na.rm = TRUE)
      med <- function(K, s) stats::median(K[, s], na.rm = TRUE)
      out$separation_pooled <- list(
        sectors            = cn_sec,
        prob_sep_by_sector = vapply(seq_len(S), function(s) mean(KM[, s] > KP[, s]), numeric(1)),
        prob_sep_joint     = mean(apply(KM > KP, 1, all)),
        kappa_m_median     = vapply(seq_len(S), function(s) med(KM, s),    numeric(1)),
        kappa_p_median     = vapply(seq_len(S), function(s) med(KP, s),    numeric(1)),
        halflife_m_median  = vapply(seq_len(S), function(s) hl(KM, s, 0.5),   numeric(1)),
        halflife_m_q2.5    = vapply(seq_len(S), function(s) hl(KM, s, 0.975), numeric(1)),
        halflife_m_q97.5   = vapply(seq_len(S), function(s) hl(KM, s, 0.025), numeric(1)),
        halflife_p_median  = vapply(seq_len(S), function(s) hl(KP, s, 0.5),   numeric(1)),
        halflife_p_q2.5    = vapply(seq_len(S), function(s) hl(KP, s, 0.975), numeric(1)),
        halflife_p_q97.5   = vapply(seq_len(S), function(s) hl(KP, s, 0.025), numeric(1)))
      if (isTRUE(keep_kappa_mixture)) {
        nr  <- nrow(KM)
        sub <- if (nr > kappa_mixture_draws)
          round(seq(1, nr, length.out = kappa_mixture_draws)) else seq_len(nr)
        colnames(KM) <- cn_sec; colnames(KP) <- cn_sec
        out$kappa_mixture <- list(kappa_m_base = KM[sub, , drop = FALSE],
                                  kappa_p      = KP[sub, , drop = FALSE])
      }
    }
  }

  class(out) <- c("ou_nested_mi", "list")
  out
}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' Normalize a possibly-per-imputation model input (D-IMPL-13.5)
#'
#' Each of X / V_value / CAPITAL_TOTAL / COM may be supplied either as a fixed
#' matrix (used verbatim for every imputation -- the legacy, bit-identical path)
#' or as a per-imputation set paired with \code{phi_draws}: a list of \code{D}
#' matrices, or a 3-D array \code{[T, S, D]}. A data frame counts as a single
#' fixed matrix (it is a list, hence the explicit guard). Returns \code{NULL}
#' for a \code{NULL} input (e.g. \code{V_value} below Level 3). The per-imputation
#' length is validated against \code{D} so a mismatch fails before any Stan fit.
#' @param x The input (NULL, matrix/data.frame, list of matrices, or [T,S,D] array).
#' @param name Argument name, for error messages.
#' @param D Number of available imputation draws (length of \code{phi_draws}).
#' @return NULL, or a list with \code{$varies} (logical) and either \code{$fixed}
#'   (the matrix, when not varying) or \code{$list} (length-\code{D} list of
#'   matrices, when varying).
#' @keywords internal
#' @noRd
.as_imp_input <- function(x, name, D) {
  if (is.null(x)) return(NULL)
  # 3-D array [T, S, D] -> per-imputation
  if (is.array(x) && length(dim(x)) == 3L) {
    d <- dim(x)
    if (d[3L] != D) {
      stop(sprintf("`%s` is a [T, S, %d] array but there are %d imputation draws; ",
                   name, d[3L], D),
           "a per-imputation input must have one slice per draw of `phi_draws`.",
           call. = FALSE)
    }
    lst <- lapply(seq_len(d[3L]), function(k) x[, , k, drop = TRUE])
    return(list(varies = TRUE, list = lst, fixed = NULL))
  }
  # List of matrices (but NOT a data frame, which is also a list) -> per-imputation
  if (is.list(x) && !is.data.frame(x)) {
    if (length(x) != D) {
      stop(sprintf("`%s` is a list of %d matrices but there are %d imputation ",
                   name, length(x), D),
           "draws; a per-imputation input must have one matrix per draw of ",
           "`phi_draws`.", call. = FALSE)
    }
    lst <- lapply(x, as.matrix)
    return(list(varies = TRUE, list = lst, fixed = NULL))
  }
  # Single matrix / data frame -> fixed across imputations (legacy path)
  list(varies = FALSE, list = NULL, fixed = as.matrix(x))
}

#' Pick the per-imputation matrix for draw \code{m} (NULL passes through)
#' @keywords internal
#' @noRd
.pick_imp <- function(imp, m) {
  if (is.null(imp)) return(NULL)
  if (imp$varies) imp$list[[m]] else imp$fixed
}

#' Coerce phi_draws (list of matrices or 3-D array) to a list of matrices
#' @keywords internal
#' @noRd
.as_phi_list <- function(phi_draws) {
  if (is.list(phi_draws)) {
    return(lapply(phi_draws, as.matrix))
  }
  if (is.array(phi_draws) && length(dim(phi_draws)) == 3L) {
    d <- dim(phi_draws)
    return(lapply(seq_len(d[3]), function(k) phi_draws[, , k, drop = TRUE]))
  }
  stop("`phi_draws` must be a list of T x S matrices or a [T, S, D] array.",
       call. = FALSE)
}

#' Canonical column order "Phi[t,s]" (t fastest, matching Stan/cmdstanr df order)
#' @keywords internal
#' @noRd
.phi_colnames <- function(Tn, S) {
  as.vector(vapply(seq_len(S), function(s) sprintf("Phi[%d,%d]", seq_len(Tn), s),
                   character(Tn)))
}

#' Rubin's rules: combine per-imputation means (Qm) and variances (Um)
#'
#' Qm, Um are M x P matrices. Returns pooled estimate, total SD, 95\% CI under the
#' Barnard-Rubin t reference distribution, df and the fraction of missing
#' information. When the between-imputation variance is ~0 the df tends to
#' infinity (the t collapses to a normal).
#' @keywords internal
#' @noRd
.rubin_combine <- function(Qm, Um, M) {
  eps  <- .Machine$double.eps
  Qbar <- colMeans(Qm, na.rm = TRUE)
  Ubar <- colMeans(Um, na.rm = TRUE)
  B    <- apply(Qm, 2, stats::var, na.rm = TRUE)        # between-imputation (denom M-1)
  Tvar <- Ubar + (1 + 1 / M) * B
  r    <- (1 + 1 / M) * B / pmax(Ubar, eps)             # relative variance increase
  df   <- (M - 1) * (1 + 1 / pmax(r, eps))^2
  df[!is.finite(df)] <- Inf
  lambda <- ((1 + 1 / M) * B) / pmax(Tvar, eps)         # proportion of variance from MI
  fmi  <- (r + 2 / (df + 3)) / (r + 1)
  fmi[!is.finite(fmi)] <- lambda[!is.finite(fmi)]
  tcrit <- ifelse(is.finite(df), stats::qt(0.975, df), stats::qnorm(0.975))
  list(estimate = Qbar, total_sd = sqrt(Tvar), within_var = Ubar, between_var = B,
       df = df, fmi = fmi,
       lo = Qbar - tcrit * sqrt(Tvar), hi = Qbar + tcrit * sqrt(Tvar))
}
