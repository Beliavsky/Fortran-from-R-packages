#' Out-of-sample forecast metrics for the 2-level nested model
#'
#' Out-of-sample analogue of \code{\link{evaluate_oos}} for the nested model.
#' The decisive difference from the single-level recursion is the attractor: in
#' the single-level model the market price reverts to a constant level with a
#' linear forcing by the exogenous production-price covariate \eqn{X}; in the
#' nested model the market price \eqn{\varphi} reverts to the \emph{latent}
#' production price \eqn{\Phi}, which is itself propagated forward through its own
#' Level-2 Ornstein-Uhlenbeck equation toward the \eqn{G'}-driven mean
#' \eqn{\mu_{s,t} = m_{0,s} + m_1 G'_t}. Both states are advanced jointly; the
#' market never sees the realized future \eqn{X}.
#'
#' @details
#' At each forecast origin \eqn{t_0} the recursion is seeded with the realized
#' standardized market price \eqn{\varphi_{t_0-1}} and the posterior-median latent
#' production price \eqn{\hat\Phi_{t_0-1}}, then for \eqn{h = 1,\dots} steps:
#' \deqn{\Phi_{t} = \Phi_{t-1} + \kappa^p_s(\mu_{s,t} - \Phi_{t-1}),\quad
#'   \mu_{s,t} = m_{0,s} + m_1 G'_t,}
#' \deqn{\varphi_{t} = \varphi_{t-1} + \kappa^m_s(-d_{t-1} + a_{3,s} d_{t-1}^3)
#'   + \gamma\, \mathrm{COM}_{t-1,s},\quad d_{t-1} = \varphi_{t-1} - \Phi_{t-1},}
#' with \eqn{\kappa^m_s = \texttt{kappa\_cap}\cdot
#' \mathrm{logit}^{-1}(\tilde\kappa_s + \beta_1 zTMG_t)} the bounded gravitation
#' speed (the same logit link used in fitting). The market step uses
#' \eqn{\Phi_{t-1}} (one-step lag), matching the Stan likelihood.
#'
#' For the Level-3 value model (\code{n_levels = 3}), the Level-2 attractor also
#' tracks the value index, \eqn{\mu_{s,t} = m_{0,s} + m_1 G'_t + m_v V_{s,t}}
#' (production prices gravitate around values, Capital III ch. IX), recovered by
#' supplying the standardized value anchor \code{V_z} and the coupling
#' \code{meds$m_v}. Without them the recursion is the plain Level-2 mean.
#'
#' Two caveats, identical in spirit to \code{\link{evaluate_oos}}:
#' \enumerate{
#'   \item \strong{Conditional forecast.} The recursion still consumes the
#'     realized future \eqn{G'_t} and \eqn{zTMG_t} (the macro drivers), so it is a
#'     forecast \emph{conditional} on those aggregate paths, not on the sectoral
#'     production prices (which are now latent and forecast endogenously).
#'   \item \strong{Plug-in medians.} States are propagated with posterior medians
#'     and ignore parameter uncertainty, the latent-\eqn{\Phi} innovation, the
#'     stochastic volatility and the Student-t tails. It is a point (mean-equation)
#'     forecast, not the full posterior predictive distribution. Fit with
#'     \code{fit_window = "train"} for genuinely held-out numbers.
#' }
#'
#' @param meds List of posterior medians with components \code{kappa_tilde}
#'   (length S, the per-sector pre-link market reversion), \code{a3_s} (length S),
#'   \code{beta1} (scalar), \code{gamma} (scalar), \code{kappa_p} (length S),
#'   \code{mu_const} (length S, the Level-2 mean intercept \eqn{m_{0,s}}) and
#'   \code{m1} (scalar, the \eqn{G'} slope). For the Level-3 value model it may
#'   also carry \code{m_v} (scalar, the value coupling); paired with \code{V_z}
#'   it adds the value term to the Level-2 attractor.
#' @param Yz Numeric matrix (T x S). Standardized market price \eqn{\varphi}.
#' @param Phi_med Numeric matrix (T x S). Posterior-median latent production-price
#'   path \eqn{\hat\Phi} (seeds the Level-2 recursion).
#' @param Gprime_z Numeric vector (length T). Standardized aggregate profit rate
#'   \eqn{G'} (training statistics), the Level-2 mean driver.
#' @param zTMG Numeric vector (length T). Standardized TMG series used in fitting.
#' @param T_train Integer. End of the training window.
#' @param COM_ts Numeric matrix (T x S). Composition of capital by sector.
#' @param K_ts Numeric matrix (T x S). Total capital (training weights for the
#'   COM standardization, matching the Stan transformed-data block).
#' @param kappa_cap Numeric > 0. Stability cap used in fitting.
#' @param com_in_mean Logical. Whether the COM term enters the mean. Default FALSE.
#' @param horizons Integer vector. Forecast horizons to evaluate.
#' @param V_z Numeric matrix (T x S) or NULL. Standardized value anchor (direct
#'   prices \eqn{c+v+p}) for the Level-3 model. When supplied together with a
#'   \code{meds$m_v} coupling, the Level-2 attractor gains the value term so
#'   \eqn{\mu_{s,t} = m_{0,s} + m_1 G'_t + m_v V_{s,t}} (faithful out-of-sample
#'   recursion for \code{n_levels = 3}). NULL (default) reproduces the Level-2
#'   \eqn{G'}-driven mean exactly.
#'
#' @return Named list with one element per horizon, each a list with \code{h},
#'   \code{RMSE}, \code{MAE} and \code{n_obs} (pooled sector-by-origin errors;
#'   \code{0} / \code{NA} when the horizon exceeds the test window). Errors are in
#'   standardized \eqn{\varphi} units, directly comparable with
#'   \code{\link{evaluate_oos}}.
#'
#' @seealso \code{\link{evaluate_oos}} (single-level), \code{\link{fit_ou_nested}}.
#' @export
evaluate_oos_nested <- function(meds, Yz, Phi_med, Gprime_z, zTMG, T_train,
                                COM_ts, K_ts, kappa_cap, com_in_mean = FALSE,
                                horizons = c(1, 4, 8), V_z = NULL) {
  Tn <- nrow(Yz); S <- ncol(Yz)

  # Level-3 value coupling (D-IMPL-10): when the fit carries the value coupling
  # m_v and the standardized value anchor V_z (direct prices c+v+p) is supplied,
  # the latent production price reverts to a mean that also tracks the value,
  # mu_{s,t} = m0_s + m1 G'_t + m_v V_{s,t} (Capital III ch. IX). Inert (the plain
  # Level-2 G'-driven mean) when either the coupling or the anchor is absent, so
  # the n_levels <= 2 forecast is bit-identical.
  has_value <- !is.null(V_z) && !is.null(meds$m_v) &&
    length(meds$m_v) == 1L && is.finite(meds$m_v)
  if (has_value) V_z <- as.matrix(V_z)

  # COM standardization on the training window, K-weighted (mirrors the Stan
  # transformed-data block so the COM term matches the fitted model).
  com_wmean <- numeric(S); com_wsd <- numeric(S)
  for (s in seq_len(S)) {
    denom <- sum(K_ts[seq_len(T_train), s], na.rm = TRUE)
    if (!is.finite(denom) || denom <= 0) denom <- 1
    w <- K_ts[seq_len(T_train), s] / denom
    com_wmean[s] <- sum(COM_ts[seq_len(T_train), s] * w, na.rm = TRUE)
    v <- sum(w * (COM_ts[seq_len(T_train), s] - com_wmean[s])^2, na.rm = TRUE)
    com_wsd[s] <- sqrt(max(v, 1e-16))
  }

  inv_logit <- function(x) 1 / (1 + exp(-x))

  res <- lapply(horizons, function(hh) {
    last_origin <- Tn - hh + 1
    if (last_origin < (T_train + 1)) {
      return(list(h = hh, RMSE = NA_real_, MAE = NA_real_, n_obs = 0L))
    }
    errs <- c()

    for (t0 in seq.int(T_train + 1, last_origin)) {
      if (t0 - 1 < 1) next
      phi_pred <- Yz[t0 - 1, ]
      Phi_pred <- Phi_med[t0 - 1, ]

      for (h in seq_len(hh)) {
        t_pred <- t0 - 1 + h
        if (t_pred > Tn) break
        ztmg <- zTMG[min(t_pred, Tn)]
        gp_t <- Gprime_z[min(t_pred, Tn)]

        phi_new <- numeric(S); Phi_new <- numeric(S)
        for (s in seq_len(S)) {
          # Level 1: market reverts to the (lagged) latent production price.
          dev     <- phi_pred[s] - Phi_pred[s]
          kappa_m <- kappa_cap * inv_logit(meds$kappa_tilde[s] + meds$beta1 * ztmg)
          com_denom <- if (com_wsd[s] > 0) com_wsd[s] else 1
          com_std <- (COM_ts[min(t_pred - 1, Tn), s] - com_wmean[s]) / com_denom
          com_term <- if (com_in_mean) meds$gamma * com_std else 0
          phi_new[s] <- phi_pred[s] +
            kappa_m * (-dev + meds$a3_s[s] * dev^3) + com_term

          # Level 2: latent production price reverts to the G'-driven mean
          # (plus the value term m_v * V when the Level-3 value model is active).
          mu_ts <- meds$mu_const[s] + meds$m1 * gp_t
          if (has_value) mu_ts <- mu_ts + meds$m_v * V_z[min(t_pred, Tn), s]
          Phi_new[s] <- Phi_pred[s] + meds$kappa_p[s] * (mu_ts - Phi_pred[s])
        }
        phi_pred <- phi_new; Phi_pred <- Phi_new
      }

      errs <- c(errs, Yz[t0 + hh - 1, ] - phi_pred)
    }

    list(
      h = hh,
      RMSE = if (length(errs)) sqrt(mean(errs^2, na.rm = TRUE)) else NA_real_,
      MAE  = if (length(errs)) mean(abs(errs), na.rm = TRUE) else NA_real_,
      n_obs = length(errs)
    )
  })

  names(res) <- paste0("h", horizons)
  res
}


#' Latent Level-2 mean trajectory mu_t = m0_s + m1 * G'_t
#'
#' Reconstructs the time-varying Level-2 mean toward which the latent production
#' price \eqn{\Phi} reverts, \eqn{\mu_{s,t} = m_{0,s} + m_1 G'_t}, with credible
#' bands obtained by evaluating the expression over the posterior draws of
#' \code{mu_const} (\eqn{m_{0,s}}) and \code{m1}. This is the structural hook to
#' the (not-yet-formalized) Level 3 of values: \eqn{\mu} is the slow attractor of
#' the production price, here driven by the aggregate profit rate \eqn{G'}.
#'
#' @param fit Fitted Stan object (\code{CmdStanMCMC} or \code{stanfit}) from a
#'   2-level fit.
#' @param Gprime_z Numeric vector (length T). Standardized aggregate profit rate
#'   used in fitting (training statistics).
#' @param sector_names Character vector (length S) or NULL. Column names for the
#'   returned matrices.
#' @param V_z Numeric matrix (T x S) or NULL. Standardized value anchor (direct
#'   prices \eqn{c+v+p}) for the Level-3 model (\code{n_levels = 3}). When supplied
#'   and the fit carries the value coupling \code{m_v}, the mean trajectory gains
#'   the value term \eqn{m_v V_{s,t}}, so \eqn{\mu_{s,t} = m_{0,s} + m_1 G'_t + m_v
#'   V_{s,t}} (production prices gravitate around values, Capital III ch. IX). NULL
#'   (default) reproduces the Level-2 \eqn{G'}-driven mean.
#'
#' @return List with \code{median}, \code{q2.5}, \code{q97.5} (each a T x S
#'   matrix) and \code{Gprime_z} (the driver), or \code{NULL} if the Level-2
#'   parameters are not present in the fit.
#'
#' @seealso \code{\link{fit_ou_nested}}.
#' @export
extract_mu_trajectory <- function(fit, Gprime_z, sector_names = NULL, V_z = NULL) {
  M0 <- .draws_matrix(fit, "mu_const")   # [draws x S]
  M1 <- .draws_matrix(fit, "m1")         # [draws x 1]
  if (is.null(M0) || is.null(M1)) return(NULL)
  Tn <- length(Gprime_z); S <- ncol(M0)
  m1_vec <- M1[, 1]
  # Level-3 value coupling: add m_v * V when both the draws and the anchor exist.
  Mv <- .draws_matrix(fit, "m_v")
  has_value <- !is.null(Mv) && !is.null(V_z) &&
    nrow(as.matrix(V_z)) == Tn && ncol(as.matrix(V_z)) == S
  if (has_value) { mv_vec <- Mv[, 1]; V_z <- as.matrix(V_z) }

  med <- matrix(NA_real_, Tn, S); lo <- med; hi <- med
  for (s in seq_len(S)) {
    # draws x T matrix of mu_{s,t} = m0_s[draw] + m1[draw] * G'_t (+ m_v[draw] * V_{s,t})
    mu_st <- outer(M0[, s], rep(1, Tn)) + outer(m1_vec, Gprime_z)
    if (has_value) mu_st <- mu_st + outer(mv_vec, V_z[, s])
    q <- apply(mu_st, 2, stats::quantile, probs = c(0.025, 0.5, 0.975),
               na.rm = TRUE)
    lo[, s] <- q[1, ]; med[, s] <- q[2, ]; hi[, s] <- q[3, ]
  }
  if (!is.null(sector_names) && length(sector_names) == S) {
    colnames(med) <- colnames(lo) <- colnames(hi) <- sector_names
  }
  list(median = med, q2.5 = lo, q97.5 = hi, Gprime_z = Gprime_z)
}


#' Median vector of a (possibly indexed) Stan parameter
#'
#' Backend-agnostic column-wise posterior median for a parameter whose draws form
#' a \code{[draws x k]} matrix (scalar -> length 1, vector -> length k).
#' @keywords internal
#' @noRd
.median_vec <- function(fit, p) {
  M <- .draws_matrix(fit, p)
  if (is.null(M)) return(NULL)
  apply(M, 2, stats::median, na.rm = TRUE)
}
