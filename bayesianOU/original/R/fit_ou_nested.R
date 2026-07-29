#' Fit the unified nonlinear OU model: single-level or 2-level nested
#'
#' Single engine behind the package. \code{n_levels = 1} reproduces the legacy
#' single-level model exactly (market price reverts to a constant level with a
#' linear forcing by the production price and a TMG interaction);
#' \code{n_levels = 2} fits the nested model in which the market price
#' \eqn{\varphi} reverts to a \emph{latent} production price \eqn{\Phi}
#' (Level 1), and \eqn{\Phi} has its own Ornstein-Uhlenbeck dynamics reverting
#' to a \eqn{G'}-driven mean (Level 2). The constructed production-price index
#' enters as a noisy measurement (anchor) of the latent \eqn{\Phi}.
#'
#' \code{fit_ou_nonlinear_tmg} is a thin backward-compatible wrapper around this
#' engine with \code{n_levels = 1}; this function is the single source of truth.
#'
#' @param Y Numeric matrix (T x S). Market prices/values \eqn{\varphi} by sector.
#' @param X Numeric matrix (T x S). Production prices \eqn{\Phi}. In single-level
#'   mode this is a linear forcing covariate; in 2-level mode it is the
#'   constructed production-price index that anchors the latent \eqn{\Phi}.
#' @param TMG Numeric vector (length T). Aggregate profit-rate series used as the
#'   TMG interaction (modulates the gravitation speed).
#' @param COM Numeric matrix (T x S). Composition of capital by sector.
#' @param CAPITAL_TOTAL Numeric matrix (T x S). Total capital by sector.
#' @param n_levels Integer, 1, 2 or 3. Model depth. \code{1} single-level
#'   (legacy); \code{2} nested (market \eqn{\varphi} reverts to latent production
#'   \eqn{\Phi}, which reverts to a \eqn{G'}-driven mean); \code{3} adds the value
#'   level (D-IMPL-10): the latent production price reverts toward a mean that also
#'   tracks the value anchor \eqn{V} (direct prices \eqn{c+v+p}), so prices of
#'   production gravitate around values (Marx, Capital III ch. IX). \code{n_levels
#'   = 3} requires \code{V_value} and \code{Gprime}; it inherits all of the 2-level
#'   machinery and adds the value coupling \eqn{m_v}. Default 1.
#' @param Gprime Numeric vector (length T) or NULL. Aggregate profit rate
#'   \eqn{G'_t = \sum EBO_t / \sum K_t} that drives the Level-2 mean. Required
#'   when \code{n_levels = 2}; ignored otherwise. In \code{k_uncertainty =
#'   "recon"} the same (raw) series also enters the price identity
#'   \eqn{\Phi = k + K G'}.
#' @param k_cost Numeric matrix (T x S) or NULL. Cost price \eqn{k = c + v} by
#'   sector (raw units). Required only when \code{k_uncertainty = "recon"}, where
#'   the latent production-price anchor is reconstructed as \eqn{\Phi = k + K G'}
#'   with the total capital advanced \eqn{K} (taken from \code{CAPITAL_TOTAL})
#'   treated as lognormally uncertain.
#' @param V_value Numeric matrix (T x S) or NULL. Value anchor (direct prices
#'   \eqn{V = c + v + p =} \code{DirectPrices_Index}, constructed directly from
#'   \code{k_cost + EBO}, indexed). Required only when \code{n_levels = 3}, where
#'   the latent production price reverts toward a mean that tracks \eqn{V} with the
#'   estimated coupling \eqn{m_v} (Capital III ch. IX). \eqn{V} is standardized on
#'   the training window and enters as a \emph{datum} (direct empirical
#'   construction), never solved simultaneously (no Leontief inverse).
#' @param level_spec List or NULL. Per-level richness toggles, each a list
#'   \code{level1}/\code{level2} of four logical switches: \code{cubic} (cubic
#'   drift), \code{sv} (stochastic volatility), \code{student_t} (Student-t
#'   tails) and \code{hierarchy} (cross-sector partial pooling). \code{NULL}
#'   selects the canonical specification (Level 1 full, Level 2 lean). The
#'   switches act only in the 2-level mode; in single-level mode Level 1 is always
#'   full (pass \code{NULL}). Use \code{\link{ou_level_spec}} for the named
#'   experiment configurations (\code{"canonical"}, \code{"both_full"},
#'   \code{"both_lean"}, \code{"n1_lean"}).
#' @param theta_separation Character. Time-scale separation \eqn{\kappa^m >
#'   \kappa^p} between market and production reversion: \code{"soft"} (a prior
#'   that favours it, falsifiable) or \code{"hard"} (imposed by construction).
#'   Only used when \code{n_levels = 2}. Default \code{"soft"}.
#' @param k_uncertainty Character. Treatment of the methodological uncertainty
#'   about the capital \eqn{K} (eq. 11), two contrastable modes. \code{"meas"}
#'   feeds the constructed production-price index as a noisy measurement of the
#'   latent \eqn{\Phi} with scale \eqn{\sigma_\Phi} (first-order measurement
#'   error). \code{"recon"} reconstructs the anchor as \eqn{\Phi = k + K G'} with
#'   \eqn{K} lognormally uncertain (prior \eqn{f_Y} around the point estimate
#'   \code{CAPITAL_TOTAL}, scale \code{priors$sigma_K_recon}); this propagates the
#'   uncertainty about \eqn{K} into the latent production price. The reconstruction
#'   is standardized so that the \eqn{\sigma_K \to 0} limit reproduces the
#'   \code{"meas"} anchor exactly (recon nests meas). \code{"recon"} requires
#'   \code{k_cost} and \code{n_levels = 2}. Default \code{"meas"}.
#' @param sigma_phi_meas_fixed Numeric > 0 or NULL. Latency dial for the anchor
#'   measurement SD \eqn{\sigma_\Phi} (2-level only; decision D-IMPL-9.4). When
#'   \code{NULL} (default) the SD is \emph{estimated} as a parameter with the
#'   half-normal prior \code{priors$sigma_phi_meas_sd} (the K-stochastic / recon
#'   regime, looser prior). When a positive number it is held \emph{fixed} to that
#'   value as a datum and no parameter is sampled: this is the
#'   \eqn{\sigma_\Phi \to 0} limit done well for the K-deterministic regime
#'   (\eqn{\Phi \approx} the constructed index), and it removes the boundary funnel
#'   (G6) the estimated SD fell into in the deterministic limit, where the
#'   anchor-minus-\eqn{\Phi} residuals vanish and the posterior of \eqn{\sigma_\Phi}
#'   piles up at 0 (Session-9 CP smoke: \code{rhat} 2.208 on \code{sigma_phi_meas}).
#'   A small value such as \code{0.05} reproduces the tight-latency K-deterministic
#'   anchor. Ignored when \code{n_levels = 1}.
#' @param kappa_cap Numeric > 0. Stability cap for the 2-level reversion speeds:
#'   both the market speed \eqn{\kappa^m} and the production speed \eqn{\kappa^p}
#'   are constrained to \eqn{(0, \texttt{kappa\_cap})} through a logit link (with
#'   the TMG modulation inside the link). The Euler-discretized OU is stable for
#'   \eqn{\kappa < 2} and monotone (no overshoot) for \eqn{\kappa < 1}; the
#'   default \code{2} admits the full stable region, oscillatory convergence
#'   included. Only used when \code{n_levels = 2}.
#' @param results_robust List. Previous results object to extend (or empty list).
#' @param priors List. Prior overrides (partial; missing entries use robust
#'   defaults). In addition to the single-level names (\code{sigma_delta},
#'   \code{beta1_mean}/\code{beta1_sd}, \code{nu_shape}/\code{nu_rate},
#'   \code{rho_mean}/\code{rho_sd}) the 2-level mode reads
#'   \code{sigma_phi_meas_sd} (half-normal scale of the anchor measurement SD,
#'   standardized units; default 0.5) and, in \code{"recon"} mode,
#'   \code{sigma_K_recon} (lognormal scale of the methodological uncertainty about
#'   \eqn{K}, in log units; default 0.10, i.e. roughly a 10\% multiplicative
#'   uncertainty on the total capital advanced).
#' @param com_in_mean Logical. Include COM effect in the mean. Default TRUE.
#' @param train_frac Numeric in (0,1). Training-window fraction. Default 0.70.
#' @param fit_window Character. \code{"train"} fits the likelihood on the
#'   training window only (honest forecasting); \code{"full"} uses all
#'   observations. Default \code{"train"}.
#' @param chains,iter,warmup,thin Integer MCMC controls. Defaults 6 / 12000 /
#'   6000 / 1.
#' @param cores Integer. Cores for parallel chains.
#' @param threads_per_chain Integer. Within-chain threads (capped so
#'   \code{parallel_chains * threads <= cores}).
#' @param hard_sum_zero Logical. Fix the TMG wedge at zero. Default TRUE.
#' @param orthogonalize_tmg Logical. Orthogonalize TMG w.r.t. the common factor.
#'   Default TRUE.
#' @param factor_from Character. Source for the common factor: \code{"X"} or
#'   \code{"Y"}. Default \code{"X"}.
#' @param use_train_loadings Logical. Compute factor loadings from training only.
#'   Default TRUE.
#' @param adapt_delta Numeric in (0,1). NUTS target acceptance. Default 0.97.
#' @param max_treedepth Integer. NUTS maximum tree depth. Default 12.
#' @param seed Integer. Random seed.
#' @param init Numeric, list or function. Initial values (passed to the backend).
#' @param moment_match Logical or NULL. Accepted for API parity; not applied to
#'   the array-based PSIS-LOO (see \code{fit_ou_nonlinear_tmg}).
#' @param verbose Logical. Print progress. Default FALSE.
#'
#' @return A list. For \code{n_levels = 1} the structure matches the legacy
#'   single-level fit, with the inner \code{factor_ou} carrying the historical S3
#'   class \code{"ou_nonlinear_tmg"}. For \code{n_levels = 2} the \emph{top-level}
#'   object carries the S3 class \code{"ou_nested_2level"} (so
#'   \code{\link[=print.ou_nested_2level]{print}},
#'   \code{\link[=summary.ou_nested_2level]{summary}} and
#'   \code{\link[=plot.ou_nested_2level]{plot}} dispatch on it) and adds
#'   \code{factor_ou$level2} (production reversion speed \code{kappa_p}, market
#'   speed at \eqn{zTMG=0} \code{kappa_m_base}, mean intercept \code{mu_const},
#'   slope \code{m1}, innovation scale \code{sigma_p}, measurement scale
#'   \code{sigma_phi_meas}, and the Level-2 richness quantities that are present
#'   only when their switch is on: the cubic restoring coefficient \code{a3_p}
#'   (\code{l2_cubic}), the Student-t degrees of freedom \code{nu_p}
#'   (\code{l2_studentt}) and the time-varying SV scale path \code{sigma_p_t} =
#'   \eqn{\exp(h_p/2)} (\code{l2_sv}); see decision D-IMPL-10);
#'   \code{phi_latent} (the latent \eqn{\Phi} trajectory
#'   with credible bands); \code{mu_path} (the \eqn{G'}-driven Level-2 mean
#'   trajectory \eqn{\mu_{s,t} = m_{0,s} + m_1 G'_t} with bands); a
#'   \code{separation} block with the posterior evidence for \eqn{\kappa^m >
#'   \kappa^p}; and \code{diagnostics$oos}, the 2-level out-of-sample metrics in
#'   which the market reverts to the latent \eqn{\Phi} propagated forward (see
#'   \code{\link{evaluate_oos_nested}}).
#'
#' @seealso \code{\link{fit_ou_nonlinear_tmg}} for the single-level wrapper and
#'   \code{\link{fit_ou_nested_mi}} for the multiple-imputation driver that
#'   couples disaggregation draws to the 2-level model via Rubin's rule.
#'
#' @examples
#' \donttest{
#' T_obs <- 24; S <- 2
#' Y <- matrix(rnorm(T_obs * S), T_obs, S)
#' X <- matrix(rnorm(T_obs * S), T_obs, S)
#' TMG <- rnorm(T_obs)
#' COM <- matrix(runif(T_obs * S), T_obs, S)
#' K <- matrix(runif(T_obs * S, 100, 1000), T_obs, S)
#' Gp <- rnorm(T_obs)
#' if (requireNamespace("cmdstanr", quietly = TRUE)) {
#'   try(fit_ou_nested(Y, X, TMG, COM, K, n_levels = 2, Gprime = Gp,
#'                     chains = 1, iter = 100, warmup = 50), silent = TRUE)
#' }
#' }
#'
#' @export
fit_ou_nested <- function(
    Y, X,
    TMG,
    COM,
    CAPITAL_TOTAL,
    n_levels = 1L,
    Gprime = NULL,
    k_cost = NULL,
    V_value = NULL,
    level_spec = NULL,
    theta_separation = c("soft", "hard"),
    k_uncertainty = c("meas", "recon"),
    sigma_phi_meas_fixed = NULL,
    kappa_cap = 2,
    results_robust = list(),
    priors = list(),
    com_in_mean = TRUE,
    train_frac = 0.70,
    fit_window = c("train", "full"),
    chains = 6,
    iter = 12000,
    warmup = 6000,
    thin = 1,
    cores = max(1, parallel::detectCores() - 1),
    threads_per_chain = 2,
    hard_sum_zero = TRUE,
    orthogonalize_tmg = TRUE,
    factor_from = c("X", "Y"),
    use_train_loadings = TRUE,
    adapt_delta = 0.97,
    max_treedepth = 12,
    seed = 1234,
    init = NULL,
    moment_match = NULL,
    verbose = FALSE
) {

  factor_from      <- match.arg(factor_from)
  fit_window       <- match.arg(fit_window)
  theta_separation <- match.arg(theta_separation)
  k_uncertainty    <- match.arg(k_uncertainty)

  n_levels <- as.integer(n_levels)
  if (length(n_levels) != 1L || !n_levels %in% c(1L, 2L, 3L)) {
    stop("`n_levels` must be 1, 2 or 3.", call. = FALSE)
  }

  # ---- Dispatch validation for options not yet realised by the Stan model ----
  level_spec <- .resolve_level_spec(level_spec, n_levels)
  if (identical(k_uncertainty, "recon") && n_levels < 2L) {
    stop("k_uncertainty = \"recon\" (reconstruccion con K incierto, ec. 11) ",
         "solo aplica con n_levels >= 2.", call. = FALSE)
  }
  if (n_levels == 3L && is.null(V_value)) {
    stop("`V_value` (indice de valor / precio directo c+v+p = DirectPrices_Index, ",
         "matriz T x S) es obligatorio para n_levels = 3 (Nivel 3, valores).",
         call. = FALSE)
  }

  # ---- Merge priors over robust defaults (partial override allowed) ----
  prior_defaults <- list(
    sigma_delta       = 0.002,
    beta1_mean        = 0,    beta1_sd = 0.5,
    nu_shape          = 2,    nu_rate  = 0.1,
    rho_mean          = 0.7,  rho_sd   = 0.2,
    sigma_phi_meas_sd = 0.5,           # half-normal scale of the anchor SD (2-level)
    sigma_K_recon     = 0.10           # lognormal scale of the K-uncertainty (recon mode)
  )
  priors <- utils::modifyList(prior_defaults, priors %||% list())

  # ---- Shared input validation ----
  stopifnot(is.matrix(Y) || is.data.frame(Y))
  stopifnot(is.matrix(X) || is.data.frame(X))
  Y <- as.matrix(Y); X <- as.matrix(X)
  stopifnot(nrow(Y) == nrow(X), ncol(Y) == ncol(X))
  stopifnot(length(TMG) == nrow(Y))

  if (!all(is.finite(Y)))   stop("`Y` contains non-finite values (NA/NaN/Inf).", call. = FALSE)
  if (!all(is.finite(X)))   stop("`X` contains non-finite values (NA/NaN/Inf).", call. = FALSE)
  if (!all(is.finite(TMG))) stop("`TMG` contains non-finite values.", call. = FALSE)
  if (!is.numeric(train_frac) || train_frac <= 0 || train_frac >= 1) {
    stop("`train_frac` must be in (0, 1).", call. = FALSE)
  }
  if (warmup < 1L || iter <= warmup) {
    stop(sprintf("`iter` (%d) must be strictly greater than `warmup` (%d).",
                 as.integer(iter), as.integer(warmup)), call. = FALSE)
  }
  if (chains < 1L) stop("`chains` must be >= 1.", call. = FALSE)
  if (thin < 1L)   stop("`thin` must be >= 1.", call. = FALSE)
  if (!is.numeric(kappa_cap) || length(kappa_cap) != 1L || kappa_cap <= 0) {
    stop("`kappa_cap` must be a single positive number.", call. = FALSE)
  }
  # ---- Measurement-SD latency dial (D-IMPL-9.4): fixed datum vs estimated ----
  if (!is.null(sigma_phi_meas_fixed)) {
    smf <- as.numeric(sigma_phi_meas_fixed)
    if (length(smf) != 1L || !is.finite(smf) || smf <= 0) {
      stop("`sigma_phi_meas_fixed` must be a single positive number (the fixed ",
           "anchor measurement SD) or NULL.", call. = FALSE)
    }
    if (n_levels != 2L) {
      warning("`sigma_phi_meas_fixed` only acts when n_levels = 2 (the measurement ",
              "block); ignored for the single-level model.", call. = FALSE)
    }
  }

  stopifnot(is.matrix(COM) || is.data.frame(COM))
  stopifnot(is.matrix(CAPITAL_TOTAL) || is.data.frame(CAPITAL_TOTAL))
  COM_ts <- as.matrix(COM)
  K_ts   <- as.matrix(CAPITAL_TOTAL)

  Tn <- nrow(Y); S <- ncol(Y)
  T_train <- max(2L, floor(Tn * train_frac))
  if (T_train >= Tn) {
    stop("`train_frac` leaves no observations for the test window; lower it.",
         call. = FALSE)
  }
  T_lik <- if (fit_window == "train") T_train else Tn

  # ---- Level-2/3 specific input: the G'-driven mean needs G'_t ----
  if (n_levels >= 2L) {
    if (is.null(Gprime)) {
      stop("`Gprime` (tasa de ganancia agregada G'_t, largo T) es obligatorio ",
           "para n_levels >= 2.", call. = FALSE)
    }
    Gprime <- as.numeric(Gprime)
    if (length(Gprime) != Tn) {
      stop(sprintf("`Gprime` must have length T = %d (got %d).", Tn, length(Gprime)),
           call. = FALSE)
    }
    if (!all(is.finite(Gprime))) stop("`Gprime` contains non-finite values.", call. = FALSE)
  }

  vmsg(sprintf("Data: T=%d, S=%d, T_train=%d, n_levels=%d, fit_window=%s (T_lik=%d)",
               Tn, S, T_train, n_levels, fit_window, T_lik), verbose)

  # ---- Align COM / CAPITAL_TOTAL columns to Y ----
  COM_ts <- .align_to_Y(COM_ts, Y, S, "COM")
  K_ts   <- .align_to_Y(K_ts,   Y, S, "CAPITAL_TOTAL")

  # ---- Standardize on training statistics ----
  vmsg("Standardizing data using training-period statistics", verbose)
  zY <- zscore_train(Y, T_train)
  zX <- zscore_train(X, T_train)

  mu_tmg <- mean(TMG[seq_len(T_train)], na.rm = TRUE)
  sd_tmg <- stats::sd(TMG[seq_len(T_train)], na.rm = TRUE)
  if (!is.finite(sd_tmg) || sd_tmg < 1e-8) sd_tmg <- 1
  zTMG <- (TMG - mu_tmg) / sd_tmg

  vmsg(sprintf("Computing common factor from %s", factor_from), verbose)
  Mz_factor <- if (factor_from == "X") zX$Mz else zY$Mz
  Ft <- compute_common_factor(Mz_factor, T_train, use_train_loadings, verbose)

  if (orthogonalize_tmg) {
    vmsg("Orthogonalizing TMG with respect to common factor", verbose)
    fit_t <- stats::lm(zTMG[seq_len(T_train)] ~ Ft[seq_len(T_train)])
    zTMG_use <- as.numeric(zTMG - cbind(1, Ft) %*% stats::coef(fit_t))
  } else {
    zTMG_use <- zTMG
  }

  sigma_delta_z <- priors$sigma_delta / sd_tmg
  soft_wedge    <- as.integer(!hard_sum_zero)

  # ---- Assemble Stan data (common block + level-specific block) ----
  stan_dat <- list(
    n_levels = n_levels,
    T = Tn, S = S, T_train = T_train, T_lik = as.integer(T_lik),
    Yz = zY$Mz, Xz = zX$Mz,
    zTMG_byK = as.vector(zTMG_use), zTMG_exo = as.vector(zTMG),
    soft_wedge = soft_wedge, sigma_delta_z = sigma_delta_z,
    COM_ts = COM_ts, K_ts = K_ts,
    com_in_mean = as.integer(isTRUE(com_in_mean)),
    mu_xz = rep(0.0, S),
    beta1_prior_mean = as.numeric(priors$beta1_mean),
    beta1_prior_sd   = as.numeric(priors$beta1_sd),
    nu_prior_shape   = as.numeric(priors$nu_shape),
    nu_prior_rate    = as.numeric(priors$nu_rate),
    rho_prior_mean   = as.numeric(priors$rho_mean),
    rho_prior_sd     = as.numeric(priors$rho_sd),
    sigma_phi_meas_prior_sd = as.numeric(priors$sigma_phi_meas_sd),
    # ---- Measurement-SD latency dial (D-IMPL-9.4). Defaults to the estimated mode
    #      (flag 0, value unused); set to the fixed mode below when the caller
    #      passes sigma_phi_meas_fixed. ----
    sigma_phi_meas_fixed = 0L,
    sigma_phi_meas_value = 0.5,
    # ---- Level-3 value anchor (D-IMPL-10); length-0 unless n_levels == 3, set below. ----
    V_anchor_z       = matrix(0.0, 0L, 0L),
    kappa_cap        = as.numeric(kappa_cap),  # 2-level stability cap (ignored if n=1)
    # ---- Reconstruction-with-uncertain-K fields (recon mode). Length-0 / off by
    #      default; overwritten below when k_uncertainty == "recon". sigma_K_recon
    #      is always a scalar (unused unless k_recon == 1). ----
    k_recon          = 0L,
    k_cost           = matrix(0.0, 0L, 0L),
    K_hat            = matrix(0.0, 0L, 0L),
    Gprime_raw       = numeric(0),
    phi_recon_center = numeric(0),
    phi_recon_scale  = numeric(0),
    sigma_K_recon    = as.numeric(priors$sigma_K_recon)
  )

  # Per-level richness switches (act only when n_levels == 2; canonical otherwise).
  stan_dat <- c(stan_dat, .level_spec_flags(level_spec))

  # Fix the anchor measurement SD to a datum when requested (K-deterministic mode,
  # D-IMPL-9.4): no sigma_phi_meas parameter is sampled, removing the boundary
  # funnel. Only meaningful for n_levels == 2 (the measurement block); harmless in
  # single-level mode (the parameter is length 0 regardless).
  if (!is.null(sigma_phi_meas_fixed)) {
    stan_dat$sigma_phi_meas_fixed <- 1L
    stan_dat$sigma_phi_meas_value <- as.numeric(sigma_phi_meas_fixed)
  }

  if (n_levels >= 2L) {
    mu_g <- mean(Gprime[seq_len(T_train)], na.rm = TRUE)
    sd_g <- stats::sd(Gprime[seq_len(T_train)], na.rm = TRUE)
    if (!is.finite(sd_g) || sd_g < 1e-8) sd_g <- 1
    Gprime_z <- (Gprime - mu_g) / sd_g

    stan_dat$Phi_anchor_z <- zX$Mz                 # constructed Phi index = noisy anchor
    stan_dat$Gprime       <- as.vector(Gprime_z)
    stan_dat$theta_sep    <- as.integer(theta_separation == "hard")

    # ---- Level-3 value anchor (D-IMPL-10): the latent production price reverts
    #      toward a mean that tracks the VALUE index V (direct prices c+v+p). V is
    #      standardized on training statistics (its own mean/sd) so the coupling
    #      m_v is a dimensionless slope; it is a DATUM (direct empirical
    #      construction), never solved simultaneously. ----
    if (n_levels == 3L) {
      V_value <- as.matrix(V_value)
      V_value <- .align_to_Y(V_value, Y, S, "V_value")
      if (nrow(V_value) != Tn) {
        stop(sprintf("`V_value` must have T = %d rows (got %d).", Tn, nrow(V_value)),
             call. = FALSE)
      }
      if (!all(is.finite(V_value))) stop("`V_value` contains non-finite values.", call. = FALSE)
      zV <- zscore_train(V_value, T_train)
      stan_dat$V_anchor_z <- zV$Mz
    }

    if (identical(k_uncertainty, "recon")) {
      # Reconstruction anchor Phi = k + K G' with K lognormal-uncertain (eq. 11).
      # K_hat is the point estimate of total capital advanced = CAPITAL_TOTAL;
      # the standardization (X train mean/sd) maps the reconstructed price into the
      # latent-Phi units so the sigma_K -> 0 limit reproduces the "meas" anchor.
      if (is.null(k_cost)) {
        stop("k_uncertainty = \"recon\" requiere `k_cost` (precio de costo c+v, ",
             "matriz T x S).", call. = FALSE)
      }
      k_cost <- as.matrix(k_cost)
      k_cost <- .align_to_Y(k_cost, Y, S, "k_cost")
      if (nrow(k_cost) != Tn) {
        stop(sprintf("`k_cost` must have T = %d rows (got %d).", Tn, nrow(k_cost)),
             call. = FALSE)
      }
      if (!all(is.finite(k_cost))) stop("`k_cost` contains non-finite values.", call. = FALSE)
      if (any(K_ts <= 0)) {
        stop("`CAPITAL_TOTAL` (the K point estimate) must be strictly positive in ",
             "recon mode (it anchors a lognormal prior on K).", call. = FALSE)
      }
      # Consistency: the constructed index X must equal the deterministic
      # reconstruction k + K_hat G'. Otherwise recon and meas are not on the same
      # scale (recon would not nest meas). Warn with the relative discrepancy.
      recon_det <- k_cost + K_ts * Gprime          # X at K = K_hat (raw units)
      rel_disc  <- max(abs(X - recon_det)) / (max(abs(X)) + 1e-12)
      if (is.finite(rel_disc) && rel_disc > 1e-3) {
        warning(sprintf(paste0(
          "recon anchor: max relative discrepancy between X and (k_cost + ",
          "CAPITAL_TOTAL * Gprime) is %.3g (> 1e-3). The constructed index X ",
          "and the reconstruction are not on the same scale, so recon will not ",
          "reduce to meas as sigma_K -> 0. Check that X = k + K G' (eq. 9)."),
          rel_disc), call. = FALSE)
      }

      stan_dat$k_recon          <- 1L
      stan_dat$k_cost           <- k_cost
      stan_dat$K_hat            <- K_ts
      stan_dat$Gprime_raw       <- as.vector(Gprime)     # RAW G' for the price identity
      stan_dat$phi_recon_center <- as.vector(zX$mu)      # X train mean (standardization)
      stan_dat$phi_recon_scale  <- as.vector(zX$sd)      # X train sd
    }
  } else {
    stan_dat$Phi_anchor_z <- matrix(0.0, 0L, 0L)   # length-0 in single-level mode
    stan_dat$Gprime       <- numeric(0)
    stan_dat$theta_sep    <- 0L
  }

  # ---- Fit ----
  # The 2-level model propagates a latent production-price path, so the default
  # wide random init (radius 2) can land on an explosive trajectory and fail at
  # the initial log density (especially with the Level-2 richness on). A smaller
  # init radius starts every parameter near the prior centre, where the latent
  # recursion is well behaved; warmup then adapts. Single-level keeps the legacy
  # random init. Only applied when the caller did not supply an explicit init.
  if (n_levels == 2L && is.null(init)) init <- 0.3

  fit <- .run_stan_ou(
    stan_dat, chains = chains, iter = iter, warmup = warmup, thin = thin,
    cores = cores, threads_per_chain = threads_per_chain,
    adapt_delta = adapt_delta, max_treedepth = max_treedepth,
    seed = seed, init = init, verbose = verbose
  )

  # ---- Post-processing branches ----
  vmsg("Extracting posterior summaries", verbose)
  summ <- extract_posterior_summary(fit)
  rhat_vec   <- as.numeric(summ$rhat)
  rhat_max   <- max(rhat_vec, na.rm = TRUE)
  rhat_share <- mean(rhat_vec > 1.01, na.rm = TRUE)

  if (isTRUE(moment_match)) {
    vmsg(paste("Note: moment_match is not applied to the array-based PSIS-LOO.",
               "Use loo::loo_moment_match() on the returned stan_fit if needed."),
         verbose)
  }

  vmsg("Computing PSIS-LOO over the fitted window (2:T_lik)", verbose)
  loo_res <- NULL; loo_pareto_k_summary <- NULL
  if (requireNamespace("loo", quietly = TRUE)) {
    loo_res <- tryCatch(.compute_loo(fit, T_lik, S),
                        error = function(e) {
                          warning("PSIS-LOO could not be computed: ",
                                  conditionMessage(e), call. = FALSE); NULL
                        })
    if (!is.null(loo_res)) loo_pareto_k_summary <- .summarize_pareto_k(loo_res)
  }

  vmsg("Computing divergence count", verbose)
  dv <- count_divergences(fit)

  out <- results_robust %||% list()

  if (n_levels == 1L) {
    vmsg("Computing out-of-sample metrics (single-level)", verbose)
    oos <- evaluate_oos(
      summ, zY$Mz, zX$Mz, zTMG_use, T_train,
      COM_ts = COM_ts, K_ts = K_ts,
      com_in_mean = isTRUE(com_in_mean), horizons = c(1, 4, 8)
    )

    out$factor_ou <- list(
      model = "ou_nonlinear_tmg", n_levels = 1L, stan_fit = fit,
      beta1 = summ$beta1, beta0_s = summ$beta0_s, kappa_s = summ$kappa_s,
      a3_s = summ$a3_s, theta_s = summ$theta_s,
      sv = list(alpha = summ$alpha_s, rho = summ$rho_s, sigma_eta = summ$sigma_eta_s),
      nu = summ$nu, gamma = summ$gamma,
      factor_ou_info = list(
        T_train = T_train, T_lik = T_lik, fit_window = fit_window,
        com_in_mean = isTRUE(com_in_mean), factor_from = factor_from,
        use_train_loadings = isTRUE(use_train_loadings)
      )
    )
    out$beta_tmg    <- build_beta_tmg_table(fit, zTMG_use, summ = summ)
    out$sv          <- list(h_summary = summarize_sv_sigmas(fit), rho_s = summ$rho_s)
    out$nonlinear   <- list(a3 = summ$a3_s, drift_decomp = drift_decomposition_grid(fit, summ))
    out$accounting  <- build_accounting_block(TMG, zTMG, zTMG_use, mu_tmg, sd_tmg,
                                              hard_sum_zero, priors$sigma_delta)
    out$diagnostics <- list(
      rhat = summ$rhat, ess = summ$ess, rhat_max = rhat_max,
      rhat_share = rhat_share, divergences = dv,
      loo = loo_res, loo_pareto_k = loo_pareto_k_summary, oos = oos
    )
    class(out$factor_ou) <- c("ou_nonlinear_tmg", "list")

  } else {
    vmsg("Extracting Level-2 (latent production) summaries", verbose)
    lvl2 <- .extract_level2_summary(fit, S, Tn, colnames(Y))
    # In the fixed-SD mode (D-IMPL-9.4) sigma_phi_meas is a datum, not a draw, so
    # the extractor returns NA; report the known fixed value instead.
    if (!is.null(sigma_phi_meas_fixed)) {
      lvl2$sigma_phi_meas <- as.numeric(sigma_phi_meas_fixed)
    }
    phi  <- .extract_phi_latent(fit, Tn, S, colnames(Y))
    sep  <- .separation_evidence(fit, S)
    mu_path <- extract_mu_trajectory(fit, Gprime_z, colnames(Y),
                                     V_z = if (n_levels == 3L) zV$Mz else NULL)

    vmsg("Computing out-of-sample metrics (2-level, phi reverts to latent Phi)",
         verbose)
    meds <- list(
      kappa_tilde = .median_vec(fit, "kappa_tilde"),
      a3_s        = summ$a3_s,
      beta1       = summ$beta1,
      gamma       = summ$gamma,
      kappa_p     = lvl2$kappa_p[, "median"],
      mu_const    = lvl2$mu_const[, "median"],
      m1          = lvl2$m1,
      m_v         = lvl2$m_v          # scalar; NA unless n_levels == 3 (Level 3)
    )
    oos <- tryCatch(
      evaluate_oos_nested(meds, zY$Mz, phi$median, Gprime_z, zTMG_use, T_train,
                          COM_ts = COM_ts, K_ts = K_ts, kappa_cap = kappa_cap,
                          com_in_mean = isTRUE(com_in_mean), horizons = c(1, 4, 8),
                          V_z = if (n_levels == 3L) zV$Mz else NULL),
      error = function(e) {
        warning("2-level OOS recursion failed: ", conditionMessage(e),
                call. = FALSE); NULL
      }
    )

    out$factor_ou <- list(
      model = if (n_levels == 3L) "ou_nested_3level" else "ou_nested_2level",
      n_levels = n_levels, stan_fit = fit,
      beta1 = summ$beta1, kappa_s = summ$kappa_s, a3_s = summ$a3_s,
      sv = list(alpha = summ$alpha_s, rho = summ$rho_s, sigma_eta = summ$sigma_eta_s),
      nu = summ$nu, gamma = summ$gamma,
      level2 = lvl2,
      factor_ou_info = list(
        T_train = T_train, T_lik = T_lik, fit_window = fit_window,
        com_in_mean = isTRUE(com_in_mean), factor_from = factor_from,
        use_train_loadings = isTRUE(use_train_loadings),
        theta_separation = theta_separation, k_uncertainty = k_uncertainty,
        sigma_phi_meas_fixed = if (!is.null(sigma_phi_meas_fixed))
          as.numeric(sigma_phi_meas_fixed) else NULL,
        kappa_cap = kappa_cap, level_spec = level_spec
      )
    )
    out$phi_latent  <- phi
    out$mu_path     <- mu_path
    out$separation  <- sep
    out$sv          <- list(h_summary = summarize_sv_sigmas(fit), rho_s = summ$rho_s)
    out$accounting  <- build_accounting_block(TMG, zTMG, zTMG_use, mu_tmg, sd_tmg,
                                              hard_sum_zero, priors$sigma_delta)
    out$diagnostics <- list(
      rhat = summ$rhat, ess = summ$ess, rhat_max = rhat_max,
      rhat_share = rhat_share, divergences = dv,
      loo = loo_res, loo_pareto_k = loo_pareto_k_summary,
      oos = oos
    )
    # The S3 class lives on the top-level result (the object the user holds), so
    # print/summary/plot dispatch naturally; the inner factor_ou keeps the plain
    # `model` string tag to avoid a duplicate class on a nested element.
    class(out) <- c("ou_nested_2level", "list")
  }

  vmsg("Model fitting complete", verbose)
  out
}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' Resolve and validate the per-level richness specification
#'
#' Each level toggles four richness dimensions: \code{cubic} (cubic drift),
#' \code{sv} (stochastic volatility), \code{student_t} (Student-t tails) and
#' \code{hierarchy} (cross-sector partial pooling). \code{NULL} selects the
#' canonical configuration (Level 1 full, Level 2 lean). The switches act only in
#' the 2-level mode; the single-level model is always the full legacy Level-1
#' specification, so for \code{n_levels = 1} only \code{NULL} (or a Level-1-full
#' spec) is accepted.
#'
#' @keywords internal
#' @noRd
.canonical_level_spec <- function() {
  list(
    level1 = list(cubic = TRUE,  sv = TRUE,  student_t = TRUE,  hierarchy = TRUE),
    level2 = list(cubic = FALSE, sv = FALSE, student_t = FALSE, hierarchy = TRUE)
  )
}

.resolve_level_spec <- function(level_spec, n_levels) {
  if (is.null(level_spec)) return(.canonical_level_spec())

  fields <- c("cubic", "sv", "student_t", "hierarchy")
  check_level <- function(L, nm) {
    if (!is.list(L) || !all(fields %in% names(L))) {
      stop(sprintf("`level_spec$%s` must be a list with logical entries %s.",
                   nm, paste(fields, collapse = ", ")), call. = FALSE)
    }
    vals <- vapply(fields, function(f) {
      v <- L[[f]]
      if (length(v) != 1L || !is.logical(v) || is.na(v)) {
        stop(sprintf("`level_spec$%s$%s` must be a single TRUE/FALSE.", nm, f),
             call. = FALSE)
      }
      isTRUE(v)
    }, logical(1))
    stats::setNames(as.list(vals), fields)
  }
  if (!is.list(level_spec) || !all(c("level1", "level2") %in% names(level_spec))) {
    stop("`level_spec` must be a list with `level1` and `level2` entries (use ",
         "ou_level_spec() for the named experiment configurations).", call. = FALSE)
  }
  spec <- list(level1 = check_level(level_spec$level1, "level1"),
               level2 = check_level(level_spec$level2, "level2"))

  if (n_levels == 1L && !all(unlist(spec$level1))) {
    stop("In single-level mode (n_levels = 1) Level 1 is always the full legacy ",
         "specification; pass level_spec = NULL.", call. = FALSE)
  }
  spec
}

#' Translate a resolved level_spec into the 8 integer Stan switches
#' @keywords internal
#' @noRd
.level_spec_flags <- function(spec) {
  i <- function(x) as.integer(isTRUE(x))
  list(
    l1_cubic = i(spec$level1$cubic), l1_sv = i(spec$level1$sv),
    l1_studentt = i(spec$level1$student_t), l1_hier = i(spec$level1$hierarchy),
    l2_cubic = i(spec$level2$cubic), l2_sv = i(spec$level2$sv),
    l2_studentt = i(spec$level2$student_t), l2_hier = i(spec$level2$hierarchy)
  )
}

#' Align a (T x S) covariate matrix to the columns of Y
#' @keywords internal
#' @noRd
.align_to_Y <- function(M, Y, S, label) {
  if (!is.null(colnames(Y)) && !is.null(colnames(M))) {
    common <- intersect(colnames(Y), colnames(M))
    if (length(common) != S) {
      missing <- setdiff(colnames(Y), colnames(M))
      extra   <- setdiff(colnames(M), colnames(Y))
      stop(sprintf(
        "Column mismatch %s vs Y. Missing in %s: %s. Extra in %s: %s",
        label, label,
        if (length(missing) == 0) "(none)" else paste(missing, collapse = ", "),
        label,
        if (length(extra)   == 0) "(none)" else paste(extra,   collapse = ", ")
      ), call. = FALSE)
    }
    return(M[, colnames(Y), drop = FALSE])
  }
  if (ncol(M) != S) stop(sprintf("Dimension mismatch %s vs Y.", label), call. = FALSE)
  M
}

#' Compile and sample the unified Stan model (cmdstanr or rstan)
#' @keywords internal
#' @noRd
.run_stan_ou <- function(stan_dat, chains, iter, warmup, thin, cores,
                         threads_per_chain, adapt_delta, max_treedepth,
                         seed, init, verbose) {
  backend <- check_stan_backend(verbose)
  if (backend == "none") {
    stop("Stan backend required. Please install cmdstanr or rstan.", call. = FALSE)
  }

  par_chains    <- max(1L, min(as.integer(chains), as.integer(cores)))
  thr_per_chain <- max(1L, min(as.integer(threads_per_chain),
                               as.integer(floor(cores / par_chains))))

  if (backend == "cmdstanr") {
    vmsg("Compiling Stan model with cmdstanr (canonical .stan file)", verbose)
    cache_dir <- tryCatch({
      d <- tools::R_user_dir("bayesianOU", "cache")
      dir.create(d, recursive = TRUE, showWarnings = FALSE)
      if (dir.exists(d)) d else tempdir()
    }, error = function(e) tempdir())
    mod <- cmdstanr::cmdstan_model(.stan_file_path(), dir = cache_dir,
                                   pedantic = FALSE,
                                   cpp_options = list(stan_threads = TRUE))
    vmsg("Running MCMC sampling", verbose)
    return(mod$sample(
      data = stan_dat, chains = chains, parallel_chains = par_chains,
      iter_warmup = warmup, iter_sampling = iter - warmup, thin = thin,
      seed = seed, refresh = if (verbose) 200 else 0,
      adapt_delta = adapt_delta, max_treedepth = max_treedepth,
      threads_per_chain = thr_per_chain, init = init
    ))
  }

  vmsg("Compiling Stan model with rstan", verbose)
  sm <- rstan::stan_model(model_code = ou_nested_stan_code())
  vmsg("Running MCMC sampling", verbose)
  # A single numeric init is an init RADIUS in this package's convention
  # (cmdstanr semantics). rstan controls the radius via init_r with
  # init = "random", so translate accordingly; lists/functions pass through.
  rstan_init <- "random"; rstan_init_r <- 2
  if (!is.null(init)) {
    if (is.numeric(init) && length(init) == 1L) rstan_init_r <- init
    else rstan_init <- init
  }
  rstan::sampling(
    sm, data = stan_dat, chains = chains, iter = iter, warmup = warmup,
    thin = thin, seed = seed,
    control = list(adapt_delta = adapt_delta, max_treedepth = max_treedepth),
    refresh = if (verbose) 200 else 0, init = rstan_init, init_r = rstan_init_r
  )
}

#' Extract Level-2 (latent production) parameter summaries
#'
#' Returns the Level-2 parameter bands and, when the corresponding richness
#' switch is on, the extra-richness quantities of the latent production price:
#' the cubic restoring coefficient \code{a3_p} (\code{l2_cubic}), the Student-t
#' degrees of freedom \code{nu_p} (\code{l2_studentt}) and the time-varying
#' stochastic-volatility scale path \code{sigma_p_t} (\code{l2_sv}). The last two
#' are derived from the posterior draws (\eqn{\nu_p = 2 + \tilde\nu_p};
#' \eqn{\sigma_p(t) = \exp(h_p/2)}) so no Stan-side change is required, and they
#' are \code{NULL} when the switch is off (the canonical configuration is
#' unchanged). See decision D-IMPL-10.
#' @keywords internal
#' @noRd
.extract_level2_summary <- function(fit, S, Tn = NULL, sector_names = NULL) {
  med_ci <- function(p) {
    M <- .draws_matrix(fit, p)
    if (is.null(M)) return(NULL)
    qs <- t(apply(M, 2, stats::quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE))
    colnames(qs) <- c("q2.5", "median", "q97.5"); qs
  }
  scal <- function(p) {
    M <- .draws_matrix(fit, p)
    if (is.null(M)) return(NA_real_)
    stats::median(M, na.rm = TRUE)
  }
  # Level-2 Student-t df nu_p = 2 + nu_p_tilde (present only when l2_studentt is
  # on; NULL otherwise). One scalar, returned as a 1-row median/CI band.
  nu_p <- {
    M <- .draws_matrix(fit, "nu_p_tilde")
    if (is.null(M)) NULL else {
      q <- stats::quantile(2 + M[, 1], probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
      matrix(q, nrow = 1, dimnames = list(NULL, c("q2.5", "median", "q97.5")))
    }
  }
  list(
    kappa_p        = med_ci("kappa_p"),
    kappa_m_base   = med_ci("kappa_m_base"),
    mu_const       = med_ci("mu_const"),
    sigma_p        = med_ci("sigma_p"),
    m1             = scal("m1"),
    # Level-3 value coupling (D-IMPL-10): NA unless n_levels == 3 (no m_v draws).
    m_v            = scal("m_v"),
    sigma_phi_meas = scal("sigma_phi_meas"),
    # Level-2 richness (NULL unless the corresponding switch is on): cubic
    # restoring coefficient, Student-t df, and the SV scale path (D-IMPL-10).
    a3_p           = med_ci("a3_p"),
    nu_p           = nu_p,
    sigma_p_t      = .extract_sigma_p_path(fit, Tn, S, sector_names)
  )
}

#' Time-varying Level-2 SV scale path sigma_p(t) = exp(h_p / 2)
#'
#' Median and 95\% bands of the latent production-price stochastic-volatility
#' scale, present only when \code{l2_sv} is on (the draws of \code{h_p} exist).
#' Mirrors \code{\link{.extract_phi_latent}} but maps each log-variance draw to a
#' scale. \code{NULL} when SV is off or \code{Tn} is unknown.
#' @keywords internal
#' @noRd
.extract_sigma_p_path <- function(fit, Tn, S, sector_names = NULL) {
  if (is.null(Tn)) return(NULL)
  M <- .draws_matrix(fit, "h_p")           # columns "h_p[t,s]"; NULL when l2_sv off
  if (is.null(M)) return(NULL)
  cn  <- colnames(M)
  tt  <- as.integer(sub("^h_p\\[(\\d+),(\\d+)\\]$", "\\1", cn))
  ss  <- as.integer(sub("^h_p\\[(\\d+),(\\d+)\\]$", "\\2", cn))
  sig <- exp(0.5 * M)                       # SV scale draws
  q   <- apply(sig, 2, stats::quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
  med <- matrix(NA_real_, Tn, S); lo <- med; hi <- med
  for (j in seq_along(cn)) {
    med[tt[j], ss[j]] <- q[2, j]; lo[tt[j], ss[j]] <- q[1, j]; hi[tt[j], ss[j]] <- q[3, j]
  }
  if (!is.null(sector_names) && length(sector_names) == S) {
    colnames(med) <- colnames(lo) <- colnames(hi) <- sector_names
  }
  list(median = med, q2.5 = lo, q97.5 = hi)
}

#' Extract the latent production-price trajectory Phi (median + 95% bands)
#' @keywords internal
#' @noRd
.extract_phi_latent <- function(fit, Tn, S, sector_names = NULL) {
  M <- .draws_matrix(fit, "Phi")           # columns "Phi[t,s]"
  if (is.null(M)) return(NULL)
  cn  <- colnames(M)
  tt  <- as.integer(sub("^Phi\\[(\\d+),(\\d+)\\]$", "\\1", cn))
  ss  <- as.integer(sub("^Phi\\[(\\d+),(\\d+)\\]$", "\\2", cn))
  med <- matrix(NA_real_, Tn, S); lo <- med; hi <- med
  q   <- apply(M, 2, stats::quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
  for (j in seq_along(cn)) {
    med[tt[j], ss[j]] <- q[2, j]; lo[tt[j], ss[j]] <- q[1, j]; hi[tt[j], ss[j]] <- q[3, j]
  }
  if (!is.null(sector_names) && length(sector_names) == S) {
    colnames(med) <- colnames(lo) <- colnames(hi) <- sector_names
  }
  list(median = med, q2.5 = lo, q97.5 = hi)
}

#' Posterior evidence for the time-scale separation kappa^m > kappa^p
#'
#' Compares the bounded market and production reversion speeds at zTMG = 0
#' (\code{kappa_m_base} vs \code{kappa_p}), both in (0, kappa_cap).
#' @keywords internal
#' @noRd
.separation_evidence <- function(fit, S) {
  Km <- .draws_matrix(fit, "kappa_m_base") # bounded market speed at zTMG = 0
  Kp <- .draws_matrix(fit, "kappa_p")      # bounded production speed
  if (is.null(Km) || is.null(Kp)) return(NULL)
  per_sector <- vapply(seq_len(S), function(s) mean(Km[, s] > Kp[, s]), numeric(1))
  joint      <- mean(apply(Km > Kp, 1, all))
  list(prob_sep_by_sector = per_sector, prob_sep_joint = joint)
}

#' Draws as a [draws x columns] matrix for a parameter, backend-agnostic
#' @keywords internal
#' @noRd
.draws_matrix <- function(fit, p) {
  if (inherits(fit, "CmdStanMCMC")) {
    df <- tryCatch(fit$draws(p, format = "df"), error = function(e) NULL)
    if (is.null(df)) return(NULL)
    keep <- grep(sprintf("^%s(\\[|$)", p), names(df), value = TRUE)
    keep <- setdiff(keep, c(".chain", ".iteration", ".draw"))
    if (length(keep) == 0) return(NULL)
    return(as.matrix(df[, keep, drop = FALSE]))
  }
  out <- tryCatch(rstan::extract(fit, pars = p, permuted = TRUE)[[1]],
                  error = function(e) NULL)
  if (is.null(out)) return(NULL)
  if (is.null(dim(out))) return(matrix(out, ncol = 1))
  if (length(dim(out)) == 2L) return(out)
  matrix(out, nrow = dim(out)[1])          # flatten higher-dim (e.g. Phi[draw,t,s])
}
