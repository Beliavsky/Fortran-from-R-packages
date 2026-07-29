// =============================================================================
// Unified nonlinear Ornstein-Uhlenbeck model: single-level and 2-level nested.
//
// SINGLE SOURCE OF TRUTH. This file subsumes the former ou_nonlinear_tmg.stan:
//   - n_levels == 1  : EXACT legacy single-level model (market price reverts to a
//                      constant theta_s, with a linear forcing by the production
//                      price X and a TMG interaction). Bit-for-bit equivalent to
//                      the retired ou_nonlinear_tmg.stan (certified by a test).
//   - n_levels >= 2  : nested model. Market price phi (= Yz) reverts to a LATENT
//                      production price Phi (Level 1); Phi has its own OU reverting
//                      to a G'-driven mean mu (Level 2); the constructed production
//                      -price index (Phi_anchor_z) enters as a noisy measurement of
//                      the latent Phi. The aggregate profit rate (TMG) modulates the
//                      Level-1 reversion speed; G' drives the Level-2 mean.
//
// All Level-2 parameters are sized to length 0 when n_levels == 1, so the
// single-level mode carries no extra dimensions (clean and exactly the old model).
//
// Discretization: Euler-Maruyama with dt = 1 (see R documentation for half-lives).
// Likelihood window: summed over t in 2..T_lik (train-only or full-information).
// =============================================================================

functions {
  // One-step-ahead Student-t log density for the MARKET price phi (= Yz), shared
  // by both modes. The mean equation branches on n_levels; everything else (SV
  // scale, COM-in-mean, TMG wedge) is common.
  real ou_nested_partial_sum(array[] int t_idx_slice,
                             int start, int end,
                             int    n_levels,
                             int    T_lik,
                             matrix Yz,
                             matrix Xz,
                             matrix Phi,            // latent Phi (n_levels==2); Xz passed otherwise
                             matrix COM_ts,
                             vector zTMG_byK,
                             int    soft_wedge,
                             vector delta_z,
                             vector com_wmean_train,
                             vector com_wsd_train,
                             vector mu_xz,
                             vector theta_s,
                             vector kappa_s,
                             vector kappa_tilde,
                             real   kappa_cap,
                             vector a3_s,
                             vector beta0_s,
                             real   beta1,
                             matrix h_eff,
                             real   nu,
                             int    com_in_mean,
                             real   gamma,
                             int    l1_cubic,
                             int    l1_studentt) {
    real lp = 0;
    int S = cols(Yz);
    for (t_idx in t_idx_slice) {
      int t = t_idx;
      if (t <= 1) continue;
      if (t > T_lik) continue;
      real ztmg_eff = zTMG_byK[t];
      if (soft_wedge == 1) ztmg_eff += delta_z[t];
      ztmg_eff = fmin(fmax(ztmg_eff, -1e6), 1e6);
      for (s in 1:S) {
        real denom_sd = com_wsd_train[s];
        denom_sd = (denom_sd > 1e-12) ? denom_sd : 1.0;
        real com_std = (COM_ts[t-1,s] - com_wmean_train[s]) / denom_sd;
        com_std = fmin(fmax(com_std, -1e6), 1e6);
        real com_term = (com_in_mean == 1) ? gamma * com_std : 0;
        real sd_safe  = fmin(fmax(exp(0.5 * h_eff[t,s]), 1e-8), 1e8);

        real mean_;
        if (n_levels == 1) {
          // ---- Legacy single-level mean (reversion to constant theta_s) ----
          real zlag  = Yz[t-1,s] - theta_s[s];
          real drift = kappa_s[s] * (theta_s[s] - Yz[t-1,s] + a3_s[s] * zlag^3);
          real betaT = beta0_s[s] + beta1 * ztmg_eff;
          mean_ = drift + betaT * (Xz[t-1,s] - mu_xz[s]) + com_term;
        } else {
          // ---- Nested Level 1: market reverts to latent production Phi ----
          // dev = phi - Phi; cubic restoring on dev (a3_s < 0, gated by l1_cubic).
          // The gravitation speed is bounded to (0, kappa_cap) via a logit link
          // with the TMG term INSIDE the link, so neither the level nor the TMG
          // modulation can push the speed out of the stable Euler region (this
          // prevents the latent path from exploding). beta1 = 0 -> no modulation.
          real dev     = Yz[t-1,s] - Phi[t-1,s];
          real kappa_m = kappa_cap * inv_logit(kappa_tilde[s] + beta1 * ztmg_eff);
          real cub     = (l1_cubic == 1) ? a3_s[s] * dev^3 : 0;
          mean_ = kappa_m * (-dev + cub) + com_term;
        }

        real resid = Yz[t,s] - Yz[t-1,s] - mean_;
        // Student-t tails are always on in single-level mode; in the 2-level mode
        // they are gated by l1_studentt (off -> Gaussian innovations).
        if (n_levels >= 2 && l1_studentt == 0)
          lp += normal_lpdf(resid | 0, sd_safe);
        else
          lp += student_t_lpdf(resid | nu, 0, sd_safe);
      }
    }
    return lp;
  }
}

data {
  // 1 = single-level (legacy); 2 = nested (market -> latent production);
  // 3 = nested + values (Level 3, D-IMPL-10): the latent production price reverts
  // toward a mean derived from the VALUE anchor V (direct prices c+v+p), so prices
  // of production gravitate around values (Marx, Capital III ch. IX). n_levels == 3
  // inherits ALL of the n_levels >= 2 machinery (latent Phi, its OU, the
  // measurement anchor) and adds the value coupling m_v * V; the Level-2 gates were
  // widened from "== 2" to ">= 2" for exactly this. Values enter by DIRECT empirical
  // construction (V is observed data), never by simultaneist/Leontief inversion.
  int<lower=1, upper=3> n_levels;

  int<lower=2> T;
  int<lower=1> S;
  int<lower=2> T_train;
  int<lower=2, upper=T> T_lik;
  matrix[T,S] Yz;                       // market price phi (standardized)
  matrix[T,S] Xz;                       // production price (single-level forcing)
  vector[T] zTMG_byK;
  vector[T] zTMG_exo;
  int<lower=0,upper=1> soft_wedge;
  real<lower=0> sigma_delta_z;
  matrix[T,S] COM_ts;
  matrix[T,S] K_ts;
  int<lower=0,upper=1> com_in_mean;
  vector[S] mu_xz;

  // Configurable priors (defaults set on the R side).
  real beta1_prior_mean;
  real<lower=0> beta1_prior_sd;
  real<lower=0> nu_prior_shape;
  real<lower=0> nu_prior_rate;
  real rho_prior_mean;
  real<lower=0> rho_prior_sd;

  // ---- Level-2 (nested) data; ignored when n_levels == 1 ----
  matrix[n_levels >= 2 ? T : 0, n_levels >= 2 ? S : 0] Phi_anchor_z; // constructed Phi index
  vector[n_levels >= 2 ? T : 0] Gprime;          // standardized aggregate profit rate

  // ---- Level-3 (values) data; length 0 unless n_levels == 3 (D-IMPL-10). The
  //      VALUE anchor V is the standardized direct-price index (c+v+p =
  //      DirectPrices_Index, constructed directly from k_cost + EBO). The latent
  //      production price reverts toward a mean that tracks V with coupling m_v
  //      (production prices gravitate around values, Capital III ch. IX). V is a
  //      datum (direct empirical construction), NOT solved simultaneously. ----
  matrix[n_levels == 3 ? T : 0, n_levels == 3 ? S : 0] V_anchor_z;   // value index (standardized)
  int<lower=0,upper=1> theta_sep;                // 1 = hard ordering kappa_p < kappa_m
  real<lower=0> sigma_phi_meas_prior_sd;         // half-normal scale for the anchor sd (estimated mode)
  // ---- Measurement-SD latency dial (D-IMPL-9.4). In K-deterministic mode the
  //      anchor measurement SD is the sigma_meas -> 0 LIMIT done well: a FIXED
  //      datum rather than an estimated parameter, which removes the boundary
  //      funnel (G6_boundary) the estimated sigma_phi_meas pinned against 0 in the
  //      deterministic limit (rhat 2.208 in the Session-9 CP smoke). When
  //      sigma_phi_meas_fixed == 1 the measurement SD equals sigma_phi_meas_value
  //      and no parameter is sampled; when 0 it is estimated with the half-normal
  //      prior above (K-stochastic / recon, looser prior). Ignored when
  //      n_levels == 1 (the measurement block is inside the n_levels >= 2 guard). ----
  int<lower=0,upper=1> sigma_phi_meas_fixed;     // 1 = fixed datum (K-det); 0 = estimated (K-est)
  real<lower=0> sigma_phi_meas_value;            // the fixed measurement SD (used iff fixed == 1)
  real<lower=0> kappa_cap;                       // stability cap for the 2-level reversion
                                                 // speeds: kappa_m, kappa_p in (0, kappa_cap).
                                                 // The Euler OU is stable for kappa < 2 and
                                                 // monotone (no overshoot) for kappa < 1.
                                                 // Unused when n_levels == 1.

  // ---- Per-level richness switches (level_spec). Only act in the 2-level mode
  //      (n_levels >= 2); the single-level branch is always the legacy full L1
  //      specification, so n_levels == 1 is untouched. Level-1 switches gate
  //      terms that already exist (cubic drift, stochastic volatility, Student-t
  //      tails, cross-sector hierarchy); the canonical configuration is L1 = all
  //      on. Level-2 switches enable extra richness on the latent production
  //      price; the canonical configuration is L2 = cubic/SV/Student-t off,
  //      hierarchy on. ----
  int<lower=0,upper=1> l1_cubic;
  int<lower=0,upper=1> l1_sv;
  int<lower=0,upper=1> l1_studentt;
  int<lower=0,upper=1> l1_hier;
  int<lower=0,upper=1> l2_cubic;
  int<lower=0,upper=1> l2_sv;
  int<lower=0,upper=1> l2_studentt;
  int<lower=0,upper=1> l2_hier;

  // ---- Reconstruction with uncertain K (eq. 11 / Hallazgo 2); k_uncertainty
  //      = "recon". When k_recon == 1 (only valid with n_levels >= 2) the anchor
  //      for the latent Phi is NOT the fixed constructed index but the
  //      reconstruction Phi = k + K G', with the total capital advanced K treated
  //      as uncertain (lognormal prior f_Y around the point estimate K_hat). The
  //      standardization (center/scale) maps the reconstructed price into the
  //      same units as the latent Phi, so the sigma_K_recon -> 0 limit reproduces
  //      the "meas" anchor exactly (recon nests meas). All fields are length 0
  //      unless k_recon == 1. ----
  int<lower=0,upper=1> k_recon;
  matrix[k_recon == 1 ? T : 0, k_recon == 1 ? S : 0] k_cost;   // cost price k = c+v (raw)
  matrix[k_recon == 1 ? T : 0, k_recon == 1 ? S : 0] K_hat;    // point estimate of K (raw, > 0)
  vector[k_recon == 1 ? T : 0] Gprime_raw;                     // raw aggregate profit rate G'_t
  vector[k_recon == 1 ? S : 0] phi_recon_center;               // standardization center (X train mean)
  vector[k_recon == 1 ? S : 0] phi_recon_scale;                // standardization scale (X train sd)
  real<lower=0> sigma_K_recon;                                 // lognormal scale of the K-uncertainty (f_Y)
}

transformed data {
  vector[S] com_wmean_train;
  vector[S] com_wsd_train;
  vector[S] COM_s;

  for (s in 1:S) {
    real denom = 0;
    for (t in 1:T_train) denom += K_ts[t, s];
    if (denom <= 0) denom = 1;
    {
      real num = 0;
      for (t in 1:T_train) num += COM_ts[t, s] * (K_ts[t, s] / denom);
      com_wmean_train[s] = num;
    }
    {
      real v = 0;
      for (t in 1:T_train) {
        real wt = K_ts[t, s] / denom;
        v += wt * square(COM_ts[t, s] - com_wmean_train[s]);
      }
      com_wsd_train[s] = sqrt(fmax(v, 1e-16));
    }
  }

  {
    real muS = mean(com_wmean_train);
    real sdS = sd(com_wmean_train);
    if (sdS <= 1e-8) sdS = 1.0;
    for (s in 1:S) COM_s[s] = (com_wmean_train[s] - muS) / sdS;
  }

  array[T_lik] int t_idx;
  for (t in 1:T_lik) t_idx[t] = t;
  int grainsize = 1;
}

parameters {
  // ---- Shared Level-1 structure (non-centered hierarchy) ----
  vector[S] theta_z;
  vector[S] kappa_z;
  vector[S] a3_z;
  vector[S] beta0_z;

  real theta0;     real theta_COM;    real<lower=1e-6> sigma_theta;
  real kappa0;     real kappa_COM;    real<lower=1e-6> sigma_kappa;
  real a3_0;                          real<lower=1e-6> sigma_a3;
  real beta00;     real beta0_COM;    real<lower=1e-6> sigma_beta0;
  real beta1;

  // SV and fat tails (shared)
  vector[S] alpha_s;
  vector<lower=-0.995, upper=0.995>[S] rho_s;
  vector<lower=1e-6>[S] sigma_eta_s;
  matrix[T,S] h_raw;
  real<lower=1e-6> nu_tilde;

  vector[soft_wedge == 1 ? T : 0] delta_z;
  real gamma;

  // ---- Level-2 (nested) parameters; length 0 when n_levels == 1 ----
  // CENTERED parametrization (CP): the latent production price is sampled
  // DIRECTLY as Phi_lat, not reconstructed from non-centered innovations. This
  // removes the multiplicative non-identification ridge sigma_p * Phi_innov that
  // the former NCP factorization induced under an informative anchor (diagnosed
  // on the real 37-sector panel: cor(sigma_p_tilde, log||Phi_innov||^2) ~ -0.98,
  // E-BFMI ~0.03-0.06, rhat 1.5-2.0 in the failing trio Phi_innov / sigma_p_tilde
  // / m0_z, while Level 1 always converged). The data identify the increment
  // product, not its scale * innovation factorisation; CP samples the increment
  // directly, so its scale is identified by the residual variance instead. The
  // first row Phi_lat[1, ] is anchored by its own prior (there is no separate
  // Phi1). See decision D-IMPL-9.1.
  matrix[n_levels >= 2 ? T : 0, n_levels >= 2 ? S : 0] Phi_lat;  // latent production price (CP)
  vector[n_levels >= 2 ? S : 0] kappa_p_z;                  // Level-2 reversion (hierarchy)
  vector[n_levels >= 2 ? S : 0] m0_z;                       // Level-2 mean intercept (hierarchy)
  vector[n_levels >= 2 ? S : 0] sigma_p_tilde;              // log Level-2 innovation scale
  array[n_levels >= 2 ? 1 : 0] real kappa_p0;
  array[n_levels >= 2 ? 1 : 0] real kappa_p_COM;
  array[n_levels >= 2 ? 1 : 0] real<lower=1e-6> sigma_kappa_p;
  array[n_levels >= 2 ? 1 : 0] real m0_0;
  array[n_levels >= 2 ? 1 : 0] real m0_COM;
  array[n_levels >= 2 ? 1 : 0] real<lower=1e-6> sigma_m0;
  array[n_levels >= 2 ? 1 : 0] real m1;                     // G'-driven mean slope
  // Level-3 value coupling (D-IMPL-10): the production-price reversion mean gains
  // m_v * V (the standardized value anchor). Free and FALSIFIABLE (neutral prior
  // centered at 0, like m1); the Capital III ch. IX hypothesis is m_v > 0 (prices
  // of production gravitate around values), adjudicated by the data, not imposed.
  // Length 0 unless n_levels == 3, so n_levels in {1,2} is bit-for-bit unchanged.
  array[n_levels == 3 ? 1 : 0] real m_v;                    // value-coupling slope (Level 3)
  array[n_levels >= 2 ? 1 : 0] real sep0;                   // hard-separation offset
  // Measurement SD: a parameter ONLY when estimated (sigma_phi_meas_fixed == 0);
  // sized to 0 in the K-deterministic fixed mode, where the SD is a datum
  // (sigma_phi_meas_value), so the boundary funnel disappears (D-IMPL-9.4).
  array[(n_levels >= 2 && sigma_phi_meas_fixed == 0) ? 1 : 0] real<lower=1e-6> sigma_phi_meas;

  // ---- Level-2 richness (production price Phi); each block has length 0 unless
  //      its switch is on, so the canonical L2 (all off) reproduces the linear
  //      Gaussian OU exactly. ----
  vector[l2_cubic == 1 ? S : 0] a3_p_tilde;                 // cubic restoring on Phi (a3_p < 0)
  vector[l2_sv == 1 ? S : 0] alpha_p;                       // L2 SV: log-variance level
  vector<lower=-0.995, upper=0.995>[l2_sv == 1 ? S : 0] rho_p;
  vector<lower=1e-6>[l2_sv == 1 ? S : 0] sigma_eta_p;
  matrix[l2_sv == 1 ? T : 0, l2_sv == 1 ? S : 0] h_p_raw;
  array[l2_studentt == 1 ? 1 : 0] real<lower=1e-6> nu_p_tilde;  // L2 Student-t df shift

  // Non-centered lognormal for the uncertain total capital K (recon mode only).
  matrix[k_recon == 1 ? T : 0, k_recon == 1 ? S : 0] z_K;
}

transformed parameters {
  // Level-1 sector parameters from the non-centered hierarchy.
  vector[S] theta_s;
  vector[S] kappa_tilde;
  vector[S] a3_tilde;
  vector[S] beta0_s;

  for (s in 1:S) {
    if (n_levels >= 2 && l1_hier == 0) {
      // Lean L1 hierarchy: independent sector effects at fixed priors (no
      // cross-sector partial pooling). The hyperparameters become prior-only.
      theta_s[s]     = theta_z[s];
      kappa_tilde[s] = -1.0      + 0.5 * kappa_z[s];
      a3_tilde[s]    = log(0.05) + 0.4 * a3_z[s];
      beta0_s[s]     =             0.5 * beta0_z[s];
    } else {
      theta_s[s]     = theta0 + theta_COM * COM_s[s] + sigma_theta * theta_z[s];
      kappa_tilde[s] = kappa0 + kappa_COM * COM_s[s] + sigma_kappa * kappa_z[s];
      a3_tilde[s]    = a3_0                          + sigma_a3   * a3_z[s];
      beta0_s[s]     = beta00 + beta0_COM * COM_s[s] + sigma_beta0 * beta0_z[s];
    }
  }

  vector<lower=0>[S] kappa_s = exp(kappa_tilde);
  vector<upper=0>[S] a3_s    = -exp(a3_tilde);
  real<lower=2> nu = 2 + nu_tilde;

  // Stationary non-centered AR(1) log-variance (shared).
  matrix[T,S] h;
  matrix[T,S] h_std;
  for (s in 1:S) {
    h_std[1,s] = h_raw[1,s] / sqrt(1 - square(rho_s[s]) + 1e-8);
    for (t in 2:T) h_std[t,s] = rho_s[s] * h_std[t-1,s] + h_raw[t,s];
  }
  for (t in 1:T) for (s in 1:S) h[t,s] = alpha_s[s] + sigma_eta_s[s] * h_std[t,s];

  // Effective Level-1 log-variance fed to the market likelihood. Stochastic
  // volatility (AR(1)) when l1_sv == 1 (always so in single-level mode); a
  // constant per-sector level (alpha_s) for the 2-level "lean" L1. The full SV
  // machinery still runs but becomes prior-only when switched off.
  matrix[T,S] h_eff;
  if (n_levels >= 2 && l1_sv == 0) {
    for (t in 1:T) for (s in 1:S) h_eff[t,s] = alpha_s[s];
  } else {
    h_eff = h;
  }

  // ---- Production-price path. Latent OU (Level 2) when nested; set equal to the
  //      exogenous Xz in single-level mode (then it is unused by the likelihood,
  //      whose Level-1 branch reverts to theta_s instead of Phi). ----
  matrix[T,S] Phi;
  // Level-2 derived quantities exposed (length 0 when n_levels == 1) so the R
  // side can extract the production reversion speed kappa_p, the Level-2 mean
  // intercept mu_const, and the Level-2 innovation scale sigma_p. The hard/soft
  // separation diagnostic compares kappa_p against the base market speed kappa_s.
  vector[n_levels >= 2 ? S : 0] kappa_p;
  vector[n_levels >= 2 ? S : 0] kappa_m_base;   // market speed at zTMG = 0, in (0, kappa_cap)
  vector[n_levels >= 2 ? S : 0] mu_const;
  vector[n_levels >= 2 ? S : 0] sigma_p;        // baseline (constant) L2 innovation scale
  // Level-2 richness, exposed only when switched on (length 0 otherwise).
  vector<upper=0>[l2_cubic == 1 ? S : 0] a3_p;  // cubic restoring on Phi
  matrix[l2_sv == 1 ? T : 0, l2_sv == 1 ? S : 0] h_p;  // L2 SV log-variance
  if (n_levels >= 2) {
    for (s in 1:S) {
      kappa_m_base[s] = kappa_cap * inv_logit(kappa_tilde[s]);
      mu_const[s] = m0_0[1] + m0_COM[1] * COM_s[s] + sigma_m0[1] * m0_z[s];
      sigma_p[s]  = exp(sigma_p_tilde[s]);
      if (theta_sep == 1) {
        // Hard ordering: kappa_p in (0, kappa_m_base) by construction.
        kappa_p[s] = kappa_m_base[s] * inv_logit(sep0[1] + sigma_kappa_p[1] * kappa_p_z[s]);
      } else {
        // Soft: bounded free speed in (0, kappa_cap); the prior centers it below
        // the market speed (kappa_p0 < kappa center) without imposing the order.
        kappa_p[s] = kappa_cap * inv_logit(kappa_p0[1] + kappa_p_COM[1] * COM_s[s]
                                           + sigma_kappa_p[1] * kappa_p_z[s]);
      }
    }
    // L2 cubic coefficient (a3_p < 0), like Level 1.
    if (l2_cubic == 1) for (s in 1:S) a3_p[s] = -exp(a3_p_tilde[s]);
    // L2 stochastic volatility: stationary non-centered AR(1) log-variance.
    if (l2_sv == 1) {
      matrix[T,S] h_p_std;
      for (s in 1:S) {
        h_p_std[1,s] = h_p_raw[1,s] / sqrt(1 - square(rho_p[s]) + 1e-8);
        for (t in 2:T) h_p_std[t,s] = rho_p[s] * h_p_std[t-1,s] + h_p_raw[t,s];
        for (t in 1:T) h_p[t,s] = alpha_p[s] + sigma_eta_p[s] * h_p_std[t,s];
      }
    }
    // CENTERED: the latent path IS the parameter Phi_lat. The OU drift and the
    // (constant or SV) innovation scale enter the increment PRIOR in the model
    // block instead of a forward recursion over non-centered innovations, so the
    // scale * innovation ridge no longer exists. Phi (transformed) aliases
    // Phi_lat to keep the reduce_sum / generated-quantities / extractor interface
    // ([T,S] matrix named "Phi") bit-identical downstream.
    Phi = Phi_lat;
  } else {
    Phi = Xz;   // single-level: production price is exogenous (unused by Level-1 drift)
  }

  // ---- Effective anchor for the latent Phi (2-level only). In "meas" mode it is
  //      the fixed constructed index (data); in "recon" mode it is the
  //      reconstruction Phi = k + K G' with K lognormal-uncertain, standardized
  //      to the latent-Phi units. With sigma_K_recon -> 0 the reconstruction
  //      collapses onto the constructed index, so recon nests meas. ----
  matrix[n_levels >= 2 ? T : 0, n_levels >= 2 ? S : 0] Phi_anchor_eff;
  if (n_levels >= 2) {
    if (k_recon == 1) {
      for (s in 1:S) for (t in 1:T) {
        real K_unc   = K_hat[t, s] * exp(sigma_K_recon * z_K[t, s]);  // lognormal(K_hat, sigma_K)
        real phi_rec = k_cost[t, s] + K_unc * Gprime_raw[t];          // Phi = k + K G' (eq. 9)
        Phi_anchor_eff[t, s] = (phi_rec - phi_recon_center[s]) / phi_recon_scale[s];
      }
    } else {
      Phi_anchor_eff = Phi_anchor_z;
    }
  }
}

model {
  // ---- Level-1 hyperpriors (centers preserve the legacy implied locations) ----
  theta0    ~ normal(0, 1);
  theta_COM ~ normal(0, 0.5);
  sigma_theta ~ normal(0, 1);

  kappa0    ~ normal(-1, 0.5);
  kappa_COM ~ normal(0, 0.5);
  sigma_kappa ~ normal(0, 0.5);

  a3_0      ~ normal(log(0.05), 0.4);
  sigma_a3  ~ normal(0, 0.3);

  beta00    ~ normal(0, 0.5);
  beta0_COM ~ normal(0, 0.5);
  sigma_beta0 ~ normal(0, 0.5);

  theta_z ~ normal(0, 1);
  kappa_z ~ normal(0, 1);
  a3_z    ~ normal(0, 1);
  beta0_z ~ normal(0, 1);

  beta1 ~ normal(beta1_prior_mean, beta1_prior_sd);

  alpha_s      ~ normal(0, 1);
  rho_s        ~ normal(rho_prior_mean, rho_prior_sd);
  sigma_eta_s  ~ normal(0, 0.5);
  to_vector(h_raw) ~ normal(0, 1);
  nu_tilde     ~ gamma(nu_prior_shape, nu_prior_rate);

  if (soft_wedge == 1) delta_z ~ normal(zTMG_exo - zTMG_byK, sigma_delta_z);
  gamma ~ normal(0, 0.5);

  // ---- Level-2 priors and measurement (nested only) ----
  if (n_levels >= 2) {
    // Effective measurement SD: a fixed datum in the K-deterministic mode, the
    // sampled parameter otherwise (D-IMPL-9.4). Resolved with an if/else (not a
    // ternary) so the sigma_phi_meas[1] index is never touched in the fixed mode,
    // where that array has length 0.
    real sigma_meas_eff;
    if (sigma_phi_meas_fixed == 1) sigma_meas_eff = sigma_phi_meas_value;
    else                           sigma_meas_eff = sigma_phi_meas[1];

    kappa_p0[1]     ~ normal(-1.5, 0.5);    // slower than the market by default
    kappa_p_COM[1]  ~ normal(0, 0.5);
    sigma_kappa_p[1] ~ normal(0, 0.5);
    m0_0[1]         ~ normal(0, 1);
    m0_COM[1]       ~ normal(0, 0.5);
    sigma_m0[1]     ~ normal(0, 1);
    m1[1]           ~ normal(0, 0.5);
    if (n_levels == 3) m_v[1] ~ normal(0, 0.5);  // Level-3 value coupling (falsifiable)
    sep0[1]         ~ normal(-1, 1);        // inv_logit(sep0) < 0.5 -> kappa_p < kappa_m
    // Prior on the measurement SD only when it is estimated; in the fixed mode it
    // is a datum and carries no prior.
    if (sigma_phi_meas_fixed == 0)
      sigma_phi_meas[1] ~ normal(0, sigma_phi_meas_prior_sd);

    kappa_p_z ~ normal(0, 1);
    m0_z      ~ normal(0, 1);
    sigma_p_tilde ~ normal(-1, 0.5);
    // Anchor the start: the first latent state is centered on the first anchor.
    to_vector(Phi_lat[1, ]) ~ normal(to_vector(Phi_anchor_eff[1, ]), 1);

    // ---- Level-2 richness priors (each block is empty unless its switch is on) ----
    if (l2_cubic == 1)    a3_p_tilde ~ normal(log(0.05), 0.4);
    if (l2_sv == 1) {
      alpha_p     ~ normal(0, 1);
      rho_p       ~ normal(rho_prior_mean, rho_prior_sd);
      sigma_eta_p ~ normal(0, 0.5);
      to_vector(h_p_raw) ~ normal(0, 1);
    }
    if (l2_studentt == 1) nu_p_tilde[1] ~ gamma(nu_prior_shape, nu_prior_rate);

    // CENTERED process prior on the latent production price: each one-step
    // increment is an OU step toward the G'-driven mean, with the (constant or
    // SV) Level-2 scale and Student-t tails when l2_studentt. This is the CP
    // analogue of the former non-centered Phi_innov prior; because the data
    // identify the increment directly (via the anchor measurement below), the
    // scale sigma_p is identified by the residual variance and there is no
    // scale * innovation ridge (D-IMPL-9.1). The cubic argument and the SV scale
    // keep the same clamps as before (Level-2 stability bound, D-IMPL-9), now
    // protecting the drift mean rather than a forward recursion.
    for (t in 2:T) {
      for (s in 1:S) {
        // Level-2 reversion mean for the latent production price. In n_levels == 3
        // it gains the value coupling m_v * V (production gravitates around values,
        // Capital III ch. IX); V_anchor_z is length 0 otherwise so the index is
        // only touched when n_levels == 3.
        real mu_ts  = mu_const[s] + m1[1] * Gprime[t];
        if (n_levels == 3) mu_ts += m_v[1] * V_anchor_z[t, s];
        real dev_p  = Phi_lat[t-1,s] - mu_ts;
        real dev_pc = fmin(fmax(dev_p, -10), 10);
        real cub_p  = (l2_cubic == 1) ? a3_p[s] * dev_pc^3 : 0;
        real sp     = (l2_sv == 1) ? fmin(fmax(exp(0.5 * h_p[t,s]), 1e-8), 1e1)
                                   : sigma_p[s];
        real m_incr = Phi_lat[t-1,s] + kappa_p[s] * (-dev_p + cub_p);
        if (l2_studentt == 1)
          Phi_lat[t,s] ~ student_t(2 + nu_p_tilde[1], m_incr, sp);
        else
          Phi_lat[t,s] ~ normal(m_incr, sp);
      }
    }

    // f_Y: lognormal prior on the uncertain total capital K (non-centered).
    if (k_recon == 1) to_vector(z_K) ~ normal(0, 1);

    // Anchor: the effective production-price index (constructed in "meas",
    // reconstructed with uncertain K in "recon") is a noisy measurement of the
    // latent Phi, with the effective measurement SD (fixed datum or estimated).
    for (t in 1:T)
      Phi_anchor_eff[t, ] ~ normal(Phi[t, ], sigma_meas_eff);
  }

  // ---- Market-price likelihood (reduce_sum), restricted to t <= T_lik ----
  target += reduce_sum(ou_nested_partial_sum, t_idx, grainsize,
                       n_levels, T_lik, Yz, Xz,
                       Phi,
                       COM_ts, zTMG_byK, soft_wedge, delta_z,
                       com_wmean_train, com_wsd_train, mu_xz,
                       theta_s, kappa_s, kappa_tilde, kappa_cap,
                       a3_s, beta0_s, beta1,
                       h_eff, nu, com_in_mean, gamma,
                       l1_cubic, l1_studentt);
}

generated quantities {
  // Pointwise one-step-ahead log density on the fitted window (2..T_lik), zero
  // elsewhere so the R side subsets consistently for PSIS-LOO.
  matrix[T,S] log_lik;
  for (t in 1:T) for (s in 1:S) log_lik[t,s] = 0;

  for (t in 2:T_lik) {
    real ztmg_eff = zTMG_byK[t];
    if (soft_wedge == 1) ztmg_eff += delta_z[t];
    ztmg_eff = fmin(fmax(ztmg_eff, -1e6), 1e6);
    for (s in 1:S) {
      real denom_sd = com_wsd_train[s];
      denom_sd = (denom_sd > 1e-12) ? denom_sd : 1.0;
      real com_std = (COM_ts[t-1,s] - com_wmean_train[s]) / denom_sd;
      com_std = fmin(fmax(com_std, -1e6), 1e6);
      real com_term = (com_in_mean == 1) ? gamma * com_std : 0;
      real sd_   = fmin(fmax(exp(0.5 * h_eff[t,s]), 1e-8), 1e8);

      real mean_;
      if (n_levels == 1) {
        real zlag  = Yz[t-1,s] - theta_s[s];
        real drift = kappa_s[s] * (theta_s[s] - Yz[t-1,s] + a3_s[s] * zlag^3);
        real betaT = beta0_s[s] + beta1 * ztmg_eff;
        mean_ = drift + betaT * (Xz[t-1,s] - mu_xz[s]) + com_term;
      } else {
        real dev     = Yz[t-1,s] - Phi[t-1,s];
        real kappa_m = kappa_cap * inv_logit(kappa_tilde[s] + beta1 * ztmg_eff);
        real cub     = (l1_cubic == 1) ? a3_s[s] * dev^3 : 0;
        mean_ = kappa_m * (-dev + cub) + com_term;
      }

      real resid = Yz[t,s] - Yz[t-1,s] - mean_;
      if (n_levels >= 2 && l1_studentt == 0)
        log_lik[t,s] = normal_lpdf(resid | 0, sd_);
      else
        log_lik[t,s] = student_t_lpdf(resid | nu, 0, sd_);
    }
  }
}
