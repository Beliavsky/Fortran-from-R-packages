# bayesianOU 0.2.0

Unified single-/two-level engine for the nested gravitation model (market price
reverting to a latent production price, which has its own G'-driven OU). Work in
progress toward 0.2.0.

## New

* **Block 9 validation (internal + external + simulation-based calibration).**
  A full validation block on the canonical K-deterministic variant, documented in
  the `nested-gravitation` vignette (sections 16-17): recovery + 90% coverage
  (`n=3`, pooled 0.91); adversarial negative controls (a unit-root DGP yields no
  false fast gravitation; an `m_v=0` DGP is covered at 0); multiple-imputation
  coverage with a quantified, conservative disaggregation bias on `kappa_m`;
  simulation-based calibration (Talts et al. 2018) showing uniform posterior ranks
  for the linear core; and external comparisons (held-out forecast vs random walk /
  AR(1) / `kappa=0`; PSIS-LOO across depth; wedge placebos) reported anti-overreach.
  Section 17 develops the logical verdict: the "negative" external results are
  deductively entailed by slow gravitation and corroborate, rather than refute, the
  structural thesis. No package code changed (validation scripts are external).
* **Unified external multiple imputation of the model inputs (`fit_ou_nested_mi`).**
  Beyond the market price `phi_draws`, the anchors `X` (production price `Phi`),
  `V_value` (value), `CAPITAL_TOTAL` (`K`) and `COM` may now each vary by
  imputation: each accepts either a fixed matrix (the legacy path, used verbatim
  for every imputation, byte-identical results) or a per-imputation set -- a list
  of `D` matrices or a `[T, S, D]` array -- paired one-to-one with `phi_draws`.
  This lets a single external generator propagate ALL the construction
  uncertainty (market disaggregation, the imputed `v`/`p` split that feeds `K` and
  `V`, and the reconstructed `Phi`) through the convergent K-DET / Variant-A
  engine, pooled by Rubin -- so the uncertainty of `V` and `K` is carried by
  EXTERNAL multiple imputation rather than an in-model latent (the in-model
  Variant B / `recon` modes are not reopened). The observed aggregates `TMG` and
  `Gprime` stay fixed. Per-imputation lengths are validated up-front (a mismatch
  fails before any Stan fit). Additive, R-side only: the Stan model is untouched
  and `n_levels` in {1,2} and the all-fixed path are byte-identical (equivalence
  and the 7 goldens re-verified, 0 regressions). Decision D-IMPL-13.5 / D-IMPL-14.1.
* **Multiple imputation extended to Level 3 (`fit_ou_nested_mi`, `n_levels = 3`).**
  The MI driver now accepts the value model: a named `V_value` argument (the
  direct-price anchor `c+v+p`, fixed across imputations like `X`), the value
  coupling `m_v` added to the default Rubin-pooled set, and the latent-`Phi`
  pooling and the `kappa^m > kappa^p` separation evidence widened from the
  2-level case to `n_levels >= 2` (mirroring the Stan `== 2 -> >= 2` refactor of
  D-IMPL-10.1). `n_levels` in {1,2} is unchanged. At `n_levels = 3` the joint
  separation typically collapses (kappa^p is re-identified once the value term
  absorbs the production mean). Additive, R-side only. Decision D-IMPL-11.1.
* **Faithful out-of-sample recursion for Level 3 (`evaluate_oos_nested`).** The
  forecast recursion now adds the value term to the Level-2 attractor,
  `mu_{s,t} = m0_s + m1 G'_t + m_v V_{s,t}`, via a new optional `V_z` argument and
  an `m_v` component in `meds`. Both gated: the term is inert (byte-identical to
  the Level-2 forecast) unless the standardized anchor `V_z` and a finite `m_v`
  are both supplied, so `n_levels <= 2` is unaffected. Decision D-IMPL-11.2.
* **Per-imputation convergence over the pooled parameters in the MI driver.**
  `fit_ou_nested_mi` now records, for each imputation, the chain-aware
  `rhat_max_pooled`, `ess_bulk_min_pooled` and `ess_tail_min_pooled` over the
  Rubin-pooled parameters -- the binding diagnostic for the reliability of the
  pooling -- plus the global `rhat_max`, `rhat_share`, `ess_bulk_min` (over all
  parameters, latent states included; reference only) and `divergences`. This
  makes per-imputation convergence verifiable to an impeccable standard
  (R-hat < 1.01, ESS-bulk and ESS-tail >= 1000, 0 divergences) on the
  quantities that actually enter the pooled estimates, rather than on a global
  extreme that a single weakly-identified latent cell can dominate. Additive
  and numerically inert.
* **Level 3 (values): `n_levels = 3`.** The nested cascade extends to a third
  level: the latent production price reverts toward a mean that also tracks the
  VALUE anchor `V` (direct prices `c+v+p` = `DirectPrices_Index`, constructed
  directly from `k_cost + EBO`), so prices of production gravitate around values
  (Marx, Capital III ch. IX). New argument `V_value` (T x S, required for
  `n_levels = 3`); the production reversion mean gains the term `m_v * V` with the
  coupling `m_v` estimated under a neutral, falsifiable prior (the ch. IX
  hypothesis is `m_v > 0`, adjudicated by the data). Values enter as a datum
  (direct empirical construction), never solved simultaneously (no Leontief
  inverse). The Level-2 gates were widened from `== 2` to `>= 2` so `n_levels = 3`
  inherits all of the 2-level machinery; `n_levels` in {1,2} is byte-for-byte
  unchanged (equivalence and all goldens re-certified, plus a new `level3`
  golden). Decision D-IMPL-10.1.
* **Fixed measurement SD in the K-deterministic mode (`sigma_phi_meas_fixed`).**
  In the deterministic limit the estimated anchor measurement SD piled up against
  its lower bound 0 (a boundary funnel, G6; `rhat` 2.208 on the real panel). The
  SD can now be held FIXED as a datum instead of estimated (`sigma_phi_meas_fixed
  = 0.05`), which is the `sigma_meas -> 0` limit done well and removes the funnel.
  On the real 37-sector panel this collapses `rhat_max` 2.208 -> 1.079 and lifts
  E-BFMI from ~0.02 to ~0.9-1.1 (the boundary funnel was also dragging the energy
  mixing). The estimated mode (looser prior) remains for the K-stochastic / recon
  regime. Decision D-IMPL-9.4.
* **Centered Level-2 latent (`Phi_lat`).** The latent production price is now
  sampled directly (centered parametrization) instead of being reconstructed
  from non-centered innovations. On the real 37-sector panel the former NCP
  factorization induced a multiplicative non-identification ridge
  (`cor(sigma_p_tilde, log||Phi_innov||^2) ~ -0.98`, E-BFMI ~0.03, rhat 1.5-2.0
  in `Phi_innov`/`sigma_p_tilde`/`m0_z`); the centered form removes it (those
  parameters now converge, the fit is ~13x faster and no longer saturates the
  tree depth). The single-level mode is byte-for-byte unchanged (equivalence and
  goldens re-certified). Decision D-IMPL-9.1.
* **Geometry-adaptive sampling engine** (`ou_geom_target()`,
  `ou_geom_metric_euclidean()`, `ou_geom_metric_riemannian()`, `ou_geom_hmc()`,
  `ou_geom_bridge()`) — a self-contained R-native HMC with a pluggable metric
  (dense Euclidean or position-dependent Riemannian SoftAbs / supplied-Fisher)
  and an E-BFMI diagnostic, for probing and remedying hard posterior geometry.
  Opt-in; the cmdstan fit path is untouched. Decision D-IMPL-9.2.
* **`fit_ou_nested()`** — single engine. `n_levels = 1` reproduces the legacy
  single-level model; `n_levels = 2` fits the nested model with the latent
  production price `Phi`, anchored by the constructed production-price index as a
  noisy measurement. Dispatch over `theta_separation` (`"soft"`/`"hard"`
  time-scale separation κ^m > κ^p), `k_uncertainty` (`"meas"` σ_Φ measurement
  error; `"recon"` reconstruction with uncertain K, eq. 11), `kappa_cap`
  (stability cap) and `level_spec` (per-level richness; all named configurations
  realised).
* **`fit_ou_nested_mi()`** — multiple-imputation driver (Rubin's rule + mixture
  posterior) coupling the disaggregation draws of the market price to the OU
  model; default `M = 25`.
* The 2-level fit returns the latent `Phi` trajectory with credible bands, the
  Level-2 parameters (`kappa_p`, `kappa_m_base`, `mu_const`, `sigma_p`, `m1`,
  `sigma_phi_meas`) and the posterior evidence for κ^m > κ^p.
* **2-level diagnostics layer.** `evaluate_oos_nested()` — out-of-sample
  recursion in which the market price reverts to the *latent* production price
  `Phi` propagated forward through its own Level-2 OU (never to a realized `X`);
  `diagnostics$oos` is now populated for 2-level fits. `extract_mu_trajectory()`
  reconstructs the `G'`-driven Level-2 mean `mu_t = m0_s + m1 G'_t` with bands
  (`fit$mu_path`). The 2-level result now carries the S3 class
  `ou_nested_2level` at the top level, with `print`/`summary`/`plot` methods
  (`plot` types `"phi"`, `"mu"`, `"separation"`); `ou_nested_mi` gains
  `print`/`summary`/`plot` (`"phi"`, `"fmi"`).
* **Per-level richness switches are now realised** (`level_spec`). Each level
  toggles `cubic` drift, `sv` (stochastic volatility), `student_t` tails and
  `hierarchy` (cross-sector partial pooling); the switches act only in the
  2-level mode (single-level is always the full legacy Level 1). The new
  `ou_level_spec()` helper builds the named experiment configurations
  `"canonical"` (N1-full/N2-lean), `"both_full"`, `"both_lean"` and `"n1_lean"`;
  the single-level arm is `n_levels = 1`. Level-2 enrichment (cubic / SV /
  Student-t on the latent production price) carries a Level-2 stability bound —
  the cubic argument and the SV scale are clamped so the *propagated* latent path
  cannot diverge (the Level-2 analogue of the 0.1.4 reversion-speed bound) — and
  the 2-level fit uses a reduced default init radius so the latent recursion
  starts in a well-behaved region. The canonical configuration reproduces the
  exact linear Gaussian Level-2 OU (the richness blocks are length 0 when off).
* **`k_uncertainty = "recon"` is now wired** (eq. 11 / Hallazgo 2). The latent-
  `Phi` anchor is reconstructed as `Phi = k + K G'` with the total capital
  advanced `K` lognormally uncertain (non-centered prior `f_Y` around the point
  estimate `CAPITAL_TOTAL`, scale `priors$sigma_K_recon`, default 0.10). The
  reconstruction is standardized so that the `sigma_K -> 0` limit reproduces the
  `"meas"` anchor exactly (recon nests meas); a consistency guard warns when the
  constructed index `X` departs from `k + K G'`. Requires the new `k_cost`
  argument (cost price `c + v`, raw units) and `n_levels = 2`.
* **Level-2 richness exposed in the summaries** (`summary`/`plot`). When the
  Level-2 switches are on, the latent production price's cubic coefficient
  `a3_p`, Student-t degrees of freedom `nu_p` and time-varying SV scale path
  `sigma_p_t = exp(h_p/2)` are added to `factor_ou$level2`, the Level-2 summary
  table and a new `plot(type = "sv_p")`. They are derived from the posterior
  draws, so the canonical configuration and the bit-exact single-level path are
  unchanged (the fields are `NULL` when the switch is off).
* **Validation layer for the experiment.** A per-configuration parameter-recovery
  test (`tests/testthat/test-recovery-by-config.R`, gated by
  `BAYESOU_RUN_RECOVERY_CONFIG`) checking 90% credible-interval coverage and LOO
  discrimination, with the heavier script `validacion/recovery_by_config.R`; and
  per-configuration bit-for-bit `log_lik` golden guards
  (`tests/testthat/test-golden-configs.R`, fixtures generated by
  `validacion/make_golden_configs.R`).

## Changed

* Single source of truth is now **`inst/stan/ou_nested.stan`** (covers both
  modes). `fit_ou_nonlinear_tmg()` is a thin wrapper of `fit_ou_nested(n_levels =
  1)`; its behaviour and result structure are unchanged. `ou_nested_stan_code()`
  replaces `ou_nonlinear_tmg_stan_code()`.
* **2-level stability fix**: the Euler reversion speeds κ^m, κ^p are bounded to
  `(0, kappa_cap)` via a logit link (TMG modulation inside the link). Unbounded
  speeds let the latent path explore the explosive region (κ > 2 ⇒ |1−κ| > 1),
  which broke convergence; the single-level path is untouched (`exp` link).

## Removed

* `inst/stan/ou_nonlinear_tmg.stan` (the legacy single-level file) and
  `ou_nonlinear_tmg_stan_code()`. The single-level mode of `ou_nested.stan` was
  certified bit-for-bit equivalent to the retired file (recomputing `log_lik`
  from the same draws via `generate_quantities`); the equivalence is frozen as a
  golden fixture and guarded by `test-equivalence.R`.

# bayesianOU 0.1.4

Robustness overhaul. Several changes alter results and are intentional; review
the methodology section of the README and re-run with your data.

## Correctness fixes

* **PSIS-LOO**: the log-likelihood was passed to `loo::loo()` as a 3-D array
  `[draws, time, sector]`, which `loo` interprets as `[iterations, chains,
  observations]` — it silently treated time as chains and only the sectors as
  observations. It is now reshaped to a proper `[draws x observations]` matrix
  over the fitted window, with a `chain_id` so `relative_eff()` is correct.
* **Train/test leakage**: the likelihood was summed over the full sample even
  when a train/test split was used, so the "out-of-sample" evaluation was
  contaminated. New `fit_window` argument (`"train"` default / `"full"`) keeps
  the likelihood, `log_lik`, and OOS evaluation coherent.
* **Real hierarchy**: `kappa0/kappa_COM/sigma_kappa`, `a3_0/sigma_a3`,
  `beta00/beta0_COM/sigma_beta0` were declared but unused — only `theta_s` was
  hierarchical. All four sector blocks now use a non-centered hierarchical
  (partial-pooling) parameterization that actually uses the hyperparameters.
* **`evaluate_oos` index bug**: for horizons longer than the test window the
  `(T_train+1):(Tn-hh+1)` sequence ran backwards and indexed out of range; it
  now guards the window and returns `NA`/`n_obs = 0` instead.
* **Factor leakage**: `use_train_loadings` now defaults to `TRUE` (loadings from
  training only), avoiding look-ahead in the orthogonalized TMG regressor.
* **`compare_models_loo`**: validates that both `loo` objects exist, are of
  class `loo`, and have matching observation counts; `deltaELPD` is now computed
  unambiguously as `elpd_new - elpd_base`.

## Priors (defaults changed — override via `priors`)

* `beta1 ~ Normal(0, 0.5)` (was `Normal(0.5, 0.25)`): neutral on the TMG-effect
  hypothesis (no sign baked into the prior).
* `nu_tilde ~ Gamma(2, 0.1)` (was `Exponential(3)`): weakly informative; the old
  prior forced `nu` into `(2,3)` (extreme heavy tails).
* `rho_s ~ Normal(0.7, 0.2)` (was `Normal(0.90, 0.05)`): less rigid SV
  persistence.
* New tunable prior entries: `beta1_mean`, `beta1_sd`, `nu_shape`, `nu_rate`,
  `rho_mean`, `rho_sd`.

## Diagnostics and robustness

* `validate_ou_fit()` is now a real validator: structural checks plus a Pareto-k
  summary with a warning when PSIS-LOO is unreliable; separates MCMC convergence
  from dynamic mean reversion.
* `extract_convergence_evidence()` no longer prints "CONVERGENCE GUARANTEED"; it
  reports dynamic mean-reversion evidence and is aliased by the clearer
  `kappa_stability_evidence()`.
* `count_divergences()` / LOO are surfaced with `loo_pareto_k` in `diagnostics`.

## Engineering

* The TMG wedge `delta_z` now has length 0 in the default hard case
  (`hard_sum_zero = TRUE`). Previously it was sampled as `T` parameters pinned by
  a `normal(0, 1e-6)` prior, whose near-zero scale throttled the NUTS step size
  and made sampling crawl. Removing them is both faster and better-mixing.
* The Stan model is now compiled into a writable user cache
  (`tools::R_user_dir`), not the (possibly read-only) package directory, and the
  compiled binary is no longer left in the source tree.
* Single source of truth for the Stan model: `inst/stan/ou_nonlinear_tmg.stan`
  is canonical and is read/compiled directly; the R string copy is gone.
* `fit_ou_nonlinear_tmg()` no longer calls `setwd()`; validates inputs
  (finiteness, `iter > warmup`, `train_frac`); caps threads to avoid
  oversubscription; exposes `train_frac`; `thin` default is now `1`; no longer
  duplicates the full draws array in the result (access via `stan_fit`).
* Declared `parallel` in Imports; removed stray `LazyData: true` (no `data/`);
  fixed repository URL capitalization.
* New tests: `evaluate_oos` index guard, PSIS-LOO reshape, and a Stan-based
  parameter-recovery test (skipped when no backend is available).
