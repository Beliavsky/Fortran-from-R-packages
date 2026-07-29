# Example runner script for bayesianOU package
# This script demonstrates the full workflow

library(bayesianOU)

# =============================================================================
# DATA PREPARATION BLOCK
# =============================================================================

# Raw matrices (WITHOUT z-score - the function handles standardization)
Y <- as.matrix(BayesianCPI[, -1])
X <- as.matrix(MarxistPricesofProduction_Index[, -1])

# Basic safeguards
stopifnot(all(is.finite(Y)), all(is.finite(X)))
stopifnot(nrow(Y) == nrow(X), ncol(Y) == ncol(X))

Tn <- nrow(Y)
S <- ncol(Y)
T_train <- max(2L, floor(Tn * 0.70))

# Align COM and CAPITAL_TOTAL to Y columns
align_cols <- function(A, B) {
  if (!is.null(colnames(A)) && !is.null(colnames(B))) {
    missing_in_B <- setdiff(colnames(A), colnames(B))
    extra_in_B <- setdiff(colnames(B), colnames(A))
    if (length(missing_in_B) > 0) {
      stop("Missing in B: ", paste(missing_in_B, collapse = ", "))
    }
    if (length(extra_in_B) > 0) {
      message("Extra columns in B will be ignored: ", 
              paste(extra_in_B, collapse = ", "))
    }
    B[, colnames(A), drop = FALSE]
  } else {
    B
  }
}

COM <- align_cols(Y, as.matrix(COM))
CAPITAL_TOTAL <- align_cols(Y, as.matrix(CAPITAL_TOTAL))

stopifnot(
  nrow(COM) == Tn, ncol(COM) == S,
  nrow(CAPITAL_TOTAL) == Tn, ncol(CAPITAL_TOTAL) == S
)

# TMG: avoid near-zero variance
if (stats::sd(TMG$TMG, na.rm = TRUE) < 1e-8) {
  warning("TMG nearly constant; centering to avoid zero variance")
  TMG$TMG <- TMG$TMG - mean(TMG$TMG, na.rm = TRUE)
}

# =============================================================================
# MODEL FITTING
# =============================================================================

factor_from <- "X"

results_nl <- fit_ou_nonlinear_tmg(
  results_robust = list(),
  Y = Y,
  X = X,
  TMG = TMG$TMG,
  COM = COM,
  CAPITAL_TOTAL = CAPITAL_TOTAL,
  model = "base",
  priors = list(sigma_delta = 0.002),
  com_in_mean = TRUE,
  chains = 4,
  iter = 8000,
  warmup = 5000,
  thin = 1,
  cores = 4,
  threads_per_chain = 4,
  hard_sum_zero = TRUE,
  orthogonalize_tmg = TRUE,
  factor_from = factor_from,
  adapt_delta = 0.999,
  max_treedepth = 18,
  seed = 1234,
  init = 0.1,
  moment_match = TRUE,
  verbose = TRUE
)

# =============================================================================
# CONVERGENCE EVIDENCE
# =============================================================================

conv_evidence <- extract_convergence_evidence(results_nl)

# =============================================================================
# OUT-OF-SAMPLE EVALUATION (h=1 only)
# =============================================================================

# Reproduce internal standardization for manual OOS
zY <- zscore_train(Y, T_train)
zX <- zscore_train(X, T_train)

mu_tmg <- mean(TMG$TMG[seq_len(T_train)], na.rm = TRUE)
sd_tmg <- stats::sd(TMG$TMG[seq_len(T_train)], na.rm = TRUE)
if (!is.finite(sd_tmg) || sd_tmg < 1e-8) sd_tmg <- 1
zTMG <- (TMG$TMG - mu_tmg) / sd_tmg

# Compute common factor
compute_F <- function(Mz) {
  sv <- svd(scale(Mz, center = TRUE, scale = FALSE))
  f <- sv$u[, 1]
  as.numeric((f - mean(f[seq_len(T_train)])) / stats::sd(f[seq_len(T_train)]))
}

Ft <- if (factor_from == "X") compute_F(zX$Mz) else compute_F(zY$Mz)

# Orthogonalize TMG
fit_t <- stats::lm(zTMG[seq_len(T_train)] ~ Ft[seq_len(T_train)])
zTMG_use <- as.numeric(zTMG - cbind(1, Ft) %*% stats::coef(fit_t))

# Extract posterior summary
summ <- extract_posterior_summary(results_nl$factor_ou$stan_fit)

# OOS h=1
oos_h1 <- evaluate_oos(
  summ, zY$Mz, zX$Mz, zTMG_use, T_train,
  COM_ts = COM,
  K_ts = CAPITAL_TOTAL,
  com_in_mean = TRUE,
  horizons = c(1)
)
print(oos_h1)

# =============================================================================
# Y-X RELATIONSHIP SUMMARY
# =============================================================================

stopifnot(!is.null(results_nl$beta_tmg$beta_point))
B <- results_nl$beta_tmg$beta_point
Tn <- nrow(B)
S <- ncol(B)

X_component <- B[2:Tn, ] * zX$Mz[1:(Tn - 1), ]
dY_obs <- zY$Mz[2:Tn, ] - zY$Mz[1:(Tn - 1), ]
corr_global <- stats::cor(
  as.vector(dY_obs),
  as.vector(X_component),
  use = "complete.obs"
)

message(sprintf(
  "\nGlobal correlation dY vs X-induced component: %.3f",
  corr_global
))

beta_sector_mean <- colMeans(B, na.rm = TRUE)
print(sort(beta_sector_mean, decreasing = TRUE)[1:min(10, S)])

# =============================================================================
# STANDARD VALIDATION
# =============================================================================

validate_ou_fit(results_nl)

message(sprintf(
  "\nDivergences: %d | Rhat_max: %.4f | share(Rhat>1.01): %.4f",
  results_nl$diagnostics$divergences,
  results_nl$diagnostics$rhat_max,
  results_nl$diagnostics$rhat_share
))

print(results_nl$diagnostics$loo)

# =============================================================================
# MCMC DIAGNOSTICS
# =============================================================================

message("\n=== MCMC DIAGNOSTICS ===")
message(sprintf("Divergences: %d", results_nl$diagnostics$divergences))
message(sprintf("Rhat_max: %.4f", results_nl$diagnostics$rhat_max))
message(sprintf("Rhat > 1.01 (proportion): %.4f", results_nl$diagnostics$rhat_share))

# =============================================================================
# KEY OU PARAMETERS
# =============================================================================

message("\n=== OU PARAMETERS (medians) ===")
message("kappa_s (mean reversion speed) - summary:")
print(summary(summ$kappa_s))

message("\ntheta_s (equilibrium) - summary:")
print(summary(summ$theta_s))

message("\na3_s (cubic nonlinearity) - summary:")
print(summary(summ$a3_s))

message("\nrho_s (SV persistence) - summary:")
print(summary(summ$rho_s))

message(sprintf("\nnu (t degrees of freedom): %.2f", summ$nu))
message(sprintf("gamma (COM effect): %.4f", summ$gamma))
message(sprintf("beta1 (global TMG effect): %.4f", summ$beta1))

# =============================================================================
# HALF-LIFE BY SECTOR
# =============================================================================

message("\n=== HALF-LIFE OF REVERSION BY SECTOR (years) ===")
half_life <- log(2) / summ$kappa_s
names(half_life) <- colnames(Y)
print(round(sort(half_life), 2))

# =============================================================================
# SECTORAL CORRELATIONS
# =============================================================================

message("\n=== SECTORAL CORRELATIONS ===")
com_mean_sector <- colMeans(COM[seq_len(T_train), ], na.rm = TRUE)
message(sprintf(
  "Corr(beta0_s, mean_COM): %.4f",
  stats::cor(summ$beta0_s, com_mean_sector, use = "complete.obs")
))
message(sprintf(
  "Corr(kappa_s, mean_COM): %.4f",
  stats::cor(summ$kappa_s, com_mean_sector, use = "complete.obs")
))

# =============================================================================
# VARIANCE DECOMPOSITION
# =============================================================================

message("\n=== X COMPONENT CONTRIBUTION TO dY ===")
var_dY <- stats::var(as.vector(dY_obs), na.rm = TRUE)
var_X_comp <- stats::var(as.vector(X_component), na.rm = TRUE)
cov_dY_X <- stats::cov(
  as.vector(dY_obs),
  as.vector(X_component),
  use = "complete.obs"
)

message(sprintf("Var(dY): %.4f", var_dY))
message(sprintf("Var(X_component): %.4f", var_X_comp))
message(sprintf("Cov(dY, X_comp): %.4f", cov_dY_X))
message(sprintf("Implied R-squared: %.4f", cov_dY_X^2 / (var_dY * var_X_comp)))

# =============================================================================
# KEY HYPOTHESIS TESTS
# =============================================================================

message("\n=== HYPOTHESIS EVIDENCE ===")

fit <- results_nl$factor_ou$stan_fit
beta1_draws <- as.vector(fit$draws("beta1", format = "matrix"))

message(sprintf("P(beta1 > 0): %.4f", mean(beta1_draws > 0)))
message(sprintf(
  "95%% CI for beta1: [%.4f, %.4f]",
  stats::quantile(beta1_draws, 0.025),
  stats::quantile(beta1_draws, 0.975)
))

message(sprintf(
  "P(all kappa_s > 0): %.4f",
  conv_evidence$prob_convergence
))

message(sprintf(
  "Proportion of sectors with a3 < 0: %.4f",
  mean(summ$a3_s < 0)
))

# =============================================================================
# CUMULATIVE EFFECT BY HORIZON
# =============================================================================

kappa_median <- stats::median(summ$kappa_s)
half_life_median <- stats::median(log(2) / summ$kappa_s)
beta1_est <- summ$beta1

message("\n=== CONVERGENCE METRICS ===")
message(sprintf("beta1 (direct X effect): %.5f", beta1_est))
message(sprintf("Median kappa: %.4f", kappa_median))
message(sprintf("Median half-life (years): %.2f", half_life_median))

long_run_effect <- beta1_est / kappa_median
message(sprintf("\nLong-run effect (beta1/kappa): %.4f", long_run_effect))

horizons <- c(2, 5, 10, 15, 20)
message("\n=== CUMULATIVE EFFECT BY HORIZON ===")
for (Th in horizons) {
  effect_T <- beta1_est * (1 - exp(-kappa_median * Th)) / kappa_median
  pct_total <- (1 - exp(-kappa_median * Th)) * 100
  message(sprintf(
    "T = %2d years: Effect = %.4f (%.1f%% of long-run)",
    Th, effect_T, pct_total
  ))
}
