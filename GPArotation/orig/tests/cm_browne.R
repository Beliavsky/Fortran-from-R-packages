# tests/cm_browne.R
# Regression tests for Cureton-Mulaik (CM) Normalization
# Based on Browne (2001), Multivariate Behavioral Research, Tables 5 & 6.

library(GPArotation)

all.ok <- TRUE
fuzz   <- 1e-3

# --- 1. Browne (2001) Table 5 Constructed Matrix ---
Lambda <- matrix(c(
  0.90, 0.00, 0.00,
  0.90, 0.00, 0.00,
  0.90, 0.00, 0.00,
  0.90, 0.00, 0.00,
  0.00, 0.80, 0.00,
  0.00, 0.80, 0.00,
  0.00, 0.80, 0.00,
  0.00, 0.00, 0.70,
  0.00, 0.00, 0.70,
  0.00, 0.00, 0.70,
  0.10, 0.10, 0.00,
  0.10, 0.00, 0.10
), nrow = 12, ncol = 3, byrow = TRUE)

# Expected CM weights w_i: Var1-10 receive floor weight (0.001); Var11-12 receive ~0.925
h <- sqrt(rowSums(Lambda^2))
exp_w <- c(rep(0.001, 10), 0.924951, 0.924951)
exp_al <- h / exp_w

# Verify NormalizingWeight output
al_calc <- GPArotation:::NormalizingWeight(Lambda, normalize = "CM")

if (max(abs(al_calc[, 1] - exp_al)) > fuzz) {
  cat("FAILED: CM NormalizingWeight calculation does not match Browne Table 5\n")
  all.ok <- FALSE
}

# --- 2. Table 6 Oblique Rotation Recovery Test (rho = 0.5) ---
rho <- 0.5
Phi_true <- matrix(rho, 3, 3)
diag(Phi_true) <- 1.0

# Implied correlation matrix and unrotated factor extraction
R_common <- Lambda %*% Phi_true %*% t(Lambda)
ev <- eigen(R_common)
A_unrotated <- ev$vectors[, 1:3] %*% diag(sqrt(ev$values[1:3]))

# Oblique CF-Varimax rotation with CM normalization
res_cm <- cfQ(A_unrotated, normalize = "CM", kappa = 1/12)
# GPFoblq(A_unrotated, method = "cf-varimax", normalize = "CM")

if (!res_cm$convergence) {
  cat("FAILED: Oblique CM rotation did not converge for rho = 0.5\n")
  all.ok <- FALSE
}

# Verify that primary factor loadings recover the simple structure for Var1-4
L_canon <- GPArotation:::.sortGPALoadings(res_cm)$loadings
if (max(abs(L_canon[1:4, 1])) < 0.85) {
  cat("FAILED: CM normalization failed to recover primary factor loadings at rho = 0.5\n")
  all.ok <- FALSE
}

cat("Cureton-Mulaik Browne (2001) tests completed.\n")
if (!all.ok) stop("some CM Browne tests FAILED")
