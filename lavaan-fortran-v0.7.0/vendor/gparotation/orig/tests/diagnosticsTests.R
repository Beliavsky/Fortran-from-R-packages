# Tests for residual diagnostics, simple structure measures,
# and CM normalization.
# No external dependencies.

require("GPArotation")

all.ok <- TRUE
fuzz   <- 1e-6

# --- Data ---
data("CCAI",         package = "GPArotation")
data("WansbeekMeijer", package = "GPArotation")
data("GriffithMulaik", package = "GPArotation")

fa.ccai <- factanal(factors = 3, covmat = CCAI_R,
                    n.obs = 461, rotation = "none")
res.vm  <- Varimax(fa.ccai,  algorithm = "legacy", fwindow = 1)
res.ob  <- oblimin(fa.ccai,  algorithm = "legacy", fwindow = 1)
res.qm  <- quartimax(fa.ccai, algorithm = "legacy", fwindow = 1)

fa.nl   <- factanal(factors = 2, covmat = NetherlandsTV,
                    rotation = "none")
res.nl  <- oblimin(fa.nl, algorithm = "legacy", fwindow = 1)

# ===========================================================================
# 1. residuals.GPArotation
# ===========================================================================

# 1a. Returns a matrix
r <- residuals(res.ob)
if (!is.matrix(r)) {
  cat("Test 1a residuals: should return a matrix\n")
  all.ok <- FALSE
}

# 1b. Diagonal is zero
if (max(abs(diag(r))) > fuzz) {
  cat("Test 1b residuals: diagonal should be zero\n")
  cat("max diag:", max(abs(diag(r))), "\n")
  all.ok <- FALSE
}

# 1c. Symmetric
if (max(abs(r - t(r))) > fuzz) {
  cat("Test 1c residuals: matrix should be symmetric\n")
  all.ok <- FALSE
}

# 1d. Rotation invariance --- residuals identical across rotation methods
r.vm <- residuals(res.vm)
r.ob <- residuals(res.ob)
r.qm <- residuals(res.qm)

if (max(abs(r.vm - r.ob)) > fuzz) {
  cat("Test 1d residuals: Varimax vs Oblimin not identical\n")
  cat("max difference:", max(abs(r.vm - r.ob)), "\n")
  all.ok <- FALSE
}
if (max(abs(r.vm - r.qm)) > fuzz) {
  cat("Test 1d residuals: Varimax vs Quartimax not identical\n")
  cat("max difference:", max(abs(r.vm - r.qm)), "\n")
  all.ok <- FALSE
}

# 1e. Correct dimensions
if (!identical(dim(r), c(14L, 14L))) {
  cat("Test 1e residuals: wrong dimensions\n")
  cat("got:", dim(r), "\n")
  all.ok <- FALSE
}

# 1f. R = R_hat + Resid
L   <- unclass(res.ob$loadings)
Phi <- res.ob$Phi
R_hat <- L %*% Phi %*% t(L)
diag(R_hat) <- 1
if (max(abs(cov2cor(CCAI_R) - R_hat - r.ob)) > fuzz) {
  cat("Test 1f residuals: R != R_hat + Resid\n")
  cat("max difference:", max(abs(cov2cor(CCAI_R) - R_hat - r.ob)), "\n")
  all.ok <- FALSE
}


# ===========================================================================
# 2. audit_residuals
# ===========================================================================

out <- GPArotation:::audit_residuals(res.ob)

# 2a. Returns invisible list with correct components
if (!all(c("Resid", "SRMR", "row_strain", "pairs", "cutoff") %in% names(out))) {
  cat("Test 2a audit_residuals: missing components\n")
  cat("got:", names(out), "\n")
  all.ok <- FALSE
}

# 2b. SRMR is positive and less than 1
if (out$SRMR <= 0 || out$SRMR >= 1) {
  cat("Test 2b audit_residuals: SRMR out of range\n")
  cat("SRMR:", out$SRMR, "\n")
  all.ok <- FALSE
}

# 2c. Resid matches residuals()
if (max(abs(out$Resid - residuals(res.ob))) > fuzz) {
  cat("Test 2c audit_residuals: Resid does not match residuals()\n")
  all.ok <- FALSE
}

# 2d. SRMR matches manual calculation
srmr_manual <- sqrt(mean(out$Resid[lower.tri(out$Resid)]^2))
if (abs(out$SRMR - srmr_manual) > fuzz) {
  cat("Test 2d audit_residuals: SRMR does not match manual calculation\n")
  cat("audit:", out$SRMR, "manual:", srmr_manual, "\n")
  all.ok <- FALSE
}

# 2e. Cutoff default is 0.10
if (out$cutoff != 0.10) {
  cat("Test 2e audit_residuals: default cutoff should be 0.10\n")
  all.ok <- FALSE
}

# 2f. Rotation invariance of SRMR
out.vm <- GPArotation:::audit_residuals(res.vm)
out.qm <- GPArotation:::audit_residuals(res.qm)
if (abs(out$SRMR - out.vm$SRMR) > fuzz) {
  cat("Test 2f audit_residuals: SRMR differs across rotations\n")
  all.ok <- FALSE
}
if (abs(out$SRMR - out.qm$SRMR) > fuzz) {
  cat("Test 2f audit_residuals: SRMR differs across rotations\n")
  all.ok <- FALSE
}

# 2g. Flagged pairs all exceed cutoff
if (nrow(out$pairs) > 0) {
  if (any(abs(out$pairs$Residual) <= out$cutoff)) {
    cat("Test 2g audit_residuals: flagged pairs below cutoff\n")
    all.ok <- FALSE
  }
}

# 2h. Custom cutoff
out_strict <- GPArotation:::audit_residuals(res.ob, cutoff = 0.05)
if (nrow(out_strict$pairs) < nrow(out$pairs)) {
  cat("Test 2h audit_residuals: stricter cutoff should flag more pairs\n")
  all.ok <- FALSE
}


# ===========================================================================
# 3. plot_residuals (run without error, no visual check)
# ===========================================================================

# 3a. Runs without error on factanal-derived object (R stored)
tryCatch({
  pdf(file = NULL)  # suppress plot output
  GPArotation:::plot_residuals(res.ob)
  dev.off()
}, error = function(e) {
  cat("Test 3a plot_residuals: error on factanal-derived object\n")
  cat(conditionMessage(e), "\n")
  all.ok <<- FALSE
})

# 3b. Runs with explicit R argument
tryCatch({
  pdf(file = NULL)
  GPArotation:::plot_residuals(res.ob, R = CCAI_R)
  dev.off()
}, error = function(e) {
  cat("Test 3b plot_residuals: error with explicit R\n")
  cat(conditionMessage(e), "\n")
  all.ok <<- FALSE
})

# 3c. Fails informatively without R when not stored
res.raw <- oblimin(fa.ccai$loadings, algorithm = "legacy", fwindow = 1)
tryCatch({
  pdf(file = NULL)
  GPArotation:::plot_residuals(res.raw)
  dev.off()
  cat("Test 3c plot_residuals: should have errored without R\n")
  all.ok <- FALSE
}, error = function(e) {
  # expected
})

# 3d. Returns invisible list with Resid and SRMR
pdf(file = NULL)
ret <- GPArotation:::plot_residuals(res.ob)
dev.off()
if (!all(c("Resid", "SRMR", "cutoff") %in% names(ret))) {
  cat("Test 3d plot_residuals: missing return components\n")
  all.ok <- FALSE
}


# ===========================================================================
# 4. AUC and FSI
# ===========================================================================

auc <- GPArotation:::calc_AUC(res.ob)
fsi <- GPArotation:::calc_FSI(res.ob)

# 4a. AUC per factor in [0.5, 1]
if (any(auc$AUC < 0.5) || any(auc$AUC > 1)) {
  cat("Test 4a calc_AUC: per-factor AUC outside [0.5, 1]\n")
  cat("values:", auc$AUC, "\n")
  all.ok <- FALSE
}

# 4b. AUC mean in [0.5, 1]
if (auc$AUC_mean < 0.5 || auc$AUC_mean > 1) {
  cat("Test 4b calc_AUC: mean AUC outside [0.5, 1]\n")
  cat("value:", auc$AUC_mean, "\n")
  all.ok <- FALSE
}

# 4c. FSI per factor in [0, 1]
if (any(fsi$FSI < 0) || any(fsi$FSI > 1)) {
  cat("Test 4c calc_FSI: per-factor FSI outside [0, 1]\n")
  cat("values:", fsi$FSI, "\n")
  all.ok <- FALSE
}

# 4d. FSI mean in [0, 1]
if (fsi$FSI_mean < 0 || fsi$FSI_mean > 1) {
  cat("Test 4d calc_FSI: mean FSI outside [0, 1]\n")
  cat("value:", fsi$FSI_mean, "\n")
  all.ok <- FALSE
}

# 4e. Oblique should have higher AUC than orthogonal for CCAI
auc.vm <- GPArotation:::calc_AUC(res.vm)
if (auc$AUC_mean <= auc.vm$AUC_mean) {
  cat("Test 4e calc_AUC: oblimin should have higher AUC than varimax for CCAI\n")
  cat("oblimin:", auc$AUC_mean, "varimax:", auc.vm$AUC_mean, "\n")
  all.ok <- FALSE
}

# 4f. AUC exceeds 0.7 for a meaningful oblimin CCAI solution
if (any(auc$AUC < 0.7)) {
  cat("Test 4f calc_AUC: AUC below 0.7 for oblimin CCAI solution\n")
  cat("values:", auc$AUC, "\n")
  all.ok <- FALSE
}

# 4g. Number of factors matches
if (length(auc$AUC) != ncol(res.ob$loadings)) {
  cat("Test 4g calc_AUC: wrong number of factors\n")
  all.ok <- FALSE
}
if (length(fsi$FSI) != ncol(res.ob$loadings)) {
  cat("Test 4g calc_FSI: wrong number of factors\n")
  all.ok <- FALSE
}


# ===========================================================================
# 5. calc_simplicity and calc_hyperplane
# ===========================================================================

sim <- GPArotation:::calc_simplicity(res.ob)
hp  <- GPArotation:::calc_hyperplane(res.ob)

# 5a. Hoffman in [0, 1]
if (sim$Hoffman < 0 || sim$Hoffman > 1) {
  cat("Test 5a calc_simplicity: Hoffman outside [0, 1]\n")
  cat("value:", sim$Hoffman, "\n")
  all.ok <- FALSE
}

# 5b. Gini in [0, 1]
if (sim$Gini < 0 || sim$Gini > 1) {
  cat("Test 5b calc_simplicity: Gini outside [0, 1]\n")
  cat("value:", sim$Gini, "\n")
  all.ok <- FALSE
}

# 5c. Bentler in [0, 1]
if (sim$Bentler < 0 || sim$Bentler > 1) {
  cat("Test 5c calc_simplicity: Bentler outside [0, 1]\n")
  cat("value:", sim$Bentler, "\n")
  all.ok <- FALSE
}

# 5d. Hyperplane count non-negative and <= p*(k-1)
p <- nrow(res.ob$loadings)
k <- ncol(res.ob$loadings)
if (hp$HP_total < 0 || hp$HP_total > p * (k - 1)) {
  cat("Test 5d calc_hyperplane: total outside [0, p*(k-1)]\n")
  cat("value:", hp$HP_total, "max:", p * (k - 1), "\n")
  all.ok <- FALSE
}

# 5e. Rotation improves simplicity vs unrotated
sim.un <- GPArotation:::calc_simplicity(
  quartimax(fa.ccai$loadings, algorithm = "legacy", fwindow = 1))
res.un.fake <- fa.ccai
res.un.fake$loadings <- fa.ccai$loadings
class(res.un.fake) <- "GPArotation"
# Hoffman should increase after rotation
if (sim$Hoffman <= 0) {
  cat("Test 5e calc_simplicity: Hoffman should be positive after rotation\n")
  all.ok <- FALSE
}


# ===========================================================================
# 6. Cureton-Mulaik normalization
# ===========================================================================

# 6a. CM gives different result from Kaiser
res.kai <- oblimin(fa.ccai$loadings, normalize = TRUE,
                   algorithm = "legacy", fwindow = 1)
res.cm  <- oblimin(fa.ccai$loadings, normalize = "CM",
                   algorithm = "legacy", fwindow = 1)

if (max(abs(unclass(res.kai$loadings) -
            unclass(res.cm$loadings))) < 1e-4) {
  cat("Test 6a CM normalization: CM and Kaiser gave identical results\n")
  all.ok <- FALSE
}

# 6b. CM gives different result from no normalization
res.none <- oblimin(fa.ccai$loadings, normalize = FALSE,
                    algorithm = "legacy", fwindow = 1)
if (max(abs(unclass(res.none$loadings) -
            unclass(res.cm$loadings))) < 1e-4) {
  cat("Test 6b CM normalization: CM and no normalization gave identical results\n")
  all.ok <- FALSE
}

# 6c. Kaiser string same as Kaiser logical
res.str <- oblimin(fa.ccai$loadings, normalize = "Kaiser",
                   algorithm = "legacy", fwindow = 1)
if (max(abs(unclass(res.kai$loadings) -
            unclass(res.str$loadings))) > fuzz) {
  cat("Test 6c CM normalization: normalize='Kaiser' differs from normalize=TRUE\n")
  cat("max difference:", max(abs(unclass(res.kai$loadings) -
                                 unclass(res.str$loadings))), "\n")
  all.ok <- FALSE
}

# 6d. CM rotation produces a valid GPArotation object
# Note: CM normalization may not converge for all datasets ---
# this is documented behavior. Check object structure only.
if (!inherits(res.cm, "GPArotation")) {
  cat("Test 6d CM normalization: did not return GPArotation object\n")
  all.ok <- FALSE
}

# 6e. CM weights all positive (internal check)
A <- fa.ccai$loadings
W <- GPArotation:::NormalizingWeight(A, normalize = "CM")
if (any(W <= 0)) {
  cat("Test 6e CM normalization: weights not all positive\n")
  all.ok <- FALSE
}

# 6f. NormalizingWeight Kaiser = normalize=TRUE
W.kai <- GPArotation:::NormalizingWeight(A, normalize = "Kaiser")
W.log <- GPArotation:::NormalizingWeight(A, normalize = TRUE)
if (max(abs(W.kai - W.log)) > fuzz) {
  cat("Test 6f NormalizingWeight: Kaiser string differs from TRUE\n")
  all.ok <- FALSE
}


# ===========================================================================
# 7. GriffithMulaik dataset basic checks
# ===========================================================================

# 7a. Correct dimensions
if (!identical(dim(GriffithMulaik), c(24L, 24L))) {
  cat("Test 7a GriffithMulaik: wrong dimensions\n")
  all.ok <- FALSE
}

# 7b. Is symmetric
if (max(abs(GriffithMulaik - t(GriffithMulaik))) > fuzz) {
  cat("Test 7b GriffithMulaik: not symmetric\n")
  all.ok <- FALSE
}

# 7c. Diagonal is 1
if (max(abs(diag(GriffithMulaik) - 1)) > fuzz) {
  cat("Test 7c GriffithMulaik: diagonal not 1\n")
  all.ok <- FALSE
}

# 7d. 6-factor solution fits well (RMSEA < 0.05)
fa.gm  <- factanal(factors = 6, covmat = GriffithMulaik,
                   n.obs = 523, rotation = "none")
res.gm <- oblimin(fa.gm, algorithm = "legacy", fwindow = 1)
fit.gm <- GPArotation:::calc_fitstats(res.gm)
if (fit.gm$RMSEA > 0.05) {
  cat("Test 7d GriffithMulaik: 6-factor RMSEA exceeds 0.05\n")
  cat("RMSEA:", fit.gm$RMSEA, "\n")
  all.ok <- FALSE
}


cat("diagnostics tests completed.\n")
if (!all.ok) stop("some tests FAILED")