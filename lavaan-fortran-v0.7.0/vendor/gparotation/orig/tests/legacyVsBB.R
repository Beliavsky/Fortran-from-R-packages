# tests/legacyVsBB.R
# Regression tests: BB and Cayley produce same solutions as legacy
# within numerical tolerance. Different optimization paths, same minimum.

library(GPArotation)

all.ok <- TRUE
fuzz   <- 1e-3  # tolerance — different paths, same minimum

data(Harman,    package = "GPArotation")
data(Thurstone, package = "GPArotation")
data(CCAI,      package = "GPArotation")

fa.ccai <- factanal(covmat = CCAI_R, factors = 3,
                    rotation = "none", n.obs = 461)

# Helper: Compare canonicalized loadings and Phi (handles sign flips and column swaps)
.check <- function(r.bb, r.leg, label, tol = fuzz) {
    # Sort both objects canonically (aligns loadings, Phi, and Th consistently)
    bb_canon  <- GPArotation:::.sortGPALoadings(r.bb)
    leg_canon <- GPArotation:::.sortGPALoadings(r.leg)
    
    # 1. Compare Loadings (using absolute values for sign invariance)
    L.bb  <- unclass(bb_canon$loadings)
    L.leg <- unclass(leg_canon$loadings)
    diff_L <- max(abs(abs(L.bb) - abs(L.leg)))
    
    if (diff_L > tol) {
        cat("FAILED Loadings:", label, "\n")
        cat("max abs difference:", diff_L, "\n")
        all.ok <<- FALSE
    }
    
    # 2. Compare Phi matrix for oblique rotations
    if (!is.null(bb_canon$Phi) && !is.null(leg_canon$Phi)) {
        diff_phi <- max(abs(abs(bb_canon$Phi) - abs(leg_canon$Phi)))
        if (diff_phi > tol) {
            cat("FAILED Phi:", label, "\n")
            cat("max Phi difference:", diff_phi, "\n")
            all.ok <<- FALSE
        }
    }
}

cat("=== BB vs Legacy ===\n")

# --- Orthogonal ---
# Test 1: varimax Harman8
r.bb  <- GPForth(Harman8, method = "varimax", algorithm = "bb",     fwindow = 10)
r.leg <- GPForth(Harman8, method = "varimax", algorithm = "legacy", fwindow = 1)
.check(r.bb, r.leg, "Test 1 - varimax Harman8")
if (!r.bb$convergence)  { cat("FAILED: Test 1 - BB did not converge\n");  all.ok <- FALSE }
if (!r.leg$convergence) { cat("FAILED: Test 1 - legacy did not converge\n"); all.ok <- FALSE }

# Test 2: quartimax box26
r.bb  <- GPForth(box26, method = "quartimax", algorithm = "bb",     fwindow = 10)
r.leg <- GPForth(box26, method = "quartimax", algorithm = "legacy", fwindow = 1)
.check(r.bb, r.leg, "Test 2 - quartimax box26")
if (!r.bb$convergence)  { cat("FAILED: Test 2 - BB did not converge\n");  all.ok <- FALSE }
if (!r.leg$convergence) { cat("FAILED: Test 2 - legacy did not converge\n"); all.ok <- FALSE }

# Test 3: geominT CCAI
r.bb  <- GPForth(fa.ccai, method = "geomin", algorithm = "bb",     fwindow = 10)
r.leg <- GPForth(fa.ccai, method = "geomin", algorithm = "legacy", fwindow = 1)
.check(r.bb, r.leg, "Test 3 - geominT CCAI")
if (!r.bb$convergence)  { cat("FAILED: Test 3 - BB did not converge\n");  all.ok <- FALSE }

# Test 4: entropy box26
r.bb  <- GPForth(box26, method = "entropy", algorithm = "bb",     fwindow = 10)
r.leg <- GPForth(box26, method = "entropy", algorithm = "legacy", fwindow = 1)
.check(r.bb, r.leg, "Test 4 - entropy box26")
if (!r.bb$convergence)  { cat("FAILED: Test 4 - BB did not converge\n");  all.ok <- FALSE }

# Test 5: tandemII box26
r.bb  <- GPForth(box26, method = "tandemII", algorithm = "bb",     fwindow = 10)
r.leg <- GPForth(box26, method = "tandemII", algorithm = "legacy", fwindow = 1)
.check(r.bb, r.leg, "Test 5 - tandemII box26")
if (!r.bb$convergence)  { cat("FAILED: Test 5 - BB did not converge\n");  all.ok <- FALSE }

cat("\n=== Cayley vs Legacy ===\n")

# Test 6: varimax Harman8
r.cay <- GPForth(Harman8, method = "varimax", algorithm = "cayley", fwindow = 10)
r.leg <- GPForth(Harman8, method = "varimax", algorithm = "legacy", fwindow = 1)
.check(r.cay, r.leg, "Test 6 - cayley varimax Harman8")
if (!r.cay$convergence) { cat("FAILED: Test 6 - Cayley did not converge\n"); all.ok <- FALSE }

# Test 7: quartimax box26
r.cay <- GPForth(box26, method = "quartimax", algorithm = "cayley", fwindow = 10)
r.leg <- GPForth(box26, method = "quartimax", algorithm = "legacy", fwindow = 1)
.check(r.cay, r.leg, "Test 7 - cayley quartimax box26")
if (!r.cay$convergence) { cat("FAILED: Test 7 - Cayley did not converge\n"); all.ok <- FALSE }

# Test 8: bentlerT CCAI
r.cay <- GPForth(fa.ccai, method = "bentler", algorithm = "cayley", fwindow = 10)
r.leg <- GPForth(fa.ccai, method = "bentler", algorithm = "legacy", fwindow = 1)
.check(r.cay, r.leg, "Test 8 - cayley bentlerT CCAI")
if (!r.cay$convergence) { cat("FAILED: Test 8 - Cayley did not converge\n"); all.ok <- FALSE }

cat("\n=== BB vs Legacy Oblique ===\n")

# Test 9: quartimin Harman8
r.bb  <- GPFoblq(Harman8, method = "quartimin", algorithm = "bb",     fwindow = 10)
r.leg <- GPFoblq(Harman8, method = "quartimin", algorithm = "legacy", fwindow = 1)
.check(r.bb, r.leg, "Test 9 - quartimin Harman8")
if (!r.bb$convergence)  { cat("FAILED: Test 9 - BB did not converge\n");  all.ok <- FALSE }
if (!r.leg$convergence) { cat("FAILED: Test 9 - legacy did not converge\n"); all.ok <- FALSE }

# Test 10: oblimin box26
r.bb  <- GPFoblq(box26, method = "oblimin", algorithm = "bb",     fwindow = 10)
r.leg <- GPFoblq(box26, method = "oblimin", algorithm = "legacy", fwindow = 1)
.check(r.bb, r.leg, "Test 10 - oblimin box26")
if (!r.bb$convergence)  { cat("FAILED: Test 10 - BB did not converge\n");  all.ok <- FALSE }

# Test 11: geomin CCAI
r.bb  <- GPFoblq(fa.ccai, method = "geomin", algorithm = "bb",     fwindow = 10)
r.leg <- GPFoblq(fa.ccai, method = "geomin", algorithm = "legacy", fwindow = 1)
.check(r.bb, r.leg, "Test 11 - geominQ CCAI")
if (!r.bb$convergence)  { cat("FAILED: Test 11 - BB did not converge\n");  all.ok <- FALSE }

# Test 12: simplimax CCAI
r.bb  <- GPFoblq(fa.ccai, method = "simplimax", algorithm = "bb",     fwindow = 10)
r.leg <- GPFoblq(fa.ccai, method = "simplimax", maxit = 5000, algorithm = "legacy", fwindow = 1)
.check(r.bb, r.leg, "Test 12 - simplimax CCAI")
if (!r.bb$convergence)  { cat("FAILED: Test 12 - BB did not converge\n");  all.ok <- FALSE }

# Test 13: bentlerQ CCAI
r.bb  <- GPFoblq(fa.ccai, method = "bentler", algorithm = "bb",     fwindow = 10)
r.leg <- GPFoblq(fa.ccai, method = "bentler", algorithm = "legacy", fwindow = 1)
.check(r.bb, r.leg, "Test 13 - bentlerQ CCAI")
if (!r.bb$convergence)  { cat("FAILED: Test 13 - BB did not converge\n");  all.ok <- FALSE }

cat("\n=== Via Named Wrappers ===\n")

# Test 14: oblimin wrapper default (bb, fwindow=10)
r.bb  <- oblimin(fa.ccai)
r.leg <- oblimin(fa.ccai, algorithm = "legacy", fwindow = 1)
.check(r.bb, r.leg, "Test 14 - oblimin wrapper default")

# Test 15: Varimax wrapper default (bb, fwindow=10)
r.bb  <- Varimax(fa.ccai)
r.leg <- Varimax(fa.ccai, algorithm = "legacy", fwindow = 1)
.check(r.bb, r.leg, "Test 15 - Varimax wrapper default")

# Test 16: geominQ wrapper default (bb, fwindow=10)
r.bb  <- geominQ(fa.ccai)
r.leg <- geominQ(fa.ccai, algorithm = "legacy", fwindow = 1)
.check(r.bb, r.leg, "Test 16 - geominQ wrapper default")

cat("\ntests completed.\n")
if (!all.ok) stop("some legacyVsBB tests FAILED")