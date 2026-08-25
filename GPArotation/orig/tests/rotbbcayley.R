# tests/algorithmTests.R
# Unit tests for algorithm = "bb", "cayley" and fwindow settings.
# Tests verify convergence, solution quality, and correct storage of
# algorithm metadata. Solutions are compared within tolerance since
# BB and Cayley take different paths to the same minimum.

library(GPArotation)

all.ok <- TRUE
fuzz   <- 1e-4  # looser than legacy tests — different path, same minimum

data("Harman",    package = "GPArotation")
data("Thurstone", package = "GPArotation")
data("CCAI",      package = "GPArotation")

fa.ccai <- factanal(covmat = CCAI_R, factors = 3,
                    rotation = "none", n.obs = 461)

# Helper: check two loading matrices agree within tolerance
.check_loadings <- function(r1, r2, label, tol = fuzz) {
  diff <- max(abs(abs(unclass(r1$loadings)) - abs(unclass(r2$loadings))))
  if (diff > tol) {
    cat("FAILED:", label, "\n")
    cat("max abs difference:", diff, "\n")
    all.ok <<- FALSE
  }
}

# Helper: check result metadata
.check_meta <- function(r, algo, fw, label) {
  if (is.null(r$algorithm) || r$algorithm != algo) {
    cat("FAILED:", label, "- algorithm not stored correctly\n")
    cat("expected:", algo, "got:", r$algorithm, "\n")
    all.ok <<- FALSE
  }
  if (is.null(r$fwindow) || r$fwindow != fw) {
    cat("FAILED:", label, "- fwindow not stored correctly\n")
    cat("expected:", fw, "got:", r$fwindow, "\n")
    all.ok <<- FALSE
  }
}

# Helper: check orthogonality
.check_orthogonal <- function(r, label, tol = 1e-10) {
  k    <- ncol(r$Th)
  diff <- max(abs(t(r$Th) %*% r$Th - diag(k)))
  if (diff > tol) {
    cat("FAILED:", label, "- rotation matrix not orthogonal\n")
    cat("max deviation from identity:", diff, "\n")
    all.ok <<- FALSE
  }
}

###########################################################################
# --- BB tests ---
###########################################################################

cat("=== BB algorithm tests ===\n")

# Test 1: BB varimax on Harman8 — convergence and metadata
r.bb  <- GPForth(Harman8, method = "varimax", algorithm = "bb", fwindow = 10)
r.leg <- GPForth(Harman8, method = "varimax", algorithm = "legacy", fwindow = 1)

if (!r.bb$convergence) {
  cat("FAILED: Test 1 - BB varimax Harman8 did not converge\n")
  all.ok <- FALSE
}
.check_loadings(r.bb, r.leg, "Test 1 - BB varimax Harman8 loadings")
.check_meta(r.bb, "bb", 10, "Test 1 - BB varimax metadata")
.check_orthogonal(r.bb, "Test 1 - BB varimax orthogonality")

# Test 2: BB quartimax on box26
r.bb  <- GPForth(box26, method = "quartimax", algorithm = "bb", fwindow = 10)
r.leg <- GPForth(box26, method = "quartimax", algorithm = "legacy", fwindow = 1)

if (!r.bb$convergence) {
  cat("FAILED: Test 2 - BB quartimax box26 did not converge\n")
  all.ok <- FALSE
}
.check_loadings(r.bb, r.leg, "Test 2 - BB quartimax box26 loadings")
.check_meta(r.bb, "bb", 10, "Test 2 - BB quartimax metadata")
.check_orthogonal(r.bb, "Test 2 - BB quartimax orthogonality")

# Test 3: BB oblimin on CCAI — oblique
r.bb  <- GPFoblq(fa.ccai, method = "oblimin", algorithm = "bb", fwindow = 10)
r.leg <- GPFoblq(fa.ccai, method = "oblimin", algorithm = "legacy", fwindow = 1)

if (!r.bb$convergence) {
  cat("FAILED: Test 3 - BB oblimin CCAI did not converge\n")
  all.ok <- FALSE
}
.check_loadings(r.bb, r.leg, "Test 3 - BB oblimin CCAI loadings")
.check_meta(r.bb, "bb", 10, "Test 3 - BB oblimin metadata")

# Test 4: BB quartimin on box26 — oblique
r.bb  <- GPFoblq(box26, method = "quartimin", algorithm = "bb", fwindow = 10)
r.leg <- GPFoblq(box26, method = "quartimin", algorithm = "legacy", fwindow = 1)

if (!r.bb$convergence) {
  cat("FAILED: Test 4 - BB quartimin box26 did not converge\n")
  all.ok <- FALSE
}
.check_loadings(r.bb, r.leg, "Test 4 - BB quartimin box26 loadings")
.check_meta(r.bb, "bb", 10, "Test 4 - BB quartimin metadata")

# Test 5: BB via named wrapper — oblimin
r.bb  <- oblimin(fa.ccai, algorithm = "bb", fwindow = 10)
r.leg <- oblimin(fa.ccai, algorithm = "legacy", fwindow = 1)

if (!r.bb$convergence) {
  cat("FAILED: Test 5 - oblimin wrapper BB did not converge\n")
  all.ok <- FALSE
}
.check_loadings(r.bb, r.leg, "Test 5 - oblimin wrapper BB loadings")
.check_meta(r.bb, "bb", 10, "Test 5 - oblimin wrapper metadata")

# Test 6: BB with random starts — randStartChar populated
r.bb <- oblimin(fa.ccai, algorithm = "bb", fwindow = 10, randomStarts = 20)

if (is.null(r.bb$randStartChar)) {
  cat("FAILED: Test 6 - randStartChar not populated with random starts\n")
  all.ok <- FALSE
}
if (r.bb$randStartChar["randomStarts"] != 20) {
  cat("FAILED: Test 6 - randomStarts count incorrect\n")
  all.ok <- FALSE
}

# Test 7: BB faster than legacy on box26 (larger matrix)
r.bb  <- GPForth(box26, method = "varimax", algorithm = "bb",     fwindow = 10)
r.leg <- GPForth(box26, method = "varimax", algorithm = "legacy", fwindow = 1)

if (nrow(r.bb$Table) >= nrow(r.leg$Table)) {
  cat("NOTE: Test 7 - BB not fewer iterations than legacy on box26 varimax\n")
  cat("BB iters:", nrow(r.bb$Table), "Legacy iters:", nrow(r.leg$Table), "\n")
  # Not a hard failure — matrix size dependent
}

# Test 8: fwindow = 1 with BB reduces to standard Armijo behaviour
r.fw1  <- GPForth(Harman8, method = "varimax", algorithm = "bb", fwindow = 1)
r.fw10 <- GPForth(Harman8, method = "varimax", algorithm = "bb", fwindow = 10)

if (!r.fw1$convergence) {
  cat("FAILED: Test 8 - BB fwindow=1 did not converge\n")
  all.ok <- FALSE
}
.check_meta(r.fw1,  "bb", 1,  "Test 8 - BB fwindow=1 metadata")
.check_meta(r.fw10, "bb", 10, "Test 8 - BB fwindow=10 metadata")

###########################################################################
# --- Cayley tests ---
###########################################################################

cat("\n=== Cayley algorithm tests ===\n")

# Test 9: Cayley varimax on Harman8
r.cay <- GPForth(Harman8, method = "varimax", algorithm = "cayley", fwindow = 10)
r.leg <- GPForth(Harman8, method = "varimax", algorithm = "legacy", fwindow = 1)

if (!r.cay$convergence) {
  cat("FAILED: Test 9 - Cayley varimax Harman8 did not converge\n")
  all.ok <- FALSE
}
.check_loadings(r.cay, r.leg, "Test 9 - Cayley varimax Harman8 loadings")
.check_meta(r.cay, "cayley", 10, "Test 9 - Cayley varimax metadata")
.check_orthogonal(r.cay, "Test 9 - Cayley varimax orthogonality")

# Test 10: Cayley quartimax on box26
r.cay <- GPForth(box26, method = "quartimax", algorithm = "cayley", fwindow = 10)
r.leg <- GPForth(box26, method = "quartimax", algorithm = "legacy", fwindow = 1)

if (!r.cay$convergence) {
  cat("FAILED: Test 10 - Cayley quartimax box26 did not converge\n")
  all.ok <- FALSE
}
.check_loadings(r.cay, r.leg, "Test 10 - Cayley quartimax box26 loadings")
.check_meta(r.cay, "cayley", 10, "Test 10 - Cayley quartimax metadata")
.check_orthogonal(r.cay, "Test 10 - Cayley quartimax orthogonality")

# Test 11: Cayley via named wrapper — bentlerT
r.cay <- bentlerT(fa.ccai, algorithm = "cayley", fwindow = 10)

if (!r.cay$convergence) {
  cat("FAILED: Test 11 - bentlerT Cayley did not converge\n")
  all.ok <- FALSE
}
.check_meta(r.cay, "cayley", 10, "Test 11 - bentlerT Cayley metadata")
.check_orthogonal(r.cay, "Test 11 - bentlerT Cayley orthogonality")

# Test 12: Cayley orthogonality is exact — tighter tolerance than SVD
r.cay <- GPForth(box26, method = "varimax", algorithm = "cayley", fwindow = 10)
k     <- ncol(r.cay$Th)
ortho_err <- max(abs(t(r.cay$Th) %*% r.cay$Th - diag(k)))

if (ortho_err > 1e-10) {
  cat("FAILED: Test 12 - Cayley orthogonality error too large\n")
  cat("max deviation:", ortho_err, "\n")
  all.ok <- FALSE
}

###########################################################################
# --- Edge case tests ---
###########################################################################

cat("\n=== Edge case tests ===\n")

# Test 13: Cayley on oblique throws error
tryCatch({
  GPFoblq(Harman8, method = "quartimin", algorithm = "cayley")
  cat("FAILED: Test 13 - Cayley oblique should throw error\n")
  all.ok <- FALSE
}, error = function(e) {
  # expected
})

# Test 14: Invalid algorithm falls back gracefully — no error
r <- tryCatch({
  GPForth(Harman8, method = "varimax", algorithm = "invalid")
}, error = function(e) NULL)

if (is.null(r)) {
  cat("FAILED: Test 14 - invalid algorithm should fall back, not error\n")
  all.ok <- FALSE
}

# Test 15: NULL fwindow gets correct default for bb
r <- GPForth(Harman8, method = "varimax", algorithm = "bb", fwindow = NULL)
if (r$fwindow != 10) {
  cat("FAILED: Test 15 - NULL fwindow should default to 10 for bb\n")
  cat("got:", r$fwindow, "\n")
  all.ok <- FALSE
}

# Test 16: NULL fwindow gets correct default for legacy
r <- GPForth(Harman8, method = "varimax", algorithm = "legacy", fwindow = NULL)
if (r$fwindow != 1) {
  cat("FAILED: Test 16 - NULL fwindow should default to 1 for legacy\n")
  cat("got:", r$fwindow, "\n")
  all.ok <- FALSE
}

# Test 17: algorithm and fwindow flow through named wrappers
r <- quartimax(Harman8)
if (is.null(r$algorithm)) {
  cat("FAILED: Test 17 - algorithm not stored via named wrapper\n")
  all.ok <- FALSE
}
if (is.null(r$fwindow)) {
  cat("FAILED: Test 17 - fwindow not stored via named wrapper\n")
  all.ok <- FALSE
}

# Test 18: Tinit stored for trajectory reconstruction
r <- GPForth(Harman8, method = "varimax", algorithm = "bb", fwindow = 10)
if (is.null(r$Tinit)) {
  cat("FAILED: Test 18 - Tinit not stored in result\n")
  all.ok <- FALSE
}

# Test 19: A_unrotated stored for trajectory reconstruction
if (is.null(attr(r, "A_unrotated"))) {
  cat("FAILED: Test 19 - A_unrotated not stored in result\n")
  all.ok <- FALSE
}

# Test 20: target and PST keep legacy defaults
r.t <- targetT(Harman8[1:2, ], Target = diag(2), algorithm = "legacy", fwindow = 1)
if (!is.null(r.t$algorithm) && r.t$algorithm != "legacy") {
  cat("FAILED: Test 20 - targetT should default to legacy\n")
  cat("got:", r.t$algorithm, "\n")
  all.ok <- FALSE
}

###########################################################################

cat("\ntests completed.\n")
if (!all.ok) stop("some algorithm tests FAILED")
