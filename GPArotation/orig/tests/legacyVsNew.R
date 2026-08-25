# Regression tests: GPForth and GPFoblq against legacy implementations.
# Tests that the 2027 updates produce results equivalent to the
# original 2008 implementations within numerical tolerance.
# Uses tolerance-based comparison (1e-5) for loadings and Th,
# identical() only for convergence indicators and Table structure.

require("GPArotation")

all.ok <- TRUE
tol    <- 1e-5  # numerical tolerance for loadings/Th comparisons

# Helper: strip all attributes except dim for bare matrix comparison
.strip <- function(x) {
  x <- unclass(x)
  attributes(x) <- list(dim = dim(x))
  x
}

# Test matrices
data(Harman,   package = "GPArotation")
data(Thurstone, package = "GPArotation")

A2 <- Harman8  # 8 x 2
A3 <- box26    # 26 x 3

# --- Test 1: GPForth vs GPForth.legacy, varimax, identity start ---
r1  <- GPForth(A2, method = "varimax", algorithm = "legacy", fwindow = 1)
r1L <- GPArotation:::GPForth.legacy(A2, method = "varimax")

if (max(abs(.strip(r1$loadings) - .strip(r1L$loadings))) > tol) {
  cat("Test 1 failed: GPForth varimax loadings not equivalent to legacy\n")
  cat("max difference:", max(abs(r1$loadings - r1L$loadings)), "\n")
  all.ok <- FALSE
}
if (max(abs(r1$Th - r1L$Th)) > tol) {
  cat("Test 1 failed: GPForth varimax Th not equivalent to legacy\n")
  cat("max difference:", max(abs(r1$Th - r1L$Th)), "\n")
  all.ok <- FALSE
}
if (!identical(r1$convergence, r1L$convergence)) {
  cat("Test 1 failed: GPForth varimax convergence not identical to legacy\n")
  all.ok <- FALSE
}

# --- Test 2: GPForth vs GPForth.legacy, quartimax, 3 factors ---
r2  <- GPForth(A3, method = "quartimax", algorithm = "legacy", fwindow = 1)
r2L <- GPArotation:::GPForth.legacy(A3, method = "quartimax")

if (max(abs(.strip(r2$loadings) - .strip(r2L$loadings))) > tol) {
  cat("Test 2 failed: GPForth quartimax loadings not equivalent to legacy\n")
  cat("max difference:", max(abs(r2$loadings - r2L$loadings)), "\n")
  all.ok <- FALSE
}
if (max(abs(r2$Th - r2L$Th)) > tol) {
  cat("Test 2 failed: GPForth quartimax Th not equivalent to legacy\n")
  cat("max difference:", max(abs(r2$Th - r2L$Th)), "\n")
  all.ok <- FALSE
}
if (!identical(r2$convergence, r2L$convergence)) {
  cat("Test 2 failed: GPForth quartimax convergence not identical to legacy\n")
  all.ok <- FALSE
}

# --- Test 3: GPForth vs GPForth.legacy, random start ---
set.seed(42)
Tmat2 <- Random.Start(2)
r3  <- GPForth(A2, Tmat = Tmat2, method = "varimax", algorithm = "legacy",
               fwindow = 1)
set.seed(42)
Tmat2 <- Random.Start(2)
r3L <- GPArotation:::GPForth.legacy(A2, Tmat = Tmat2, method = "varimax")

if (max(abs(.strip(r3$loadings) - .strip(r3L$loadings))) > tol) {
  cat("Test 3 failed: GPForth random start loadings not equivalent to legacy\n")
  cat("max difference:", max(abs(r3$loadings - r3L$loadings)), "\n")
  all.ok <- FALSE
}
if (max(abs(r3$Th - r3L$Th)) > tol) {
  cat("Test 3 failed: GPForth random start Th not equivalent to legacy\n")
  cat("max difference:", max(abs(r3$Th - r3L$Th)), "\n")
  all.ok <- FALSE
}

# --- Test 4: GPForth vs GPForth.legacy, with normalization ---
r4  <- GPForth(A2, method = "varimax", normalize = TRUE,
               algorithm = "legacy", fwindow = 1)
r4L <- GPArotation:::GPForth.legacy(A2, method = "varimax", normalize = TRUE)

if (max(abs(.strip(r4$loadings) - .strip(r4L$loadings))) > tol) {
  cat("Test 4 failed: GPForth normalized loadings not equivalent to legacy\n")
  cat("max difference:", max(abs(r4$loadings - r4L$loadings)), "\n")
  all.ok <- FALSE
}
if (max(abs(r4$Th - r4L$Th)) > tol) {
  cat("Test 4 failed: GPForth normalized Th not equivalent to legacy\n")
  cat("max difference:", max(abs(r4$Th - r4L$Th)), "\n")
  all.ok <- FALSE
}

# --- Test 5: GPForth vs GPForth.legacy, methodArgs ---
r5  <- GPForth(A2, method = "cf", methodArgs = list(kappa = 0.3),
               algorithm = "legacy", fwindow = 1)
r5L <- GPArotation:::GPForth.legacy(A2, method = "cf",
                                     methodArgs = list(kappa = 0.3))

if (max(abs(.strip(r5$loadings) - .strip(r5L$loadings))) > tol) {
  cat("Test 5 failed: GPForth cf kappa=0.3 loadings not equivalent to legacy\n")
  cat("max difference:", max(abs(r5$loadings - r5L$loadings)), "\n")
  all.ok <- FALSE
}
if (max(abs(r5$Th - r5L$Th)) > tol) {
  cat("Test 5 failed: GPForth cf kappa=0.3 Th not equivalent to legacy\n")
  cat("max difference:", max(abs(r5$Th - r5L$Th)), "\n")
  all.ok <- FALSE
}

# --- Test 6: GPFoblq vs GPFoblq.legacy, quartimin, identity start ---
r6  <- GPFoblq(A2, method = "quartimin", algorithm = "legacy", fwindow = 1)
r6L <- GPArotation:::GPFoblq.legacy(A2, method = "quartimin")

if (max(abs(.strip(r6$loadings) - .strip(r6L$loadings))) > tol) {
  cat("Test 6 failed: GPFoblq quartimin loadings not equivalent to legacy\n")
  cat("max difference:", max(abs(r6$loadings - r6L$loadings)), "\n")
  all.ok <- FALSE
}
if (max(abs(r6$Phi - r6L$Phi)) > tol) {
  cat("Test 6 failed: GPFoblq quartimin Phi not equivalent to legacy\n")
  cat("max difference:", max(abs(r6$Phi - r6L$Phi)), "\n")
  all.ok <- FALSE
}
if (max(abs(r6$Th - r6L$Th)) > tol) {
  cat("Test 6 failed: GPFoblq quartimin Th not equivalent to legacy\n")
  cat("max difference:", max(abs(r6$Th - r6L$Th)), "\n")
  all.ok <- FALSE
}
if (!identical(r6$convergence, r6L$convergence)) {
  cat("Test 6 failed: GPFoblq quartimin convergence not identical to legacy\n")
  all.ok <- FALSE
}

# --- Test 7: GPFoblq vs GPFoblq.legacy, oblimin, 3 factors ---
r7  <- GPFoblq(A3, method = "oblimin", algorithm = "legacy", fwindow = 1)
r7L <- GPArotation:::GPFoblq.legacy(A3, method = "oblimin")

if (max(abs(.strip(r7$loadings) - .strip(r7L$loadings))) > tol) {
  cat("Test 7 failed: GPFoblq oblimin loadings not equivalent to legacy\n")
  cat("max difference:", max(abs(r7$loadings - r7L$loadings)), "\n")
  all.ok <- FALSE
}
if (max(abs(r7$Phi - r7L$Phi)) > tol) {
  cat("Test 7 failed: GPFoblq oblimin Phi not equivalent to legacy\n")
  cat("max difference:", max(abs(r7$Phi - r7L$Phi)), "\n")
  all.ok <- FALSE
}
if (!identical(r7$convergence, r7L$convergence)) {
  cat("Test 7 failed: GPFoblq oblimin convergence not identical to legacy\n")
  all.ok <- FALSE
}

# --- Test 8: GPFoblq vs GPFoblq.legacy, random start ---
set.seed(42)
Tmat3 <- Random.Start(3)
r8  <- GPFoblq(A3, Tmat = Tmat3, method = "quartimin",
               algorithm = "legacy", fwindow = 1)
set.seed(42)
Tmat3 <- Random.Start(3)
r8L <- GPArotation:::GPFoblq.legacy(A3, Tmat = Tmat3, method = "quartimin")

if (max(abs(.strip(r8$loadings) - .strip(r8L$loadings))) > tol) {
  cat("Test 8 failed: GPFoblq random start loadings not equivalent to legacy\n")
  cat("max difference:", max(abs(r8$loadings - r8L$loadings)), "\n")
  all.ok <- FALSE
}
if (max(abs(r8$Phi - r8L$Phi)) > tol) {
  cat("Test 8 failed: GPFoblq random start Phi not equivalent to legacy\n")
  cat("max difference:", max(abs(r8$Phi - r8L$Phi)), "\n")
  all.ok <- FALSE
}
if (max(abs(r8$Th - r8L$Th)) > tol) {
  cat("Test 8 failed: GPFoblq random start Th not equivalent to legacy\n")
  cat("max difference:", max(abs(r8$Th - r8L$Th)), "\n")
  all.ok <- FALSE
}

# --- Test 9: GPFoblq vs GPFoblq.legacy, with normalization ---
r9  <- GPFoblq(A2, method = "quartimin", normalize = TRUE,
               algorithm = "legacy", fwindow = 1)
r9L <- GPArotation:::GPFoblq.legacy(A2, method = "quartimin",
                                     normalize = TRUE)

if (max(abs(.strip(r9$loadings) - .strip(r9L$loadings))) > tol) {
  cat("Test 9 failed: GPFoblq normalized loadings not equivalent to legacy\n")
  cat("max difference:", max(abs(r9$loadings - r9L$loadings)), "\n")
  all.ok <- FALSE
}
if (max(abs(r9$Phi - r9L$Phi)) > tol) {
  cat("Test 9 failed: GPFoblq normalized Phi not equivalent to legacy\n")
  cat("max difference:", max(abs(r9$Phi - r9L$Phi)), "\n")
  all.ok <- FALSE
}

# --- Test 10: convergence indicators agree ---
if (!identical(r1$convergence, r1L$convergence)) {
  cat("Test 10 failed: GPForth convergence indicator not identical to legacy\n")
  all.ok <- FALSE
}
if (!identical(r6$convergence, r6L$convergence)) {
  cat("Test 10 failed: GPFoblq convergence indicator not identical to legacy\n")
  all.ok <- FALSE
}

# --- Test 11: Table structure for non-converged case ---
r11  <- GPForth(A3, method = "simplimax", maxit = 3,
                algorithm = "legacy", fwindow = 1)
r11L <- GPArotation:::GPForth.legacy(A3, method = "simplimax", maxit = 3)

if (!identical(nrow(r11$Table), nrow(r11L$Table))) {
  cat("Test 11 failed: GPForth Table nrow not identical to legacy\n")
  cat("new:", nrow(r11$Table), "legacy:", nrow(r11L$Table), "\n")
  all.ok <- FALSE
}
if (abs(unname(r11$Table[nrow(r11$Table), 2]) -
        unname(r11L$Table[nrow(r11L$Table), 2])) > tol) {
  cat("Test 11 failed: GPForth Table last row f value not equivalent to legacy\n")
  all.ok <- FALSE
}

# --- Test 12: Table structure for non-converged oblique case ---
r12  <- GPFoblq(A3, method = "simplimax", maxit = 3,
                algorithm = "legacy", fwindow = 1)
r12L <- GPArotation:::GPFoblq.legacy(A3, method = "simplimax", maxit = 3)

if (abs(unname(r12$Table[nrow(r12$Table), 2]) -
        unname(r12L$Table[nrow(r12L$Table), 2])) > tol) {
  cat("Test 12 failed: GPFoblq Table last row f value not equivalent to legacy\n")
  cat("new:", nrow(r12$Table), "legacy:", nrow(r12L$Table), "\n")
  all.ok <- FALSE
}

cat("Legacy regression tests completed.\n")
if (!all.ok) stop("some tests FAILED")