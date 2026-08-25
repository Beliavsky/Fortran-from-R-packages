# tests/normalizationTests.R
# Tests for NormalizingWeight normalization options

library(GPArotation)

all.ok <- TRUE
data("Harman", package = "GPArotation")

# Helper
.check <- function(condition, label) {
  if (!condition) {
    cat("FAILED:", label, "\n")
    all.ok <<- FALSE
  }
}

# --- Test 1: normalize = FALSE returns all ones ---
W <- GPArotation:::NormalizingWeight(Harman8, normalize = FALSE)
.check(all(W == 1), "Test 1 - FALSE returns all ones")

# --- Test 2: normalize = TRUE gives Kaiser weights ---
W.true <- GPArotation:::NormalizingWeight(Harman8, normalize = TRUE)
W.kai  <- GPArotation:::NormalizingWeight(Harman8, normalize = "Kaiser")
.check(max(abs(W.true - W.kai)) < 1e-15,
       "Test 2 - TRUE and Kaiser identical")

# --- Test 3: Kaiser rows sum to 1 after normalization ---
A_norm <- Harman8 / W.kai
.check(max(abs(rowSums(A_norm^2) - 1)) < 1e-10,
       "Test 3 - Kaiser normalized rows have unit norm")

# --- Test 4: CM weights are positive ---
W.cm <- GPArotation:::NormalizingWeight(Harman8, normalize = "CM")
.check(all(W.cm > 0), "Test 4 - CM weights all positive")

# --- Test 5: CM and Kaiser give different weights ---
.check(max(abs(W.true - W.cm)) > 1e-6,
       "Test 5 - CM and Kaiser weights differ")

# --- Test 6: custom vector ---
w <- rep(2, nrow(Harman8))
W.vec <- GPArotation:::NormalizingWeight(Harman8, normalize = w)
.check(all(W.vec == 2), "Test 6 - custom vector applied correctly")

# --- Test 7: custom function ---
W.fun <- GPArotation:::NormalizingWeight(Harman8,
           normalize = function(A) sqrt(rowSums(A^2)))
.check(max(abs(W.fun - W.kai)) < 1e-15,
       "Test 7 - custom function matches Kaiser")

# --- Test 8: wrong length vector throws error ---
tryCatch({
  GPArotation:::NormalizingWeight(Harman8, normalize = rep(1, 5))
  cat("FAILED: Test 8 - wrong length should throw error\n")
  all.ok <- FALSE
}, error = function(e) {
  # expected
})

# --- Test 9: unrecognized string throws error ---
tryCatch({
  GPArotation:::NormalizingWeight(Harman8, normalize = "unknown")
  cat("FAILED: Test 9 - unrecognized string should throw error\n")
  all.ok <- FALSE
}, error = function(e) {
  # expected
})

# --- Test 10: rotation with CM produces valid solution ---
res.cm  <- oblimin(Harman8, normalize = "CM")
res.kai <- oblimin(Harman8, normalize = TRUE)
.check(inherits(res.cm, "GPArotation"),
       "Test 10 - CM rotation returns GPArotation object")
.check(res.cm$convergence,
       "Test 11 - CM rotation converges")

cat("Normalization tests completed.\n")
if (!all.ok) stop("some normalization tests FAILED")
