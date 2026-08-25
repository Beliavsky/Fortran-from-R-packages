# Tests to ensure single-factor models produce appropriate errors.

require("GPArotation")

all.ok <- TRUE

expected_msg <- "Rotation does not make sense for single-factor models."

# 1-factor model as a vector
xv <- runif(5)

# --- Test 1: wrapper function rejects a vector ---
y <- try(GPArotation::quartimin(xv), silent = TRUE)
if (!inherits(y, "try-error")) {
  cat("Error messages: test 1 failed\n")
  all.ok <- FALSE
}

# --- Test 2: GPForth rejects a vector ---
y <- try(GPForth(xv, method = "quartimax"), silent = TRUE)
if (!inherits(y, "try-error")) {
  cat("Error messages: test 2 failed\n")
  all.ok <- FALSE
}

# --- Test 3: GPFoblq rejects a vector ---
y <- try(GPFoblq(xv, method = "quartimin"), silent = TRUE)
if (!inherits(y, "try-error")) {
  cat("Error messages: test 3 failed\n")
  all.ok <- FALSE
}

# 1-factor model as a single-column matrix
xw <- matrix(xv)

# --- Test 4: GPForth rejects a single-column matrix with correct message ---
y <- try(GPForth(xw, method = "quartimax"), silent = TRUE)
if (!inherits(y, "try-error") ||
    !grepl(expected_msg, attr(y, "condition")$message)) {
  cat("Error messages: test 4 failed\n")
  all.ok <- FALSE
}

# --- Test 5: GPFoblq rejects a single-column matrix with correct message ---
y <- try(GPFoblq(xw, method = "quartimin"), silent = TRUE)
if (!inherits(y, "try-error") ||
    !grepl(expected_msg, attr(y, "condition")$message)) {
  cat("Error messages: test 5 failed\n")
  all.ok <- FALSE
}

# --- Test 6: wrapper rejects a single-column matrix with correct message ---
y <- try(GPArotation::quartimin(xw), silent = TRUE)
if (!inherits(y, "try-error") ||
    !grepl(expected_msg, attr(y, "condition")$message)) {
  cat("Error messages: test 6 failed\n")
  all.ok <- FALSE
}

cat("tests completed.\n")

if (!all.ok) stop("some tests FAILED")