# Copyright (c) Rob Carnell 2026

fixsmall <- function(x) {
  z <- ifelse(x < .Machine$double.eps^(2/3), "OK", "NZ")
  z <- ifelse(is.na(z), "NA", z)
  z
}

dslnex <- function(x) {
  y <- numeric(2)
  y[1] <- x[1]^2 + x[2]^2 - 2
  y[2] <- exp(x[1] - 1) + x[2]^3 - 2
  y
}

# --- Tests -------------------------------------------------------------------

test_that("testnslv works and normalizes Fnorm for a well-behaved start", {

  xstart <- c(0.5, 0.5)

  z <- testnslv(xstart, dslnex)

  # Ensure structure is present
  expect_true(is.list(z))
  expect_true(is.data.frame(z$out))
  expect_true("Fnorm" %in% colnames(z$out))

  # Normalize Fnorm and ensure only OK/NZ/NA appear
  zfn <- z$out[, "Fnorm"]
  z$out[, "Fnorm"] <- fixsmall(zfn)

  expect_true(all(z$out[, "Fnorm"] %in% c("OK", "NZ", "NA")))
  expect_true(any(z$out[, "Fnorm"] == "OK"))

  # All methods should have converged (termcd == 1)
  expect_true(all(z$out[, "termcd"] == 1))
})

test_that("testnslv handles a problematic start and records the Newton/none error", {

  xstart <- c(2.0, 0.5)

  expect_output(z <- testnslv(xstart, dslnex),
                regexp = "Error")

  # Ensure structure is present
  expect_true(is.data.frame(z$out))
  expect_true("Fnorm" %in% colnames(z$out))

  # Normalize Fnorm
  zfn <- z$out[, "Fnorm"]
  z$out[, "Fnorm"] <- fixsmall(zfn)

  # Check that the Newton/none row contains NA and ERROR as in the example
  idx <- which(z$out[, "Method"] == "Newton" & z$out[, "Global"] == "none")
  expect_length(idx, 1)

  expect_true(is.na(z$out[idx, "termcd"]))
  expect_identical(z$out[idx, "Message"], "ERROR")
  expect_identical(z$out[idx, "Fnorm"], "NA")

  # Other methods should still have valid termcd values
  ok_rows <- setdiff(seq_len(nrow(z$out)), idx)
  expect_true(all(z$out[ok_rows, "termcd"] %in% c(1, 4)))
})

