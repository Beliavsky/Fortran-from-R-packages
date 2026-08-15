# Copyright (c) 2013, Avraham Adler All rights reserved
# SPDX-License-Identifier: BSD-2-Clause

# For CRAN
oldThreads <- getDelapThreads()
setDelapThreads(2L)

tol <- sqrt(.Machine$double.eps)
VAL <- data.frame(read.csv(file = file.path(".", "RawTest.csv"), header = TRUE))
nanWarn <- "NaNs produced"

# Singleton function accuracy
expect_equal(pdelap(0:36, 2, 1, 5), VAL$PDELAP_2, tolerance = tol)

# Singleton log.p
expect_equal(pdelap(0:36, 4, 5, 1, log.p = TRUE), log(pdelap(0:36, 4, 5, 1)),
             tolerance = tol)

# Singleton lower.tail
expect_equal(pdelap(0:100, 8, 10, 6, lower.tail = FALSE),
             1 - pdelap(0:100, 8, 10, 6), tolerance = tol)
expect_equal(pdelap(6, 2.9647, 0.005 / 2.9647, 0.0057, lower.tail = FALSE), 0,
             tolerance = tol)

# Singleton lower.tail & log.p
expect_equal(pdelap(0:100, 8, 10, 6, lower.tail = FALSE, log.p = TRUE),
             log(1 - pdelap(0:100, 8, 10, 6)), tolerance = tol)

# Singleton NaN
expect_warning(pdelap(1, 0, 1, 2), nanWarn)
expect_warning(pdelap(1, -85, 1, 2), nanWarn)
expect_warning(pdelap(0, 1, 0, 2), nanWarn)
expect_warning(pdelap(0, 1, -1e4, 2), nanWarn)
expect_warning(pdelap(0, 1, 4, 0), nanWarn)
expect_warning(pdelap(0, 1, 4, -1e-3), nanWarn)
expect_warning(ddelap(NaN, 1, 4, 6), nanWarn)
expect_warning(ddelap(NA, 1, 4, 12), nanWarn)
tst <- suppressWarnings(pdelap(c(NA, 4, NaN), 0.5, 4, 0.2))
expect_equal(tst, c(NaN, pdelap(4, 0.5, 4, 0.2), NaN), tolerance = tol)
expect_identical(suppressWarnings(pdelap(c(NA, NaN), 0.5, 4, 0.2)), rep(NaN, 2))

# Vector function accuracy
expect_equal(pdelap(0:36, c(1, 2, 3), c(4, 1, 2), c(2, 5, 7)),
             VAL$PDELAP_Triple, tolerance = tol)

# Vector log.p
expect_equal(pdelap(0:36, c(1, 2, 3), c(4, 1, 2), c(2, 5, 7), log.p = TRUE),
             log(VAL$PDELAP_Triple), tolerance = tol)

# Vector lower.tail
expect_equal(pdelap(0:100, c(4, 9, 2), c(6, 12, 8), c(7, 14, 9),
                    lower.tail = FALSE),
             1 - pdelap(0:100, c(4, 9, 2), c(6, 12, 8), c(7, 14, 9)),
             tolerance = tol)

# Vector lower.tail & log.p
expect_equal(pdelap(0:100, c(4, 9, 2), c(6, 12, 8), c(7, 14, 9),
                    lower.tail = FALSE, log.p = TRUE),
             log(1 - pdelap(0:100, c(4, 9, 2), c(6, 12, 8), c(7, 14, 9))),
             tolerance = tol)

# Vector NaN
expect_warning(pdelap(1:10, 0, 1, 2), nanWarn)
expect_warning(pdelap(1:3, c(0, 1, 2), c(1, 0, 2), c(1, 2, 0)), nanWarn)
expect_warning(pdelap(1:3, c(-5e5, 1, 2), c(1, -2, 2), c(1, 2, -8e-4)), nanWarn)
expect_warning(pdelap(c(NA, 2), c(2, 1, 2), c(1, 3, 2), c(1, 2, 4)))
expect_warning(pdelap(c(3, NaN), c(2, 1, 2), c(1, 3, 2), c(1, 2, 4)))
tst <- suppressWarnings(pdelap(c(NA, 0, NaN), c(1, 2, 3), c(4, 1, 2),
                               c(2, 5, 7)))
expect_equal(tst, c(NaN, pdelap(0, 2, 1, 5), NaN), tolerance = tol)
tst <- suppressWarnings(pdelap(c(0, 0, 0), c(NA, 2, 3), c(4, 1, 2),
                               c(2, 5, NaN)))
expect_equal(tst, c(NaN, pdelap(0, 2, 1, 5), NaN), tolerance = tol)

# Negative values due to floating point issues are 0
expect_equal(pdelap(500, 13.08251, 0.02414521, 0.04421658, FALSE, FALSE), 0,
             tolerance = tol)

# Non-double parameters converted
expect_equal(pdelap(2L, 1L, 2L, 3L), pdelap(2L, 1, 2, 3), tolerance = tol)

# Floating point issues do not lead to CDF > 1
# print(pdelap(1000, 8, 15, 100), digits = 17) used to be 1.0000000000001035
expect_true(pdelap(1000, 8, 15, 100) <= 1)

# Infinite values
expect_identical(pdelap(Inf, 1L, 2L, 3L), 1)
expect_identical(pdelap(c(Inf, Inf), c(1L, 2L), 2L, 3L), c(1, 1))
expect_warning(pdelap(-Inf, 1L, 2L, 3L), nanWarn)

# Restore original thread count
setDelapThreads(oldThreads)
