test_that("nzm docstring form overflows R integer at n=46342", {

  # This checks an overflow but gives an error in R:
  # /home/hornik/src/R/src/main/arithmetic.c:357:6: runtime error:
  # signed integer overflow: 46355 * 46342 cannot be represented in
  # type 'int', so skip on CRAN

  # BUG: The nzm formula shown in the n1qn1() docstring example uses integer
  # arithmetic: nzm=as.integer(n*(n+13L)/2L)
  # For n >= 46342, n * (n + 13L) overflows R's signed 32-bit integer
  # (INT_MAX = 2,147,483,647; 46342 * 46355 = 2,148,183,410 > INT_MAX),
  # producing NA_integer_.  The function body itself uses double arithmetic
  # (n * (n + 13) / 2 with unadorned 13 and 2), so the wrapper is safer,
  # but callers who copy the docstring pattern will hit this overflow.
  n_overflow <- 46342L
  nzm_bad <- suppressWarnings(as.integer(n_overflow * (n_overflow + 13L) / 2L))
  expect_true(is.na(nzm_bad),
    label = "nzm computed with integer arithmetic overflows to NA for n=46342")
})

test_that("C-level nd overflows int32 for n=46341: out-of-bounds zm[] access", {
  # BUG: n1qn1_all.c:144  nd = *n * (*n + 1) / 2 + 1
  # Both operands are C int (signed 32-bit).  sqrt(INT_MAX) ~ 46340.95,
  # so for n >= 46341: *n * (*n+1) overflows int32.
  # The resulting wrong nd is used on line 151 as an offset into zm[],
  # producing an out-of-bounds read/write into the workspace array.
  # Demonstrate the same overflow boundary in R integer arithmetic:
  n_val <- 46341L
  product <- suppressWarnings(n_val * (n_val + 1L))
  expect_true(is.na(product) || product < 0L,
              label = "46341 * 46342 overflows signed int32")
  skip("Confirmed segfault (exit 139): n=46341 overflows nd in n1qn1_all.c:144, then accesses zm[] out-of-bounds")
  expect_error(n1qn1(function(x) sum(x^2), function(x) 2 * x, double(n_val)),
               "n is too large")
})

test_that("as.integer(nzm) overflows to NA for n=65536: error before C allocation", {
  # BUG: R/n1qn1.R:140  nzm <- as.integer(ceiling(n * (n + 13) / 2))
  # For n >= 65536, even with double intermediate arithmetic, the result
  # exceeds .Machine$integer.max (2,147,483,647), so as.integer() returns
  # NA_integer_ with a warning.  R then tries double(NA) for zm, which
  # throws "vector size cannot be NA" before reaching C.
  # If the NA somehow reached C (e.g. via direct .Call), INTEGER(nzmSEXP)[0]
  # would yield INT_MIN = -2147483648, and new double[INT_MIN] would convert
  # to size_t ~1.8e19 -> bad_alloc or segfault on overcommit systems.
  n_val <- 65536L
  nzm_dbl <- n_val * (n_val + 13) / 2   # double arithmetic: ~2.148e9
  expect_gt(nzm_dbl, .Machine$integer.max,
    label = "nzm for n=65536 exceeds INT_MAX before as.integer()")
  nzm_int <- suppressWarnings(as.integer(ceiling(nzm_dbl)))
  expect_true(is.na(nzm_int),
    label = "as.integer(nzm) for n=65536 overflows to NA_integer_")
  # R catches the NA before any large allocation; no skip needed.
  # suppressWarnings absorbs the "NAs introduced by coercion" warning from
  # as.integer() inside n1qn1() before the error is thrown:
  expect_error(
    suppressWarnings(n1qn1(function(x) sum(x^2), function(x) 2 * x, double(n_val))),
    regexp = "vector size cannot be NA"
  )
})
