library(vctrs)

stopifnot(typeof(vec_ptype_common(TRUE, 1L, 2.0)) == "double")
stopifnot(identical(vec_recycle(7L, 0L), integer()))

x <- c(1L, NA_integer_)
stopifnot(identical(vec_equal(x, x), c(TRUE, NA)))
stopifnot(identical(vec_equal(x, x, na_equal = TRUE), c(TRUE, TRUE)))

locations <- vec_match(c(3L, 1L, 9L), c(1L, 2L, 3L))
stopifnot(identical(locations, c(3L, 1L, NA_integer_)))

cat("Reference vctrs invariants passed.\n")
