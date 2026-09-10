# Independent behavior checks for the deterministic Fortran forcats subset.
suppressPackageStartupMessages(library(forcats))

f <- factor(c("b", "b", "a", "c", "c", "c"), levels = c("a", "b", "c", "d"))

stopifnot(identical(levels(fct_inorder(f)), c("b", "a", "c", "d")))
stopifnot(identical(levels(fct_infreq(f)), c("c", "b", "a", "d")))
stopifnot(identical(levels(fct_rev(f)), c("d", "c", "b", "a")))
stopifnot(identical(levels(fct_shift(f, -1)), c("d", "a", "b", "c")))
stopifnot(identical(levels(fct_relevel(f, "c", "a", after = 1)), c("b", "c", "a", "d")))
stopifnot(identical(levels(fct_expand(f, "e", "b", "f", after = 1)), c("a", "e", "f", "b", "c", "d")))
stopifnot(identical(levels(fct_drop(f)), c("a", "b", "c")))

recode <- fct_recode(f, alpha = "a", alpha = "c")
stopifnot(identical(levels(recode), c("alpha", "b", "d")))
other <- fct_other(f, keep = c("b", "c"))
stopifnot(identical(levels(other), c("b", "c", "Other")))
stopifnot(identical(as.vector(fct_match(f, c("a", "c"))), c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE)))

g <- factor(c("a", NA, "b"), levels = c("a", "b"))
explicit <- fct_na_value_to_level(g, "Missing")
stopifnot(identical(levels(explicit), c("a", "b", "Missing")))
stopifnot(!anyNA(explicit))
implicit <- fct_na_level_to_value(explicit, "Missing")
stopifnot(identical(levels(implicit), c("a", "b")), is.na(implicit[[2]]))

combined <- fct_c(factor(c("a", "b")), factor(c("b", "c")))
stopifnot(identical(levels(combined), c("a", "b", "c")))
stopifnot(identical(as.integer(combined), c(1L, 2L, 2L, 3L)))

fruit <- factor(c("apple", "kiwi", "apple", "apple"))
colour <- factor(c("green", "green", "red", "green"))
crossed <- fct_cross(fruit, colour)
stopifnot(identical(levels(crossed), c("apple:green", "apple:red", "kiwi:green")))
stopifnot(identical(as.integer(crossed), c(1L, 3L, 2L, 1L)))

h <- factor(rep(letters[1:5], c(4, 2, 1, 1, 1)), levels = c(letters[1:5], "z"))
stopifnot(identical(levels(fct_lump_min(h, 2)), c("a", "b", "Other")))
stopifnot(identical(levels(fct_lump_n(h, 2)), c("a", "b", "Other")))
stopifnot(identical(levels(fct_lump_prop(h, 0.15)), c("a", "b", "Other")))

x <- c(4, 6, 8, 10, 1, 3, 20, 12, 16)
stopifnot(identical(levels(fct_reorder(h, x, mean)), c("b", "a", "d", "e", "c", "z")))

cat("forcats reference invariants passed\n")
