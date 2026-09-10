# Independent behavior checks for the deterministic Fortran readr subset.
suppressPackageStartupMessages(library(readr))

x <- read_csv("data/sample.csv", show_col_types = FALSE)
stopifnot(nrow(x) == 3L, ncol(x) == 4L)
stopifnot(is.integer(x$id), is.logical(x$active), is.double(x$score), is.character(x$name))
stopifnot(identical(x$id, 1:3))
stopifnot(identical(x$active, c(TRUE, FALSE, TRUE)))
stopifnot(is.na(x$score[[2]]), identical(x$name[[2]], "Smith, Bob"))
stopifnot(identical(x$name[[3]], 'quote "inside"'))

stopifnot(identical(parse_logical(c("TRUE", "F", "NA")), c(TRUE, FALSE, NA)))
stopifnot(identical(parse_integer(c("-2", "+3", "NA")), c(-2L, 3L, NA_integer_)))
stopifnot(isTRUE(all.equal(parse_double(c("1,25", "-2,5"), locale = locale(decimal_mark = ",")), c(1.25, -2.5))))
stopifnot(isTRUE(all.equal(parse_number(c("$1,234.50", "EUR -20.5")), c(1234.5, -20.5))))
stopifnot(identical(guess_parser(c("1", "20")), "integer"))
stopifnot(identical(guess_parser(c("1.5", "2")), "double"))

f <- suppressWarnings(parse_factor(c("low", "high", "bad"), levels = c("low", "high"), ordered = TRUE))
stopifnot(is.ordered(f), is.na(f[[3]]), nrow(problems(f)) == 1L)

formatted <- format_csv(x)
stopifnot(grepl('2,FALSE,NA,"Smith, Bob"', formatted, fixed = TRUE))

cat("readr reference invariants passed\n")
