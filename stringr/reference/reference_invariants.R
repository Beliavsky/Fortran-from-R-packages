# SPDX-License-Identifier: MIT

stopifnot(nchar("abc") == 3L)
stopifnot(grepl("pha", "alphabet", fixed = TRUE))
stopifnot(startsWith("alphabet", "alpha"))
stopifnot(endsWith("alphabet", "bet"))
stopifnot(lengths(regmatches("banana", gregexpr("an", "banana", fixed = TRUE))) == 2L)
stopifnot(tolower("AbC") == "abc")
stopifnot(toupper("AbC") == "ABC")
stopifnot(trimws("  abc  ") == "abc")
stopifnot(gsub("an", "X", "banana", fixed = TRUE) == "bXXa")
stopifnot(paste(c("a", "b", "c"), collapse = ",") == "a,b,c")
stopifnot(identical(strsplit("a--b--c", "--", fixed = TRUE)[[1]], c("a", "b", "c")))

if (requireNamespace("stringr", quietly = TRUE)) {
  stopifnot(stringr::str_squish("  a   b  ") == "a b")
  stopifnot(stringr::str_pad("a", 3, side = "left", pad = "0") == "00a")
  stopifnot(stringr::str_trunc("abcdef", 5) == "ab...")
  stopifnot(stringr::str_to_snake("Hello world") == "hello_world")
}

cat("Reference stringr invariants passed.\n")
