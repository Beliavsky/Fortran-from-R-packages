suppressPackageStartupMessages(library(tidyr))

wide <- data.frame(id = c(1L, 2L), x = c(10, 20), y = c(11, 21))
long <- pivot_longer(wide, c(x, y))
stopifnot(identical(long$id, c(1L, 1L, 2L, 2L)))
stopifnot(identical(long$name, c("x", "y", "x", "y")))
stopifnot(identical(long$value, c(10, 11, 20, 21)))

slow <- pivot_longer(wide, c(x, y), cols_vary = "slowest")
stopifnot(identical(slow$id, c(1L, 2L, 1L, 2L)))
stopifnot(identical(slow$name, c("x", "x", "y", "y")))
stopifnot(identical(slow$value, c(10, 20, 11, 21)))
round_trip <- pivot_wider(long, names_from = name, values_from = value)
stopifnot(identical(round_trip$id, wide$id))
stopifnot(identical(round_trip$x, wide$x))
stopifnot(identical(round_trip$y, wide$y))

missing_data <- data.frame(
  a = c(1L, NA_integer_, 3L),
  b = c("x", "y", NA_character_),
  w = c(1L, 0L, 2L)
)
stopifnot(nrow(drop_na(missing_data)) == 1L)
stopifnot(identical(drop_na(missing_data, a)$a, c(1L, 3L)))
stopifnot(identical(replace_na(missing_data, list(a = 9L))$a, c(1L, 9L, 3L)))
filled <- fill(missing_data, a, b, .direction = "downup")
stopifnot(identical(filled$a, c(1L, 1L, 3L)))
stopifnot(identical(filled$b, c("x", "y", "y")))
stopifnot(identical(uncount(missing_data, w)$a, c(1L, 3L, 3L)))

parts <- separate_wider_delim(data.frame(id = 1:2, code = c("a-b", "c-d")), code, "-",
  names = c("left", "right"))
stopifnot(identical(parts$left, c("a", "c")))
stopifnot(identical(parts$right, c("b", "d")))
joined <- unite(parts, code, left, right, sep = "-")
stopifnot(identical(joined$code, c("a-b", "c-d")))

grid <- expand_grid(x = c(2L, 1L), y = c("b", "a"))
stopifnot(identical(grid$x, c(2L, 2L, 1L, 1L)))
stopifnot(identical(grid$y, c("b", "a", "b", "a")))
sorted_grid <- crossing(x = c(2L, 1L, 2L), y = c("b", "a"))
stopifnot(identical(sorted_grid$x, c(1L, 1L, 2L, 2L)))
stopifnot(identical(sorted_grid$y, c("a", "b", "a", "b")))

cat("Reference tidyr invariants passed.\n")
