suppressPackageStartupMessages(library(dplyr))

data <- tibble(
  x = c(1, 4, 6, NA_real_),
  y = c(5L, 12L, 8L, 1L),
  sector = c("energy", "tech", "energy", "energy")
)
filtered <- filter(data, x > 3 & y <= 10)
stopifnot(identical(filtered$y, 8L))
mutated <- mutate(data, z = x + y * 2)
stopifnot(abs(mutated$z[3] - 22) < 1e-12, is.na(mutated$z[4]))
stopifnot(identical(arrange(data, desc(y))$y, c(12L, 8L, 5L, 1L)))
stopifnot(nrow(distinct(data, sector, .keep_all = TRUE)) == 2L)

x <- tibble(id = c(1L, 2L, 2L), label = c("a", "b", "c"))
y <- tibble(id = c(2L, 2L, 3L), score = c(10, 20, 30))
stopifnot(nrow(inner_join(x, y, by = "id", relationship = "many-to-many")) == 4L)
stopifnot(nrow(left_join(x, y, by = "id", relationship = "many-to-many")) == 5L)
stopifnot(nrow(right_join(x, y, by = "id", relationship = "many-to-many")) == 5L)
stopifnot(nrow(full_join(x, y, by = "id", relationship = "many-to-many")) == 6L)
stopifnot(nrow(semi_join(x, y, by = "id")) == 2L)
stopifnot(nrow(anti_join(x, y, by = "id")) == 1L)

grouped <- tibble(g = c("a", "a", "b"), x = c(1L, 3L, 4L)) |>
  group_by(g) |>
  summarise(n = n(), total = sum(x), average = mean(x), .groups = "drop")
stopifnot(identical(grouped$n, c(2L, 1L)))
stopifnot(identical(grouped$total, c(4L, 4L)))
stopifnot(identical(grouped$average, c(2, 4)))

stopifnot(identical(lag(c(1L, 3L, 4L)), c(NA_integer_, 1L, 3L)))
stopifnot(identical(lead(c(1L, 3L, 4L)), c(3L, 4L, NA_integer_)))
stopifnot(identical(coalesce(c(1L, NA_integer_, 3L), 9L), c(1L, 9L, 3L)))
stopifnot(identical(consecutive_id(c(1L, 1L, 2L, 2L), c("a", "b", "b", "b")),
  c(1L, 2L, 3L, 3L)))

cat("Reference dplyr invariants passed.\n")
