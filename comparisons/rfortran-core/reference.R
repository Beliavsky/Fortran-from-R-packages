args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[1] else "r_results.csv"

x <- as.numeric(1:100)
points <- c(0.1, 0.5, 1, 2.5, 10, 100)
wcheck <- function(value) sum(as.numeric(value) * seq_along(value))
roll_right <- function(value, width) {
  as.numeric(stats::filter(value, rep(1 / width, width), sides = 1))
}

rows <- list()
add <- function(name, fun, reps = 20000L, atol = 1e-11, rtol = 1e-11) {
  tm <- system.time(for (i in seq_len(reps)) value <- as.numeric(fun()))[["elapsed"]]
  rows[[length(rows) + 1L]] <<- data.frame(
    case = name, value = value, seconds = tm, abs_tol = atol, rel_tol = rtol
  )
}

add("digamma_checksum", function() wcheck(digamma(points)))
add("trigamma_checksum", function() wcheck(trigamma(points)))
add("log_beta_checksum", function() wcheck(lbeta(points, rev(points))))
add("rolling_mean_valid", function() wcheck(roll_right(x, 5)[5:100]))
add("rolling_mean_right", function() sum(roll_right(x, 5)[5:100] * (5:100)))

write.csv(do.call(rbind, rows), out, row.names = FALSE, quote = FALSE)
