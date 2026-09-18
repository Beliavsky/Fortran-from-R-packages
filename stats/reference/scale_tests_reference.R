# Deterministic references for the median-centered Fligner-Killeen test.
options(digits = 17)
x <- c(4, 5, 6, 5, 7, 8, 3, 4, 5, 6, 6, 7, 8, 9, 11)
g <- rep(c(-7, 0, 1000), each = 5)
emit <- function(label, x, g) {
    fit <- fligner.test(x, g)
    cat(label, fit$statistic, fit$parameter, fit$p.value, "\n")
}
emit("ties", x, g)
emit("unequal", c(1, 2, 4, 10, 12, 18, 25, -4, -3, 0, 7, 9),
     c(rep(1, 3), rep(2, 4), rep(3, 5)))
emit("singleton", c(1, 3, 6, 9, 11), c(1, 2, 2, 2, 2))
emit("zero_spread", c(1, 1, 1, 4, 4, 4), c(1, 1, 1, 2, 2, 2))
