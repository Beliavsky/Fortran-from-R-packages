# Deterministic references for unequal and pooled variance one-way inference.
options(digits = 17)
x <- c(1, 2, 4, 10, 12, 18, 25, -4, -3, 0, 7, 9)
g <- c(rep(-7, 3), rep(0, 4), rep(1000, 5))
for (equal in c(FALSE, TRUE)) {
    fit <- oneway.test(x ~ factor(g), var.equal = equal)
    cat("three groups", equal, fit$statistic, fit$parameter, fit$p.value, "\n")
    fit <- oneway.test(x[1:7] ~ factor(g[1:7]), var.equal = equal)
    cat("two groups", equal, fit$statistic, fit$parameter, fit$p.value, "\n")
}
