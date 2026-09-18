# Deterministic portmanteau test references.
options(digits = 17)
x <- c(2, 3, 5, 4, 6, 8, 7, 9, 10, 8, 7, 6, 4, 5, 3, 2, 4, 6, 5, 7)
for (method in c("Box-Pierce", "Ljung-Box")) {
    for (lag in c(1, 5, 19)) {
        for (fitdf in unique(c(0, min(2, lag - 1)))) {
            fit <- Box.test(x, lag = lag, type = method, fitdf = fitdf)
            cat(method, lag, fitdf, fit$statistic, fit$parameter, fit$p.value, "\n")
        }
    }
}
