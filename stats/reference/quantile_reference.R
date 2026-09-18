# Deterministic numerical references for all nine unweighted estimators.
options(digits = 17)
x <- c(8, 1, 2, 2, 16, 4, 32, 7)
p <- c(0, 0.0625, 0.125, 0.25, 0.3, 0.5, 0.75, 0.9375, 1)
for (kind in 1:9) {
    cat('type', kind, '\n')
    cat(quantile(x, p, type = kind, names = FALSE), '\n')
    cat('IQR', IQR(x, type = kind), '\n')
}
cat('type 3 nearest-even:', quantile(1:4, c(0.375, 0.625, 0.875),
                                    type = 3, names = FALSE), '\n')
