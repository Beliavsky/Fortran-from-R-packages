plot.loadings <- function(x, n = 5, 
      k = ncol(x), mfrow = c(k, 1), ...) {
    if (k > 1) {
        opar <- graphics::par(mfrow = mfrow)
        on.exit(graphics::par(opar))
    }
    f <- function(i) 
      graphics::barplot(sort(x[, i],
            decreasing = TRUE)[1:n],
            ylim = 0:1, main = colnames(x)[i],
            ...)
    r <- lapply(1:k, f)
}

