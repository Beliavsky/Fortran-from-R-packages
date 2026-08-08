next_subset <- function(a, n, k) {
    if (a[1L] == n - k +1)
        return(NULL)
    k1 <- 1:k
    h <- min(which(a[k + 1 - k1] != n + 1 - k1))
    m1 <- a[k + 1 - h]
    h1 <- 1:h
    a[k + h1 - h] <- m1 + h1
    a
}
