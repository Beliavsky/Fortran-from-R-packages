maInfty <- function(ar, ma, m = 1000) {
  p <- length(ar)
  q <- length(ma)
  if (m - q > 0) {
    times <- m - q
  } else {
    times <- 0
  }
  ma.s <- c(ma, rep(0, times = times))
  c.out <- c(rep(0, times = p - 1), 1, rep(NA, times = m))
  lc <- length(c.out)

  if (m >= 1) {
    for (i in (p + 1):(m + p)) {
      c.out[i] <- ar %*% c.out[(i - 1):(i - p)] + ma.s[i - p]
    }
  }

  c.out[p:lc]
}
