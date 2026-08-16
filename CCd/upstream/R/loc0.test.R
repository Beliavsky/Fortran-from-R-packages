loc0.test <- function(y, tol = 1e-7) {
  lik0 <- CCd::cc.mle0(y, tol = tol)$loglik
  lik1 <- CCd::cc.mle(y)$loglik
  stat <- 2 * lik1 - 2 * lik0
  pvalue <- pchisq(stat, 1, lower.tail = FALSE)
  res <- c(stat, pvalue)
  names(res) <- c("test statistic", "p-value")
  res
}
