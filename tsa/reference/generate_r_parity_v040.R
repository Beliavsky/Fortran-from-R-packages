# TSA-fortran v0.4.0 R parity fixture generator
# SPDX-License-Identifier: GPL-2.0-or-later
#
# Run with an R installation that has TSA available:
#   Rscript reference/generate_r_parity_v040.R
#
# The deterministic datasets mirror the Fortran v0.4.0 xreg/Hessian tests.

options(digits = 17)

n <- 80L
i <- seq_len(n)
t <- i / n
x1 <- t + 0.15 * sin(0.37 * i)
x2 <- 1.0002 * x1 + 0.012 * cos(0.29 * i)
xreg <- cbind(x1 = x1, x2 = x2)
y <- 2.25 * x1 - 1.70 * x2 + 0.04 * sin(0.83 * i)

cat("stats::arima all-free xreg\n")
fit <- stats::arima(y, order = c(0L, 0L, 0L), xreg = xreg,
                    include.mean = FALSE, method = "ML")
print(fit$coef)
print(fit$sigma2)
print(fit$var.coef)
print(fit$loglik)

cat("\nstats::arima fixed second xreg coefficient\n")
fit.fixed <- stats::arima(y, order = c(0L, 0L, 0L), xreg = xreg,
                          include.mean = FALSE, method = "ML",
                          fixed = c(NA_real_, -1.70))
print(fit.fixed$coef)
print(fit.fixed$sigma2)
print(fit.fixed$var.coef)

cat("\nstats::arima missing y/xreg\n")
y.miss <- y
xreg.miss <- xreg
y.miss[13L] <- NA_real_
xreg.miss[27L, 1L] <- NA_real_
xreg.miss[51L, 2L] <- NA_real_
fit.miss <- stats::arima(y.miss, order = c(0L, 0L, 0L), xreg = xreg.miss,
                         include.mean = FALSE, method = "ML")
print(fit.miss$coef)
print(fit.miss$sigma2)
print(fit.miss$var.coef)
print(fit.miss$loglik)

cat("\nTSA::arimax ordinary xreg delegation\n")
if (requireNamespace("TSA", quietly = TRUE)) {
    fit.tsa <- TSA::arimax(y, order = c(0L, 0L, 0L), xreg = xreg,
                           include.mean = FALSE, method = "ML")
    print(fit.tsa$coef)
    print(fit.tsa$sigma2)
    print(fit.tsa$var.coef)
    print(fit.tsa$loglik)
} else {
    cat("TSA not installed; skipped.\n")
}
