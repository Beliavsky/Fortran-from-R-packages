# Deterministic R-side parity fixture for TSA-fortran v0.5.0.
# Uses base stats only; no plotting.

n <- 500L
u <- ((seq_len(n) * 0.6180339887498949) %% 1) - 0.5
x <- numeric(n)
x[1:2] <- c(0.1, -0.1)
for (i in 3:n) x[i] <- 0.65*x[i-1] - 0.20*x[i-2] + 0.35*u[i]

cat("# fixed-order AR fits\n")
for (method in c("yule-walker", "burg", "ols", "mle")) {
  fit <- ar(x, aic=FALSE, order.max=2, method=method)
  sp <- spec.ar(x, n.freq=33, order=2, method=method, plot=FALSE)
  cat(method, "order", fit$order, "var.pred", format(fit$var.pred, digits=17), "\n")
  cat("coef", paste(format(fit$ar, digits=17), collapse=" "), "\n")
  cat("spec", paste(format(sp$spec, digits=17), collapse=" "), "\n")
}

fit_b2 <- ar.burg(x, aic=FALSE, order.max=2, var.method=2)
sp_b2 <- spec.ar(x, n.freq=33, order=2, method="burg", var.method=2, plot=FALSE)
cat("burg2 var.pred", format(fit_b2$var.pred, digits=17), "\n")
cat("burg2 coef", paste(format(fit_b2$ar, digits=17), collapse=" "), "\n")
cat("burg2 spec", paste(format(sp_b2$spec, digits=17), collapse=" "), "\n")

cat("# AIC-selected orders\n")
for (method in c("yule-walker", "burg", "ols", "mle")) {
  fit <- ar(x, aic=TRUE, order.max=if (method == "mle") 4 else 6, method=method)
  cat(method, fit$order, "\n")
}

cat("# compact kernel / periodogram\n")
k <- kernel("modified.daniell", m=2)
pg <- spec.pgram(x, kernel=k, taper=0.1, detrend=TRUE, plot=FALSE)
cat("kernel.coef", paste(format(k$coef, digits=17), collapse=" "), "\n")
cat("pgram", paste(format(pg$spec, digits=17), collapse=" "), "\n")
cat("df", format(pg$df, digits=17), "bandwidth", format(pg$bandwidth, digits=17), "\n")

xm <- cbind(x, 0.8*x + 0.2*c(x[-1], x[1]))
pm <- spec.pgram(xm, spans=5, taper=c(0,0.2), detrend=TRUE, plot=FALSE)
cat("matrix.spec1", paste(format(pm$spec[,1], digits=17), collapse=" "), "\n")
cat("matrix.spec2", paste(format(pm$spec[,2], digits=17), collapse=" "), "\n")
cat("matrix.coh", paste(format(pm$coh[,1], digits=17), collapse=" "), "\n")
cat("matrix.phase", paste(format(pm$phase[,1], digits=17), collapse=" "), "\n")
