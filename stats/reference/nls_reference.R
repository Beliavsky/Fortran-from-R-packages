# Regenerate deterministic reference values for nonlinear least squares.
options(digits = 17)

x <- 0:8
y <- c(2.1, 3, 4.7, 6.8, 10.2, 14.9, 22.3, 32.8, 49.1)
fit <- nls(
  y ~ a * exp(b * x),
  start = list(a = 2, b = 0.35),
  control = nls.control(tol = 1e-8, maxiter = 200)
)
print(coef(fit))
print(c(rss = deviance(fit), sigma = summary(fit)$sigma,
        iterations = fit$convInfo$finIter, tolerance = fit$convInfo$finTol))
print(fitted(fit))
print(residuals(fit))
print(vcov(fit))
print(confint(fit, level = 0.95))

x <- 0:10
y <- 1.5 + 4 * exp(-0.35 * x) +
  c(0.02, -0.03, 0.01, 0.04, -0.02, 0.01, -0.01, 0.02, -0.02, 0.01, -0.01)
fit <- nls(
  y ~ cbind(1, exp(-k * x)),
  start = list(k = 0.2),
  algorithm = "plinear"
)
print(coef(fit))
print(c(rss = deviance(fit), sigma = summary(fit)$sigma,
        iterations = fit$convInfo$finIter, tolerance = fit$convInfo$finTol))
print(fitted(fit))
print(residuals(fit))
print(vcov(fit))

x <- 1:10
y <- c(1.7, 2.8, 3.6, 4.2, 4.8, 5.1, 5.5, 5.7, 6, 6.2)
fit <- nls(
  y ~ vm * x / (k + x),
  start = list(vm = 8, k = 3),
  control = nls.control(tol = 1e-8, maxiter = 200)
)
print(coef(fit))
print(c(rss = deviance(fit), sigma = summary(fit)$sigma,
        iterations = fit$convInfo$finIter, tolerance = fit$convInfo$finTol))
print(fitted(fit))
print(residuals(fit))
print(vcov(fit))

bounded_fit <- nls(
  y ~ vm * x / (k + x),
  start = list(vm = 8, k = 3),
  algorithm = "port",
  lower = c(0, 0),
  upper = c(8.5, Inf),
  control = nls.control(tol = 1e-10, maxiter = 200)
)
print(coef(bounded_fit))
print(c(rss = deviance(bounded_fit), sigma = summary(bounded_fit)$sigma,
        iterations = bounded_fit$convInfo$finIter))
print(fitted(bounded_fit))
print(vcov(bounded_fit))

# Pure curve-value fixtures for the translated self-start model helpers.
print(c(
  SSasymp = SSasymp(2, 10, 2, log(0.4)),
  SSgompertz = SSgompertz(2, 10, 2, 0.7),
  SSlogis = SSlogis(2, 10, 3, 0.8),
  SSmicmen = SSmicmen(2, 10, 3),
  SSweibull = SSweibull(2, 10, 8, log(0.4), 1.5)
))
