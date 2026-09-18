# Regenerate deterministic reference values for selected R stats smoothers.
options(digits = 17)

x <- c(1, 4, 2, 8, 5, 3, 9, 7, 6)
y <- c(1, 2.5, 1.2, 4.8, 3.7, 3.1, 5.9, 5.2, 4.5)
x_points <- c(1, 3, 5, 7, 9)

print(runmed(y, 3))
print(runmed(y, 5, endrule = "median"))
print(runmed(y, 5, endrule = "keep"))
print(runmed(y, 5, endrule = "constant"))
print(ksmooth(x, y, kernel = "normal", bandwidth = 2, x.points = x_points))
print(ksmooth(x, y, kernel = "box", bandwidth = 4, x.points = x_points))
print(lowess(x, y, f = 2 / 3, iter = 3, delta = 0))

x_new <- c(0.5, 2.5, 5.5, 9.5)
for (local_degree in 0:2) {
  fit <- loess(
    y ~ x,
    span = 0.75,
    degree = local_degree,
    family = "gaussian",
    control = loess.control(surface = "direct", statistics = "exact")
  )
  print(fitted(fit))
  print(predict(fit, x_new))
}

fit <- loess(
  y ~ x,
  span = 0.75,
  degree = 1,
  family = "symmetric",
  control = loess.control(surface = "direct", statistics = "exact", iterations = 4)
)
print(fitted(fit))
print(fit$robust)
print(predict(fit, x_new))

prior_weights <- c(1, 2, 0.5, 1, 3, 1.5, 0.75, 2.5, 1)
fit <- loess(
  y ~ x,
  weights = prior_weights,
  span = 0.75,
  degree = 1,
  family = "gaussian",
  control = loess.control(surface = "direct", statistics = "exact")
)
print(fitted(fit))
print(predict(fit, x_new))

fit <- smooth.spline(x, y, spar = 0.5, all.knots = TRUE)
print(c(spar = fit$spar, lambda = fit$lambda, ratio = fit$ratio, df = fit$df))
print(fit$y)
print(predict(fit, x_new)$y)
print(predict(fit, c(0.5, 1, 2.5, 5.5, 9, 9.5), deriv = 1)$y)
print(predict(fit, c(0.5, 1, 2.5, 5.5, 9, 9.5), deriv = 2)$y)

fit <- smooth.spline(x, y, w = prior_weights, spar = 0.5, all.knots = TRUE)
print(c(spar = fit$spar, lambda = fit$lambda, ratio = fit$ratio, df = fit$df))
print(fit$y)

fit <- smooth.spline(x, y, df = 4, all.knots = TRUE)
print(c(spar = fit$spar, lambda = fit$lambda, ratio = fit$ratio, df = fit$df))
print(fit$y)
print(predict(fit, x_new)$y)

fit <- smooth.spline(x, y, all.knots = TRUE)
print(c(spar = fit$spar, lambda = fit$lambda, ratio = fit$ratio, df = fit$df,
        gcv = fit$cv.crit))
print(fit$y)
print(predict(fit, x_new)$y)

fit <- smooth.spline(x, y, cv = TRUE, all.knots = TRUE)
print(c(spar = fit$spar, lambda = fit$lambda, ratio = fit$ratio, df = fit$df,
        cv = fit$cv.crit))
print(fit$y)
print(predict(fit, x_new)$y)

fit <- smooth.spline(x, y, w = prior_weights, spar = 0.5, cv = TRUE,
                     all.knots = TRUE)
print(c(spar = fit$spar, lambda = fit$lambda, ratio = fit$ratio, df = fit$df,
        cv = fit$cv.crit))
print(fit$y)
print(predict(fit, x_new)$y)

duplicate_x <- c(1, 1, 2, 3, 3, 3, 4, 5, 6, 7)
duplicate_y <- c(1, 1.4, 1.8, 3, 3.4, 2.8, 4.2, 4.8, 5.5, 7)
duplicate_weights <- c(1, 2, 1, 0.5, 2, 1.5, 1, 1, 2, 1)
fit <- smooth.spline(
  duplicate_x,
  duplicate_y,
  w = duplicate_weights,
  spar = 0.5,
  all.knots = TRUE
)
print(c(spar = fit$spar, lambda = fit$lambda, ratio = fit$ratio, df = fit$df))
print(fit$x)
print(fit$yin)
print(fit$w)
print(fit$y)
print(predict(fit, c(0.5, 1, 2.5, 7.5))$y)

fit <- suppressWarnings(smooth.spline(
  duplicate_x,
  duplicate_y,
  w = duplicate_weights,
  spar = 0.5,
  cv = TRUE,
  all.knots = TRUE
))
print(c(spar = fit$spar, lambda = fit$lambda, ratio = fit$ratio, df = fit$df,
        cv = fit$cv.crit))
print(fit$y)

near_duplicate_x <- c(1, 1.0000002, 2, 3, 4, 5)
near_duplicate_y <- c(1, 1.4, 2, 3.2, 3.8, 5.1)
fit <- smooth.spline(
  near_duplicate_x,
  near_duplicate_y,
  spar = 0.5,
  all.knots = TRUE
)
print(c(tol = fit$tol, spar = fit$spar, lambda = fit$lambda,
        ratio = fit$ratio, df = fit$df))
print(fit$x)
print(fit$yin)
print(fit$w)
print(fit$y)
print(predict(fit, c(0.5, 1, 2.5, 5.5))$y)

large_x <- 1:64
large_y <- sin(large_x / 7) + 0.015 * large_x + 0.2 * cos(large_x / 3)
fit <- smooth.spline(large_x, large_y, spar = 0.5)
print(c(spar = fit$spar, lambda = fit$lambda, ratio = fit$ratio, df = fit$df,
        knot_vector_length = length(fit$fit$knot)))
print(fit$y[c(1, 8, 16, 32, 48, 64)])
print(predict(fit, c(0.5, 10.5, 33.3, 64.5))$y)

fit <- smooth.spline(large_x, large_y, spar = 0.5, nknots = 20)
print(c(spar = fit$spar, lambda = fit$lambda, ratio = fit$ratio, df = fit$df,
        knot_vector_length = length(fit$fit$knot)))

fit <- smooth.spline(large_x, large_y, spar = 0.5, all.knots = TRUE)
print(c(spar = fit$spar, lambda = fit$lambda, ratio = fit$ratio, df = fit$df,
        knot_vector_length = length(fit$fit$knot)))
