# Regenerate deterministic fixtures in test/test_interpolation.f90.

x <- c(4, 1, 2, 2, 6)
y <- c(8, 1, 4, 6, 12)
xout <- c(0, 1, 1.5, 2, 3, 6, 7)

linear <- approx(x, y, xout = xout, ties = mean, rule = 2)
cat("APPROX_LINEAR", format(linear$y, digits = 17), "\n")
constant <- approx(x, y, xout = xout, method = "constant", f = 0, ties = mean, rule = 2)
cat("APPROX_CONSTANT", format(constant$y, digits = 17), "\n")
mixed <- approx(x, y, xout = xout, method = "constant", f = 0.25, ties = mean, rule = 2)
cat("APPROX_MIXED", format(mixed$y, digits = 17), "\n")

iso <- isoreg(1:8, c(3, 2, 4, 5, 1, 6, 8, 7))
cat("ISOREG_FITTED", format(iso$yf, digits = 17), "\n")
cat("ISOREG_KNOTS", iso$iKnots, "\n")
