# Regenerate deterministic fixtures in test/test_spline_interpolation.f90.

x <- c(0, 1, 2, 4, 5)
y <- c(0, 1, 0, 2, 1)
xout <- c(-1, 0.5, 1.5, 3, 5, 6)

for (method in c("fmm", "natural")) {
  fit <- splinefun(x, y, method = method)
  cat(toupper(method), "VALUE", format(fit(xout), digits = 17), "\n")
  cat(toupper(method), "D1", format(fit(xout, deriv = 1), digits = 17), "\n")
  cat(toupper(method), "D2", format(fit(xout, deriv = 2), digits = 17), "\n")
}

ties <- spline(c(2, 0, 1, 1, 3), c(4, 0, 1, 3, 9), xout = c(0.5, 1, 1.5),
               method = "natural", ties = mean)
cat("TIED_NATURAL", format(ties$y, digits = 17), "\n")
