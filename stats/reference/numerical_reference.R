# Deterministic R 4.6.0 fixtures for scalar numerical methods.
root <- uniroot(function(x) cos(x) - x, c(0, 1), tol = 1e-12)
extended <- uniroot(function(x) x - 10, c(0, 1), extendInt = "yes", tol = 1e-12)
minimum <- optimize(function(x) (x - 2.25)^2 + 1.5, c(-4, 6), tol = 1e-10)
maximum <- optimize(sin, c(0, 4), maximum = TRUE, tol = 1e-10)
gaussian <- integrate(function(x) exp(-x*x), -2, 2,
                      abs.tol = 1e-11, rel.tol = 1e-11)
sine <- integrate(sin, 0, pi, abs.tol = 1e-12, rel.tol = 1e-12)
right_tail <- integrate(function(x) exp(-x), 0, Inf,
                        abs.tol = 1e-10, rel.tol = 1e-10)
left_tail <- integrate(function(x) exp(x), -Inf, 0,
                       abs.tol = 1e-10, rel.tol = 1e-10)
whole_gaussian <- integrate(function(x) exp(-x*x), -Inf, Inf,
                            abs.tol = 1e-10, rel.tol = 1e-10)

cat(sprintf("UNIROOT %.17g %.17g\n", root$root, root$f.root))
cat(sprintf("UNIROOT_EXTENDED %.17g %.17g\n", extended$root, extended$f.root))
cat(sprintf("OPTIMIZE_MIN %.17g %.17g\n", minimum$minimum, minimum$objective))
cat(sprintf("OPTIMIZE_MAX %.17g %.17g\n", maximum$maximum, maximum$objective))
cat(sprintf("INTEGRATE_GAUSSIAN %.17g %.17g\n", gaussian$value, gaussian$abs.error))
cat(sprintf("INTEGRATE_SINE %.17g %.17g\n", sine$value, sine$abs.error))
cat(sprintf("INTEGRATE_RIGHT %.17g %.17g\n", right_tail$value, right_tail$abs.error))
cat(sprintf("INTEGRATE_LEFT %.17g %.17g\n", left_tail$value, left_tail$abs.error))
cat(sprintf("INTEGRATE_WHOLE %.17g %.17g\n", whole_gaussian$value, whole_gaussian$abs.error))
