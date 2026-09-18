# Generate deterministic reference values for test/test_t_tests.f90.
cases <- list(
  one = t.test(c(2.1, 2.5, 2.3, 2.8, 3.0), mu = 2),
  welch = t.test(c(1.2, 2.4, 3.1, 4.8, 5.0), c(0.7, 1.5, 2.0, 2.1)),
  equal = t.test(c(1.2, 2.4, 3.1, 4.8, 5.0), c(0.7, 1.5, 2.0, 2.1), var.equal = TRUE),
  paired = t.test(c(5.1, 6.2, 7.0, 8.4, 9.1), c(4.8, 6.0, 6.5, 8.0, 8.8), paired = TRUE)
)

for (case_name in names(cases)) {
  result <- cases[[case_name]]
  cat(case_name, "\n")
  cat(sprintf("stat=%.17g df=%.17g p=%.17g\n", result$statistic, result$parameter, result$p.value))
  cat(sprintf("estimate=%s\n", paste(sprintf("%.17g", result$estimate), collapse = ",")))
  cat(sprintf("ci=%.17g,%.17g stderr=%.17g\n", result$conf.int[1], result$conf.int[2], result$stderr))
}
