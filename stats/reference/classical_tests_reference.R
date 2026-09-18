# Generate deterministic reference values for test/test_classical_tests.f90.
print_result <- function(name, result) {
  cat(name, "\n")
  cat(sprintf("stat=%.17g df=%s p=%.17g\n", result$statistic,
              if (is.null(result$parameter)) "NA" else sprintf("%.17g", result$parameter),
              result$p.value))
  if (!is.null(result$estimate)) {
    cat(sprintf("estimate=%s\n", paste(sprintf("%.17g", result$estimate), collapse = ",")))
  }
  if (!is.null(result$conf.int)) {
    cat(sprintf("ci=%.17g,%.17g\n", result$conf.int[1], result$conf.int[2]))
  }
}

print_result("chisq_gof", chisq.test(c(10, 20, 30), p = c(0.2, 0.3, 0.5)))
contingency <- matrix(c(10, 20, 20, 10), nrow = 2)
print_result("chisq_matrix_corrected", chisq.test(contingency))
print_result("chisq_matrix_uncorrected", chisq.test(contingency, correct = FALSE))
print_result("prop_one_corrected", prop.test(60, 100, p = 0.5))
print_result("prop_one_uncorrected", prop.test(60, 100, p = 0.5, correct = FALSE))
print_result("prop_two_corrected", prop.test(c(30, 45), c(50, 80)))
print_result("prop_two_uncorrected", prop.test(c(30, 45), c(50, 80), correct = FALSE))
x <- c(1, 2, 4, 7, 11, 16)
y <- c(2, 1, 5, 8, 9, 15)
print_result("cor_pearson", cor.test(x, y, method = "pearson"))
print_result("cor_spearman_asymptotic", cor.test(x, y, method = "spearman", exact = FALSE))
