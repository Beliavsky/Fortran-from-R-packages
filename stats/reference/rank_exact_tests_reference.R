# Generate deterministic reference values for test/test_rank_exact_tests.f90.
print_result <- function(name, result) {
  cat(name, "\n")
  if (!is.null(result$statistic)) cat(sprintf("stat=%.17g\n", result$statistic))
  cat(sprintf("p=%.17g\n", result$p.value))
  if (!is.null(result$parameter)) cat(sprintf("df=%.17g\n", result$parameter))
  if (!is.null(result$estimate)) cat(sprintf("estimate=%.17g\n", result$estimate))
}

fisher_table <- matrix(c(1, 11, 9, 3), nrow = 2)
print_result("fisher", fisher.test(fisher_table))
print_result("wilcox_rank_sum", wilcox.test(c(1, 2, 3, 4, 5), c(6, 7, 8, 9),
                                             exact = FALSE, correct = TRUE))
print_result("wilcox_signed_rank", wilcox.test(c(5, 7, 8, 6, 9), c(3, 6, 7, 7, 5),
                                                paired = TRUE, exact = FALSE, correct = TRUE))
values <- c(1, 2, 2, 4, 5, 6, 6, 8, 9)
groups <- factor(rep(c("a", "b", "c"), each = 3))
print_result("kruskal", kruskal.test(values, groups))
print_result("ks_normal_asymptotic", ks.test(c(0.1, 0.4, 0.8, 1.2, 1.5), "pnorm",
                                              mean = 0.5, sd = 0.7, exact = FALSE))
print_result("ks_two_sample_asymptotic", ks.test(c(0.1, 0.4, 0.8, 1.2, 1.5),
                                                  c(0.2, 0.3, 0.9, 1.1, 1.7, 2.0),
                                                  exact = FALSE))
