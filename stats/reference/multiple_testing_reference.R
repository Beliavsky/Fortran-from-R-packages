# Deterministic references, using stats::p.adjust.
options(digits = 17)
p <- c(0.04, 0.001, 0.03, 0.2, 0.03, 0.8)
for (method in c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")) {
    cat(method, "default:", p.adjust(p, method), "\n")
    cat(method, "n=10:", p.adjust(p, method, n = 10), "\n")
}
cat("missing:", p.adjust(c(0.01, NA, 0.04), "BH"), "\n")
cat("endpoints:", p.adjust(c(0, 1, 0.05), "holm"), "\n")
cat("hommel distinguishing:", p.adjust(c(0.01, 0.04, 0.03, 0.06), "hommel"), "\n")
cat("hommel missing:", p.adjust(c(0.01, NA, 0.04), "hommel", n = 5), "\n")
cat("hommel pair:", p.adjust(c(0.04, 0.01), "hommel"), "\n")
