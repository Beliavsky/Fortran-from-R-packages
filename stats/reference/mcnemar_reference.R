# Deterministic matrix-interface references for stats::mcnemar.test.
options(digits = 17)
tables <- list(binary = matrix(c(20, 12, 5, 30), 2),
               equal = matrix(c(20, 7, 7, 30), 2),
               fractional = matrix(c(20, 1.25, 1.5, 30), 2),
               three = matrix(c(10, 4, 2, 7, 12, 6, 3, 9, 8), 3),
               empty_pair = diag(c(10, 20)))
for (name in names(tables)) {
    for (correct in c(TRUE, FALSE)) {
        test <- mcnemar.test(tables[[name]], correct = correct)
        cat(name, correct, test$statistic, test$parameter, test$p.value, "\n")
    }
}
