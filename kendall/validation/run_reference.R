# SPDX-License-Identifier: GPL-2.0-or-later
# Run from the Kendall package directory. Requires R, Kendall 2.2.2 and gfortran.
# Compiles only maintained sources and the sibling rfortran-core into a temporary directory.
main <- function() {
stopifnot(as.character(packageVersion("Kendall")) == "2.2.2")
root <- normalizePath(".", mustWork = TRUE)
work <- tempfile("kendall-reference-")
dir.create(work)
on.exit(unlink(work, recursive = TRUE), add = TRUE)

compiler <- Sys.getenv("FC", "gfortran")
core <- normalizePath("../rfortran-core/src", mustWork = TRUE)
sources <- file.path(core, c("r_kinds.f90", "r_special.f90", "r_transforms.f90",
                           "r_ordering.f90", "r_distributions.f90"))
sources <- c(sources, file.path(root, "src/kendall_module.f90"),
             file.path(root, "validation/reference_driver.f90"))
exe <- file.path(work, if (.Platform$OS.type == "windows") "reference.exe" else "reference")
old <- setwd(work)
code <- system2(compiler, c("-std=f2008", "-O2", "-J", shQuote(work), "-I", shQuote(work),
                           shQuote(sources), "-o", shQuote(exe)))
setwd(old)
stopifnot(code == 0)

run_fortran <- function(x, y = x, operation = 1L, period = 1L, method = 1L) {
    input <- file.path(work, "input.txt")
    con <- file(input, "w")
    writeLines(paste(operation, length(x), period, method), con)
    write.table(cbind(x, y), con, row.names = FALSE, col.names = FALSE, quote = FALSE,
                na = "NaN", sep = " ")
    close(con)
    output <- system2(exe, shQuote(input), stdout = TRUE)
    result <- scan(text = output, quiet = TRUE)
    names(result) <- c("status", "n_used", "exact", "tau", "sl", "S", "D", "varS")
    result
}

# Direct double-precision implementation for an independent R reference to ordinary results.
reference <- function(x, y) {
    keep <- !is.na(x) & !is.na(y)
    x <- x[keep]
    y <- y[keep]
    n <- length(x)
    score <- sum(sign(outer(x, x, "-")) * sign(outer(y, y, "-"))) / 2
    # outer-product signs above always yield symmetric products, each pair appearing twice.
    tx <- as.numeric(table(x))
    ty <- as.numeric(table(y))
    cn <- n * (n - 1)
    vx <- sum(tx * (tx - 1) * (2 * tx + 5))
    vy <- sum(ty * (ty - 1) * (2 * ty + 5))
    px <- sum(tx * (tx - 1))
    py <- sum(ty * (ty - 1))
    bx <- sum(tx * (tx - 1) * (tx - 2))
    by <- sum(ty * (ty - 1) * (ty - 2))
    d <- sqrt((cn - px) * (cn - py)) / 2
    varS <- (cn * (2 * n + 5) - vx - vy) / 18 + bx * by / (9 * cn * (n - 2)) + px * py / (2 * cn)
    if (d == 0) return(c(tau = NaN, sl = NaN, S = 0, D = 0, varS = 0))
    c(tau = score / d, sl = 2 * pnorm(-max(0, abs(score) - 1) / sqrt(varS)),
      S = score, D = d, varS = varS)
}
check <- function(actual, expected, atol = 3e-13, rtol = 3e-13) {
    stopifnot(all(is.finite(actual)), all(abs(actual - expected) <= atol + rtol * abs(expected)))
}
quiet_upstream <- function(expression) {
    capture.output(value <- expression)
    unlist(value[c("tau", "sl", "S", "D", "varS")], use.names = TRUE)
}
fields <- c("tau", "sl", "S", "D", "varS")
set.seed(271828)
max_upstream_scaled_error <- 0
cases <- 0L
for (n in rep(c(4, 5, 6, 7, 8, 10, 17, 25, 50, 100), each = 12)) {
    x <- sample(-20:20, n, replace = TRUE) / 2
    y <- sample(-20:20, n, replace = TRUE) / 2
    if (n >= 8) {
        x[2] <- NA_real_
        y[n - 1] <- NA_real_
    }
    actual <- run_fortran(x, y)
    stopifnot(actual["status"] == 0)
    check(actual[fields], reference(x, y))
    upstream <- quiet_upstream(Kendall::Kendall(x, y))
    check(actual[fields], upstream, atol = 4e-7, rtol = 2e-6)
    max_upstream_scaled_error <- max(max_upstream_scaled_error,
                                    max(abs(actual[fields] - upstream) / pmax(1, abs(upstream))))
    trend <- run_fortran(y, operation = 2L)
    check(trend[fields], reference(seq_along(y), y))
    check(trend[fields], quiet_upstream(Kendall::MannKendall(y)), atol = 4e-7, rtol = 2e-6)
    cases <- cases + 2L
}

exact_cases <- 0L
for (n in c(3:30, 40, 50, 75, 100)) {
    x <- seq_len(n)
    for (y in list(sample(x), x, rev(x))) {
        actual <- run_fortran(x, y, method = 2L)
        stopifnot(actual["status"] == 0, actual["exact"] == 1)
        if (abs(actual["S"]) == n * (n - 1) / 2) {
            expected <- 2 / exp(lgamma(n + 1))
            check(actual["sl"], expected, atol = 0, rtol = 3e-13)
        } else {
            expected <- cor.test(x, y, method = "kendall", exact = TRUE)$p.value
            check(actual["sl"], expected, atol = 2e-14, rtol = 2e-12)
        }
        exact_cases <- exact_cases + 1L
    }
}

seasonal_cases <- 0L
for (period in c(1, 2, 3, 4, 12)) {
    for (cycles in c(5, 8, 12)) {
        for (offset in c(0, 1)) {
            n <- cycles * period + offset
            x <- sample(-30:30, n, replace = TRUE)
            x[period + 1] <- NA_real_
            actual <- run_fortran(x, operation = 3L, period = period)
            stopifnot(actual["status"] == 0)
            parts <- lapply(seq_len(period), function(s) {
                v <- x[seq(s, n, by = period)]
                reference(seq_along(v), v)
            })
            sums <- Reduce(`+`, lapply(parts, function(p) p[c("S", "D", "varS")]))
            expected <- c(tau = unname(sums["S"] / sums["D"]),
                          sl = 2 * pnorm(-abs(sums["S"]) / sqrt(sums["varS"])), sums)
            names(expected) <- fields
            check(actual[fields], expected)
            upstream <- quiet_upstream(Kendall::SeasonalMannKendall(ts(x, frequency = period, start = c(2000, period))))
            check(actual[fields], upstream, atol = 4e-7, rtol = 2e-6)
            seasonal_cases <- seasonal_cases + 1L
        }
    }
}
# Mixed constant seasons are valid: only their zero score, denominator and variance are used.
x <- as.vector(rbind(1:6, rep(5, 6), c(8, 6, 7, 3, 4, 2)))
actual <- run_fortran(x, operation = 3L, period = 3L)
check(actual[fields], quiet_upstream(Kendall::SeasonalMannKendall(ts(x, frequency = 3))),
      atol = 4e-7, rtol = 2e-6)
seasonal_cases <- seasonal_cases + 1L

# Deliberate differences: upstream's exact branch is unreachable, short samples are not tested,
# single-precision conversion collapses close values, and its lower-tail subtraction loses tiny p-values.
for (n in c(3, 4, 5, 100)) {
    actual <- run_fortran(seq_len(n))
    exact <- run_fortran(seq_len(n), method = 2L)
    upstream <- quiet_upstream(Kendall::Kendall(seq_len(n), seq_len(n)))
    cat(sprintf("monotonic n=%d: upstream p=%.17g, dp normal p=%.17g, exact p=%.17g\n",
                n, upstream["sl"], actual["sl"], exact["sl"]))
}
x <- 1 + seq_len(5) * 1e-10
actual <- run_fortran(x, operation = 2L)
upstream <- quiet_upstream(Kendall::MannKendall(x))
stopifnot(actual["status"] == 0, actual["tau"] == 1, upstream["D"] == 0)
cat(sprintf("close-value n=5: upstream D=%.17g, dp D=%.17g\n", upstream["D"], actual["D"]))
cat(sprintf("PASS: %d Kendall/Mann-Kendall cases, %d exact comparisons, %d seasonal cases.\n",
            cases, exact_cases, seasonal_cases))
cat(sprintf("Maximum upstream-scaled error in ordinary pairwise cases: %.6g\n", max_upstream_scaled_error))
cat(sprintf("R %s; Kendall %s; compiler: %s\n", getRversion(), packageVersion("Kendall"), compiler))
}
main()
