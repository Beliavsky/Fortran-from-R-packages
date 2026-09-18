# Deterministic R 4.6.0 fixtures for additional hypothesis tests.
binomial <- binom.test(7, 20, p = 0.2, conf.level = 0.9)
poisson <- poisson.test(12, T = 4, r = 2, conf.level = 0.9)
x <- c(4.2, 5.1, 3.8, 6.0, 5.5, 4.7)
y <- c(2.1, 2.8, 3.5, 2.9, 3.2, 2.6, 3.0)
variance <- var.test(x, y, ratio = 1.25, conf.level = 0.9)
values <- c(4, 5, 6, 5, 7, 8, 3, 4, 5, 6, 6, 7, 8, 9, 11)
groups <- rep(1:3, each = 5)
bartlett <- bartlett.test(values, groups)
blocks <- matrix(c(8, 7, 9, 6,
                   5, 5, 7, 4,
                   9, 8, 10, 7,
                   6, 7, 8, 5,
                   7, 6, 8, 6), nrow = 5, byrow = TRUE)
friedman <- friedman.test(blocks)

cat(sprintf("BINOM %.17g %.17g %.17g %.17g\n", binomial$p.value,
            binomial$estimate, binomial$conf.int[1], binomial$conf.int[2]))
cat(sprintf("POISSON %.17g %.17g %.17g %.17g\n", poisson$p.value,
            poisson$estimate, poisson$conf.int[1], poisson$conf.int[2]))
cat(sprintf("VAR %.17g %.17g %.17g %.17g\n", variance$statistic,
            variance$p.value, variance$conf.int[1], variance$conf.int[2]))
cat(sprintf("BARTLETT %.17g %d %.17g\n", bartlett$statistic,
            bartlett$parameter, bartlett$p.value))
cat(sprintf("FRIEDMAN %.17g %d %.17g\n", friedman$statistic,
            friedman$parameter, friedman$p.value))
