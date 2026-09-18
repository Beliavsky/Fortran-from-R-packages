# Regenerate deterministic weighted and rank-deficient linear-model fixtures.
options(digits = 17)

x <- c(0, 1, 2, 3, 4, 5)
z <- c(1, 0, 1, 0, 1, 0)
y <- c(-1.8, 4.9, 4.4, 10.2, 10.1, 15.3)
w <- c(1, 2, 1, 3, 2, 1)
fit <- lm.wfit(cbind(1, x, z), y, w)
print(fit$coefficients)
print(fit$fitted.values)
print(fit$residuals)
print(c(rank = fit$rank, df = sum(w > 0) - fit$rank,
        rss = sum(w * fit$residuals^2)))
print(solve(crossprod(sqrt(w) * cbind(1, x, z))))

rank_fit <- lm.fit(cbind(1, x, 2 * x), y)
print(rank_fit$coefficients)
print(rank_fit$fitted.values)
print(c(rank = rank_fit$rank, df = length(y) - rank_fit$rank,
        rss = sum(rank_fit$residuals^2)))

group <- factor(c(1, 1, 2, 2, 3, 3))
one_way <- summary(aov(y ~ group))[[1]]
print(one_way)
print(c(ss_between = one_way[1, "Sum Sq"], ss_within = one_way[2, "Sum Sq"],
        f = one_way[1, "F value"], p = one_way[1, "Pr(>F)"]))

reduced <- lm(y ~ x)
full <- lm(y ~ x + z)
comparison <- anova(reduced, full)
print(comparison)
print(c(ss = comparison[2, "Sum of Sq"], f = comparison[2, "F"],
        p = comparison[2, "Pr(>F)"]))
