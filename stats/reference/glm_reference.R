# Regenerate deterministic reference values embedded in
# test/test_generalized_linear_models.f90 using base R's stats package.

x <- c(-2, -1, 0, 1, 2, 3)
y <- c(0, 0, 1, 0, 1, 1)
fit <- glm(y ~ x, family = binomial())
cat("BIN_COEF", format(coef(fit), digits = 17), "\n")
cat("BIN_FIT", format(fitted(fit), digits = 17), "\n")
cat("BIN_SE", format(coef(summary(fit))[, 2], digits = 17), "\n")
cat("BIN_PEAR", format(residuals(fit, type = "pearson"), digits = 17), "\n")

fit <- glm(y ~ x, family = binomial(link = "probit"))
cat("PROBIT_COEF", format(coef(fit), digits = 17), "\n")
cat("PROBIT_FIT", format(fitted(fit), digits = 17), "\n")

fit <- glm(y ~ x, family = binomial(link = "cloglog"))
cat("CLOGLOG_COEF", format(coef(fit), digits = 17), "\n")
cat("CLOGLOG_FIT", format(fitted(fit), digits = 17), "\n")

x <- 0:5
y <- c(1, 2, 1, 4, 6, 8)
offset <- log(c(1, 1.5, 2, 1, 2.5, 3))
fit <- glm(y ~ x + offset(offset), family = poisson())
cat("POI_COEF", format(coef(fit), digits = 17), "\n")
cat("POI_FIT", format(fitted(fit), digits = 17), "\n")
cat("POI_SE", format(coef(summary(fit))[, 2], digits = 17), "\n")
cat("POI_PEAR", format(residuals(fit, type = "pearson"), digits = 17), "\n")

fit <- glm(y ~ x, family = quasipoisson())
cat("QPOI_COEF", format(coef(fit), digits = 17), "\n")
cat("QPOI_DISP", format(summary(fit)$dispersion, digits = 17), "\n")
cat("QPOI_SE", format(coef(summary(fit))[, 2], digits = 17), "\n")

fit <- glm(y ~ x, family = poisson(link = "sqrt"))
cat("SQRTPOI_COEF", format(coef(fit), digits = 17), "\n")
cat("SQRTPOI_FIT", format(fitted(fit), digits = 17), "\n")

# Weighted Gaussian identity fit with an exactly aliased predictor.
x <- 0:5
y <- c(1.2, 2.0, 2.9, 4.1, 5.2, 5.8)
weights <- c(1, 2, 1, 3, 2, 1)
fit <- glm(y ~ x + I(2 * x), family = gaussian(), weights = weights)
cat("GAU_COEF", format(coef(fit), digits = 17), "\n")
cat("GAU_FIT", format(fitted(fit), digits = 17), "\n")
cat("GAU_DISP", format(summary(fit)$dispersion, digits = 17), "\n")
cat("GAU_DEV", format(deviance(fit), digits = 17), "\n")
cat("GAU_SE", format(coef(summary(fit))[, 2], digits = 17), "\n")

# Gamma regression with the commonly used log link.
x <- c(-1.5, -0.5, 0, 0.5, 1.0, 1.5, 2.0, 2.5)
y <- c(1.1, 1.5, 1.9, 2.2, 3.2, 3.8, 5.4, 6.8)
fit <- glm(y ~ x, family = Gamma(link = "log"))
cat("GAM_COEF", format(coef(fit), digits = 17), "\n")
cat("GAM_FIT", format(fitted(fit), digits = 17), "\n")
cat("GAM_DISP", format(summary(fit)$dispersion, digits = 17), "\n")
cat("GAM_DEV", format(deviance(fit), digits = 17), "\n")
cat("GAM_SE", format(coef(summary(fit))[, 2], digits = 17), "\n")

fit <- glm(y ~ x, family = Gamma())
cat("GAMI_COEF", format(coef(fit), digits = 17), "\n")
cat("GAMI_FIT", format(fitted(fit), digits = 17), "\n")
cat("GAMI_DISP", format(summary(fit)$dispersion, digits = 17), "\n")
cat("GAMI_DEV", format(deviance(fit), digits = 17), "\n")

# Prior-weighted binomial proportions.
x <- c(-2, -1, 0, 1, 2)
y <- c(0.1, 0.2, 0.45, 0.7, 0.9)
weights <- c(10, 15, 20, 15, 10)
fit <- glm(y ~ x, family = binomial(), weights = weights)
cat("WBIN_COEF", format(coef(fit), digits = 17), "\n")
cat("WBIN_FIT", format(fitted(fit), digits = 17), "\n")
cat("WBIN_DEV", format(deviance(fit), digits = 17), "\n")
cat("WBIN_SE", format(coef(summary(fit))[, 2], digits = 17), "\n")

fit <- glm(y ~ x, family = quasibinomial(), weights = weights)
cat("QBIN_COEF", format(coef(fit), digits = 17), "\n")
cat("QBIN_DISP", format(summary(fit)$dispersion, digits = 17), "\n")
cat("QBIN_SE", format(coef(summary(fit))[, 2], digits = 17), "\n")
