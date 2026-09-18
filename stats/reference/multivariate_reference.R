# Generate deterministic reference values for test/test_multivariate.f90.
x <- matrix(c(
  2, 0, 1,
  3, 1, 4,
  4, 1, 2,
  5, 3, 5,
  7, 4, 3,
  8, 6, 7
), ncol = 3, byrow = TRUE)

print_vector <- function(name, value) {
  cat(name, "=", paste(sprintf("%.17g", value), collapse = ","), "\n")
}
print_matrix <- function(name, value) {
  cat(name, "\n")
  for (i in seq_len(nrow(value))) {
    cat(paste(sprintf("%.17g", value[i, ]), collapse = ","), "\n")
  }
}
canonicalize <- function(fit) {
  for (j in seq_len(ncol(fit$rotation))) {
    pivot <- which.max(abs(fit$rotation[, j]))
    if (fit$rotation[pivot, j] < 0) {
      fit$rotation[, j] <- -fit$rotation[, j]
      fit$x[, j] <- -fit$x[, j]
    }
  }
  fit
}

scaled <- scale(x)
print_vector("scale_center", attr(scaled, "scaled:center"))
print_vector("scale_scale", attr(scaled, "scaled:scale"))
print_matrix("scaled_values", scaled)

fit <- canonicalize(prcomp(x, center = TRUE, scale. = FALSE))
print_vector("pca_sdev", fit$sdev)
print_vector("pca_center", fit$center)
print_matrix("pca_rotation", fit$rotation)
print_matrix("pca_scores", fit$x)

center <- colMeans(x)
covariance <- cov(x)
print_matrix("covariance", covariance)
print_matrix("inverse_covariance", solve(covariance))
print_vector("mahalanobis", mahalanobis(x, center, covariance))
print_vector("mahalanobis_inverted", mahalanobis(x, center, solve(covariance), inverted = TRUE))

# Multivariate linear-model and nested MANOVA fixtures.
t <- 0:9
z <- c(0, 1, 0, 1, 0, 1, 0, 1, 0, 1)
responses <- cbind(
  1 + 0.5 * t + 1.2 * z + c(.1, -.1, .2, 0, -.2, .1, 0, -.1, .1, -.1),
  2 - 0.3 * t + 0.8 * z + c(0, .2, -.1, .1, 0, -.2, .1, 0, -.1, .1),
  -1 + 0.2 * t - 0.5 * z + c(.2, 0, -.2, .1, -.1, .2, -.1, 0, .1, -.2)
)
mlm <- lm(responses ~ t + z)
print_matrix("mlm_coefficients", coef(mlm))
print_matrix("mlm_residual_sscp", crossprod(residuals(mlm)))
for (test in c("Pillai", "Wilks", "Hotelling-Lawley", "Roy")) {
  print_vector(paste0("manova_", test),
               unname(summary(manova(responses ~ t + z), test = test)$stats["z", ]))
}

# Canonical-correlation fixture.
cx <- cbind(t, c(2, 1, 4, 3, 7, 5, 8, 6, 10, 9))
cy <- cbind(
  c(1.1, 1.8, 2.9, 3.7, 5.2, 5.8, 7.1, 7.9, 9.2, 9.8),
  c(4, 2, 5, 3, 8, 6, 9, 7, 11, 10)
)
cc <- cancor(cx, cy)
print_vector("cancor_correlations", cc$cor)
print_matrix("cancor_xcoef", cc$xcoef)
print_matrix("cancor_ycoef", cc$ycoef)
print_vector("cancor_xcenter", cc$xcenter)
print_vector("cancor_ycenter", cc$ycenter)
