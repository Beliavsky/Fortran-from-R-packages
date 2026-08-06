# Canonical partial-correlation parameterization (Joe 2006): bijection with
# the set of positive-definite correlation matrices, in lower.tri packing.

test_that("every z in (-1,1)^C maps to a PD correlation matrix", {
  set.seed(1)
  for (K in c(3L, 5L, 10L)) {
    C <- K * (K - 1L) / 2L
    for (rep in 1:200) {
      R <- .rsdc_pc_to_corr(runif(C, -1, 1), K)
      expect_true(min(eigen(R, symmetric = TRUE, only.values = TRUE)$values) > 0)
      expect_equal(diag(R), rep(1, K))
      expect_identical(R, t(R))
    }
  }
})

test_that("round-trips are exact in both directions", {
  set.seed(2)
  for (K in c(2L, 5L, 10L)) {
    C <- K * (K - 1L) / 2L
    z0 <- runif(C, -0.95, 0.95)
    expect_equal(.rsdc_corr_to_pc(.rsdc_pc_to_corr(z0, K)), z0, tolerance = 1e-12)
    A <- matrix(rnorm(K * K, 0, 0.4), K) + diag(K)
    R0 <- stats::cov2cor(crossprod(A))
    expect_equal(.rsdc_pc_to_corr(.rsdc_corr_to_pc(R0), K), R0, tolerance = 1e-12)
  }
})

test_that("partials match the independent Schur-complement definition", {
  # Partial correlation of (i, j) given 1..i-1, computed from the Schur
  # complement of the conditioning block -- no shared code with the recursion.
  pc_schur <- function(R, i, j) {
    if (i == 1) return(R[i, j])
    S <- seq_len(i - 1)
    Sc <- R[c(i, j), c(i, j)] -
          R[c(i, j), S, drop = FALSE] %*% solve(R[S, S]) %*% R[S, c(i, j), drop = FALSE]
    Sc[1, 2] / sqrt(Sc[1, 1] * Sc[2, 2])
  }
  set.seed(3)
  K <- 6L
  for (rep in 1:20) {
    z0 <- runif(K * (K - 1L) / 2L, -0.8, 0.8)
    R <- .rsdc_pc_to_corr(z0, K)
    Z <- matrix(0, K, K); Z[lower.tri(Z)] <- z0
    for (j in 2:K) for (i in 1:(j - 1))
      expect_equal(pc_schur(R, i, j), Z[j, i], tolerance = 1e-10)
  }
})

test_that("the map is the identity at K = 2", {
  expect_equal(.rsdc_pc_to_corr(0.37, 2L)[2, 1], 0.37)
  R <- matrix(c(1, -0.6, -0.6, 1), 2)
  expect_equal(.rsdc_corr_to_pc(R), -0.6)
})
