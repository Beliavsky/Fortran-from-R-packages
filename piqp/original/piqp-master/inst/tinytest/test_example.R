# This file is part of PIQP-R.
#
# Copyright (c) 2023 piqp-r authors
#
# This source code is licensed under the BSD 2-Clause License found in the
# LICENSE file in the root directory of this source tree.

## A simple QP example
TOL <- 1e-6

P <- matrix(c(6, 0, 0, 4), nrow = 2)
c <- c(-1, -4)
A <- matrix(c(1, -2), nrow = 1)
b <- 0
G <- matrix(c(1, 0, -1, 0), nrow = 2)
h_u <- c(1, 1)
x_l <- c(-Inf, -1)
x_u <- c(Inf, 1)

s <- solve_piqp(P, c, A, b, G, h_u = h_u, x_l = x_l, x_u = x_u, backend = "dense")

expect_equal(s$status, 1, tolerance = TOL)
expect_equal(s$x, c(0.4285714, 0.2142857), tolerance = TOL)
expect_equal(s$y, -1.5714286, tolerance = TOL)
expect_equal(s$z_l, c(0, 0), tolerance = TOL)
expect_equal(s$z_bl, c(0, 0), tolerance = TOL)
expect_equal(s$z_bu, c(0, 0), tolerance = TOL)


## Triplet form
if (requireNamespace("slam", quietly = TRUE)) {

  A1 <- slam::as.simple_triplet_matrix(A)
  P1 <- slam::as.simple_triplet_matrix(P)
  G1 <- slam::as.simple_triplet_matrix(G)

  s1 <- solve_piqp(P1, c, A1, b, G1, h_u = h_u, x_l = x_l, x_u = x_u, backend = "sparse")

  expect_equal(s1$status, 1, tolerance = TOL)
  expect_equal(s1$x, c(0.4285714, 0.2142857), tolerance = TOL)
  expect_equal(s1$y, -1.5714286, tolerance = TOL)
  expect_equal(s1$z_l, c(0, 0), tolerance = TOL)
  expect_equal(s1$z_bl, c(0, 0), tolerance = TOL)
  expect_equal(s1$z_bu, c(0, 0), tolerance = TOL)

}

## CSC form
AA <- as(as(A, "generalMatrix"), "CsparseMatrix")
PP <- as(as(P, "generalMatrix"), "CsparseMatrix")
GG <- as(as(G, "generalMatrix"), "CsparseMatrix")

s2 <- solve_piqp(PP, c, AA, b, GG, h_u = h_u, x_l = x_l, x_u = x_u, backend = "sparse")

expect_equal(s2$status, 1, tolerance = TOL)
expect_equal(s2$x, c(0.4285714, 0.2142857), tolerance = TOL)
expect_equal(s2$y, -1.5714286, tolerance = TOL)
expect_equal(s2$z_l, c(0, 0), tolerance = TOL)
expect_equal(s2$z_bl, c(0, 0), tolerance = TOL)
expect_equal(s2$z_bu, c(0, 0), tolerance = TOL)

## Modify P, A, h_u and x_u

P[1, 1] <- 8
A[1, 2] <- -3;
h_u[1] <- 2
x_u[2] <- 2

ss <- solve_piqp(P, c, A, b, G, h_u = h_u, x_l = x_l, x_u = x_u, backend = "dense")

expect_equal(ss$status, 1, tolerance = TOL)
expect_equal(ss$x, c( 0.2763157, 0.0921052), tolerance = TOL)
expect_equal(ss$y, -1.2105263, tolerance = TOL)
expect_equal(ss$z_l, c(0, 0), tolerance = TOL)
expect_equal(ss$z_bl, c(0, 0), tolerance = TOL)
expect_equal(ss$z_bu, c(0, 0), tolerance = TOL)


## Triplet form
if (requireNamespace("slam", quietly = TRUE)) {

  A1 <- slam::as.simple_triplet_matrix(A)
  P1 <- slam::as.simple_triplet_matrix(P)
  G1 <- slam::as.simple_triplet_matrix(G)

  ss2 <- solve_piqp(P1, c, A1, b, G1, h_u = h_u, x_l = x_l, x_u = x_u, backend = "sparse")

  expect_equal(ss2$status, 1, tolerance = TOL)
  expect_equal(ss2$x, c( 0.2763157, 0.0921052), tolerance = TOL)
  expect_equal(ss2$y, -1.2105263, tolerance = TOL)
  expect_equal(ss2$z_l, c(0, 0), tolerance = TOL)
  expect_equal(ss2$z_bl, c(0, 0), tolerance = TOL)
  expect_equal(ss2$z_bu, c(0, 0), tolerance = TOL)

}

## CSC form
AA <- as(as(A, "generalMatrix"), "CsparseMatrix")
PP <- as(as(P, "generalMatrix"), "CsparseMatrix")
GG <- as(as(G, "generalMatrix"), "CsparseMatrix")

ss3 <- solve_piqp(PP, c, AA, b, GG, h_u = h_u, x_l = x_l, x_u = x_u, backend = "sparse")

expect_equal(ss3$status, 1, tolerance = TOL)
expect_equal(ss3$x, c( 0.2763157, 0.0921052), tolerance = TOL)
expect_equal(ss3$y, -1.2105263, tolerance = TOL)
expect_equal(ss3$z_l, c(0, 0), tolerance = TOL)
expect_equal(ss3$z_bl, c(0, 0), tolerance = TOL)
expect_equal(ss3$z_bu, c(0, 0), tolerance = TOL)
