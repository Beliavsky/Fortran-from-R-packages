# Random.Start generates a single random orthogonal matrix of size k x k.
#
# References:
#   Stewart, G. W. (1980). The efficient generation of random orthogonal matrices
#     with an application to condition estimators. SIAM Journal on Numerical
#     Analysis, 17(3), 403-409. http://www.jstor.org/stable/2156882
#
#   Mezzadri, F. (2007). How to generate random matrices from the classical
#     compact groups. Notices of the American Mathematical Society, 54(5),
#     592-604. https://arxiv.org/abs/math-ph/0609050
#
# Updated in GPArotation 2024.2-1 following a suggestion by Yves Rosseel.
# The original qr.Q(qr(matrix(rnorm(k*k), k))) does not guarantee uniform
# sampling from the Haar measure; the sign correction on R's diagonal does.

Random.Start <- function(k = 2L) {
  # 1. Generate a k x k matrix of standard normal variates and QR decompose
  Z      <- matrix(rnorm(k * k), nrow = k, ncol = k)
  qr.out <- qr(Z)

  # 2. Extract Q and the signs of R's diagonal
  Q <- qr.Q(qr.out)
  d <- sign(diag(qr.R(qr.out)))

  # 3. Scale columns of Q by the signs - ensures uniform sampling (Haar measure)
  Q %*% diag(d)
}