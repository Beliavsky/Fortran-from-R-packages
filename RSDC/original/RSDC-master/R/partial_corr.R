# Unconstrained parameterization of correlation matrices via canonical
# partial correlations (C-vine), following Joe (2006) \doi{10.1016/j.jmva.2005.05.010}.
#
# KEY PROPERTY (Joe 2006, Theorem/Section 2): z is a vector of
# C = K(K-1)/2 numbers, each free in (-1, 1), and the map z <-> R below is a
# bijection with the set of positive-definite correlation matrices. The
# joint positive-definiteness constraint that makes box-constrained global
# search infeasible in the natural parameterization is absorbed into the
# geometry: the whole search box becomes feasible.
#
# Convention: z packs the strict lower triangle in column-major
# (`lower.tri`) order, the same order the estimators pack correlations.
# Entry (j, i) with j > i holds z_{ij} = the partial correlation of
# variables (i, j) given variables 1, ..., i-1; z_{1j} is the plain
# correlation r_{1j}. At K = 2 the map is the identity.

#' Canonical partial correlations to correlation matrix
#'
#' @param zvec Numeric vector of K(K-1)/2 partial correlations in (-1, 1),
#'   strict lower triangle in column-major order.
#' @param K Integer, matrix dimension.
#' @returns A K x K positive-definite correlation matrix.
#' @references Joe (2006) \doi{10.1016/j.jmva.2005.05.010}
#' @noRd
.rsdc_pc_to_corr <- function(zvec, K) {
  Z <- matrix(0, K, K)
  Z[lower.tri(Z)] <- zvec
  R <- diag(K)
  for (j in 2:K) {
    for (i in 1:(j - 1)) {
      # De-partialize z_{ij} by removing conditioning variables i-1, ..., 1:
      # rho_{ij.1..l-1} = rho_{il} rho_{jl} + rho_{ij.1..l} sqrt((1-rho_il^2)(1-rho_jl^2))
      rho <- Z[j, i]
      if (i > 1) for (l in (i - 1):1)
        rho <- Z[i, l] * Z[j, l] + rho * sqrt((1 - Z[i, l]^2) * (1 - Z[j, l]^2))
      R[j, i] <- R[i, j] <- rho
    }
  }
  R
}

#' Correlation matrix to canonical partial correlations (inverse map)
#'
#' @param R A K x K positive-definite correlation matrix.
#' @returns Numeric vector of K(K-1)/2 partial correlations in (-1, 1),
#'   strict lower triangle in column-major order.
#' @references Joe (2006) \doi{10.1016/j.jmva.2005.05.010}
#' @noRd
.rsdc_corr_to_pc <- function(R) {
  K <- ncol(R)
  Z <- matrix(0, K, K)
  S <- R
  # Peel off conditioning variables 1, 2, ... recursively: at stage i,
  # S[a, b] holds rho_{ab.1..i-1} for a, b >= i.
  for (i in 1:(K - 1)) {
    for (j in (i + 1):K) Z[j, i] <- S[j, i]
    if (i < K - 1) {
      S2 <- S
      for (a in (i + 1):K) for (b in (i + 1):K)
        S2[a, b] <- (S[a, b] - S[a, i] * S[b, i]) /
                    sqrt((1 - S[a, i]^2) * (1 - S[b, i]^2))
      S <- S2
    }
  }
  Z[lower.tri(Z)]
}
