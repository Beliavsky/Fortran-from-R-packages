toy_residuals <- function(T = 30, K = 3) {
  # simple AR(0) gaussian noise with small corr structure
  set.seed(42)
  Z <- matrix(rnorm(T * K), T, K)
  colnames(Z) <- paste0("A", seq_len(K))
  Z
}

toy_cor_pairs <- function(K) combn(K, 2, simplify = FALSE)

toy_sigma_matrix <- function(T = 30, K = 3) {
  # pretend these are forecasted vols (positive)
  S <- matrix(0.2 + 0.05*abs(sin(seq_len(T)/7)), T, K)
  colnames(S) <- paste0("A", seq_len(K))
  S
}

# Simple 2-regime correlation parameterization for K=3
toy_rho_matrix <- function(K = 3) {
  # lower-tri order: (2,1), (3,1), (3,2)
  rbind(
    c(0.10, 0.05, 0.00),
    c(0.60, 0.40, 0.30)
  )
}

toy_P <- function() {
  # N=2, moderately persistent
  matrix(c(0.9, 0.1,
           0.2, 0.8), nrow = 2, byrow = TRUE)
}
