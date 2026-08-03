generate_random_returns <- function(
  portfolio,
  n,
  seed = NULL
) {

  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  n_assets <- NROW(portfolio)

  standard_deviations <- portfolio$standard_deviation

  which_sd_is_0 <- which(standard_deviations == 0)
  standard_deviations[which_sd_is_0] <- 0.0000001

  covariance_matrix <- 
    calc_covariance_matrix(
      correlations        = portfolio$correlations,
      standard_deviations = standard_deviations
    )

  cholesky_decomposition <- chol(covariance_matrix)

  independent_normal_draws <- 
    matrix(
      stats::rnorm(n_assets * n), 
      nrow = n, 
      ncol = n_assets
    )

  correlated_normal_draws <-
    independent_normal_draws %*% t(cholesky_decomposition)

  simulated_returns <- 
    matrix(
      rep(portfolio$expected_return, each = n), 
      nrow = n, 
      ncol = n_assets
    ) + correlated_normal_draws

  colnames(simulated_returns) <- portfolio$name
  simulated_returns[, which_sd_is_0] <- portfolio$expected_return[which_sd_is_0]
  
  simulated_returns |> 
    dplyr::as_tibble()
}
