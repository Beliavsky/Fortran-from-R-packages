test_that("rsdc_forecast returns expected shapes (const)", {
  K <- 3; T <- 25
  y <- toy_residuals(T, K)
  S <- toy_sigma_matrix(T, K)

  out <- rsdc_forecast(method = "const",
                       N = 1,
                       residuals = y,
                       X = NULL,
                       final_params = list(
                         correlations = matrix(c(0.2, 0.1, 0.05), nrow = 1),
                         log_likelihood = -100.0
                       ),
                       sigma_matrix = S,
                       value_cols = colnames(S),
                       out_of_sample = FALSE)

  expect_equal(dim(out$predicted_correlations), c(T, K*(K-1)/2))
  expect_equal(length(out$cov_matrices), T)
  expect_true(is.finite(out$BIC))
})

test_that("rsdc_forecast 'const' uses final_params$correlations, not empirical cor", {
  K <- 3; T <- 20
  y <- toy_residuals(T, K)
  S <- toy_sigma_matrix(T, K)
  rho_fixed <- matrix(c(0.9, 0.9, 0.9), nrow = 1)
  out <- rsdc_forecast(method = "const", N = 1, residuals = y,
                       final_params = list(correlations = rho_fixed, log_likelihood = -100),
                       sigma_matrix = S, value_cols = colnames(S))
  expect_true(all(abs(out$predicted_correlations - 0.9) < 1e-10))
})

test_that("rsdc_forecast returns expected shapes (noX and tvtp)", {
  K <- 3; T <- 25; N <- 2
  y <- toy_residuals(T, K)
  S <- toy_sigma_matrix(T, K)
  rho <- toy_rho_matrix(K)

  # noX branch (uses provided P through final_params)
  P <- toy_P()
  out_noX <- rsdc_forecast(method = "noX",
                           N = N,
                           residuals = y,
                           X = NULL,
                           final_params = list(
                             correlations = rho,
                             transition_matrix = P,
                             log_likelihood = -200
                           ),
                           sigma_matrix = S,
                           value_cols = colnames(S),
                           out_of_sample = FALSE)

  expect_equal(dim(out_noX$predicted_correlations), c(T, K*(K-1)/2))
  expect_equal(dim(out_noX$smoothed_probs), c(N, T))
  expect_true(is.finite(out_noX$BIC))

  # tvtp branch (uses X and beta)
  X <- cbind(1, scale(seq_len(T)))
  beta <- rbind(c(1.2, 0.0), c(0.8, -0.1))

  out_tvtp <- rsdc_forecast(method = "tvtp",
                            N = N,
                            residuals = y,
                            X = X,
                            final_params = list(
                              correlations = rho,
                              beta = beta,
                              log_likelihood = -210
                            ),
                            sigma_matrix = S,
                            value_cols = colnames(S),
                            out_of_sample = FALSE)

  expect_equal(dim(out_tvtp$predicted_correlations), c(T, K*(K-1)/2))
  expect_equal(dim(out_tvtp$smoothed_probs), c(N, T))
  expect_true(is.finite(out_tvtp$BIC))
})

test_that("rsdc_forecast out_of_sample=TRUE shapes match OOS horizon (noX)", {
  K <- 3; T <- 30; N <- 2
  y <- toy_residuals(T, K)
  S <- toy_sigma_matrix(T, K)
  rho <- toy_rho_matrix(K)
  P   <- toy_P()
  n_oos <- T - round(0.7 * T)

  out <- rsdc_forecast(method = "noX", N = N,
                       residuals = y, X = NULL,
                       final_params = list(correlations = rho,
                                           transition_matrix = P,
                                           log_likelihood = -200,
                                           log_likelihood_oos = -80),
                       sigma_matrix = S, value_cols = colnames(S),
                       out_of_sample = TRUE)

  expect_equal(nrow(out$predicted_correlations), n_oos)
  expect_equal(nrow(out$sigma_matrix), n_oos)
  expect_equal(length(out$cov_matrices), n_oos)
  expect_equal(dim(out$smoothed_probs), c(N, n_oos))
  expect_true(is.finite(out$BIC))
})

test_that("rsdc_forecast out_of_sample=TRUE shapes match OOS horizon (tvtp)", {
  K <- 3; T <- 30; N <- 2
  y <- toy_residuals(T, K)
  S <- toy_sigma_matrix(T, K)
  rho <- toy_rho_matrix(K)
  X   <- cbind(1, scale(seq_len(T)))
  beta <- rbind(c(1.2, 0.0), c(0.8, -0.1))
  n_oos <- T - round(0.7 * T)

  out <- rsdc_forecast(method = "tvtp", N = N,
                       residuals = y, X = X,
                       final_params = list(correlations = rho,
                                           beta = beta,
                                           log_likelihood = -210,
                                           log_likelihood_oos = -85),
                       sigma_matrix = S, value_cols = colnames(S),
                       out_of_sample = TRUE)

  expect_equal(nrow(out$predicted_correlations), n_oos)
  expect_equal(nrow(out$sigma_matrix), n_oos)
  expect_equal(length(out$cov_matrices), n_oos)
  expect_equal(dim(out$smoothed_probs), c(N, n_oos))
  expect_true(is.finite(out$BIC))
})

test_that("rsdc_forecast BIC is a scalar NA (not numeric(0)) when log_likelihood is missing", {
  K <- 3; T <- 20
  y <- toy_residuals(T, K)
  S <- toy_sigma_matrix(T, K)
  expect_warning(
    out <- rsdc_forecast(method = "const", N = 1, residuals = y,
                         final_params = list(
                           correlations = matrix(c(0.2, 0.1, 0.05), nrow = 1)),
                         sigma_matrix = S, value_cols = colnames(S)),
    "BIC is NA")
  expect_length(out$BIC, 1L)
  expect_true(is.na(out$BIC))
})
