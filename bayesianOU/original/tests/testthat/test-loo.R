# Regression test for the PSIS-LOO reshape bug: a 3-D array [draws, time, sector]
# was being passed to loo::loo(), which reads it as [iter, chain, obs] and only
# counts S observations. The helper must reshape to a proper [draws x obs]
# matrix over the window 2:T_lik, i.e. (T_lik - 1) * S observations.

test_that(".loo_from_loglik builds a loo over (T_lik-1)*S observations", {
  skip_if_not_installed("loo")

  set.seed(1)
  draws <- 200L; T <- 12L; S <- 3L; n_chains <- 2L
  vars <- as.vector(vapply(seq_len(S),
                           function(s) sprintf("log_lik[%d,%d]", seq_len(T), s),
                           character(T)))
  ll_mat <- matrix(stats::rnorm(draws * length(vars), -1, 0.5),
                   nrow = draws, dimnames = list(NULL, vars))
  chain_id <- rep(seq_len(n_chains), each = draws / n_chains)

  T_lik <- 9L
  res <- .loo_from_loglik(ll_mat, chain_id, T_lik, S)

  expect_s3_class(res, "loo")
  # The fitted window is 2:T_lik -> (T_lik - 1) observations per sector.
  expect_equal(nrow(res$pointwise), (T_lik - 1L) * S)
})

test_that(".loo_from_loglik errors when expected columns are missing", {
  skip_if_not_installed("loo")
  draws <- 50L; T <- 5L; S <- 2L
  vars <- as.vector(vapply(seq_len(S),
                           function(s) sprintf("log_lik[%d,%d]", seq_len(T), s),
                           character(T)))
  ll_mat <- matrix(stats::rnorm(draws * length(vars)),
                   nrow = draws, dimnames = list(NULL, vars))
  chain_id <- rep(1L, draws)
  # Ask for a window beyond what the matrix contains.
  expect_error(.loo_from_loglik(ll_mat, chain_id, T_lik = 8L, S = S),
               "Missing")
})
