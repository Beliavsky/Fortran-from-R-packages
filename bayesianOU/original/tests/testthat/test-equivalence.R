# Permanent bit-exact equivalence guard for the single-level mode.
#
# The unified ou_nested.stan subsumed and retired the legacy ou_nonlinear_tmg.stan.
# Before retiring it we certified that ou_nested.stan with n_levels = 1 reproduces
# the legacy log_lik bit-for-bit, recomputing both via generate_quantities over the
# SAME input draws (so the comparison reflects model-code equivalence, not the
# CSV-serialization rounding that contaminates the in-sampling log_lik). The legacy
# draws, the nested data and the reference log_lik are frozen in fixtures/, so this
# guard stays runnable after the legacy .stan is gone.
#
# Needs cmdstanr + a one-off model compilation; skipped on CRAN and when no Stan
# backend is available. The compiled binary is cached across runs.

test_that("ou_nested.stan (n_levels=1) reproduces the legacy log_lik bit-for-bit", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")
  skip_if_not_installed("posterior")
  skip_if(check_stan_backend() == "none", "No Stan backend available.")

  draws  <- readRDS(test_path("fixtures", "equiv_legacy_draws.rds"))
  data_n <- readRDS(test_path("fixtures", "equiv_nested_data.rds"))
  golden <- readRDS(test_path("fixtures", "equiv_golden_loglik.rds"))

  cache <- tryCatch({
    d <- tools::R_user_dir("bayesianOU", "cache")
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    if (dir.exists(d)) d else tempdir()
  }, error = function(e) tempdir())

  mod <- cmdstanr::cmdstan_model(.stan_file_path(), dir = cache,
                                 cpp_options = list(stan_threads = TRUE))
  gq <- mod$generate_quantities(fitted_params = draws, data = data_n,
                                threads_per_chain = 1)

  ll <- gq$draws("log_lik", format = "matrix")
  ll <- ll[, colnames(golden), drop = FALSE]

  # Bit-for-bit: identical numeric values, not merely all.equal within tolerance.
  expect_identical(dim(ll), dim(golden))
  expect_identical(unname(ll), unname(golden))
})
