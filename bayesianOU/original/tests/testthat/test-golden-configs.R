# Permanent per-configuration log_lik regression guards for ou_nested.stan.
#
# For each experiment configuration (single, canonical, both_full, both_lean,
# n1_lean, fixed_meas) a golden fixture freezes (i) a small set of posterior
# draws, (ii) the
# exact stan data (with that configuration's level_spec flags), and (iii) the
# reference log_lik recomputed via generate_quantities over those draws. This
# test recomputes log_lik with generate_quantities on the frozen draws + data and
# asserts it matches the golden bit-for-bit, catching any change to a
# configuration's Stan code path. Recomputing both via generate_quantities (not
# comparing against an in-sampling log_lik) isolates code equivalence from the
# CSV-serialization rounding (decision D-IMPL-1).
#
# Regenerate the fixtures with validacion/make_golden_configs.R. Needs cmdstanr +
# a one-off compilation; skipped on CRAN and when no Stan backend is available.

test_that("each level_spec configuration reproduces its golden log_lik bit-for-bit", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")
  skip_if_not_installed("posterior")
  skip_if(check_stan_backend() == "none", "No Stan backend available.")

  cache <- tryCatch({
    d <- tools::R_user_dir("bayesianOU", "cache")
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
    if (dir.exists(d)) d else tempdir()
  }, error = function(e) tempdir())

  mod <- cmdstanr::cmdstan_model(.stan_file_path(), dir = cache,
                                 cpp_options = list(stan_threads = TRUE))

  configs <- c("single", "canonical", "both_full", "both_lean", "n1_lean",
               "fixed_meas", "level3")
  for (cfg in configs) {
    fx <- test_path("fixtures", sprintf("golden_%s.rds", cfg))
    skip_if_not(file.exists(fx), sprintf("Missing golden fixture for %s.", cfg))
    g <- readRDS(fx)

    gq <- mod$generate_quantities(fitted_params = g$draws, data = g$data,
                                  threads_per_chain = 1)
    ll <- gq$draws("log_lik", format = "matrix")
    ll <- ll[, colnames(g$loglik), drop = FALSE]

    expect_identical(dim(ll), dim(g$loglik),
                     info = sprintf("config %s: log_lik dimensions", cfg))
    expect_identical(unname(ll), unname(g$loglik),
                     info = sprintf("config %s: log_lik not bit-for-bit", cfg))
  }
})
