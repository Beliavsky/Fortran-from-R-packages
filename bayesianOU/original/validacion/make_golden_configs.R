# Regenerate the per-configuration log_lik regression goldens for ou_nested.stan.
#
# Each golden freezes (i) a small set of posterior draws of the CURRENT model,
# (ii) the exact stan data (with that configuration's level_spec flags) and
# (iii) the reference log_lik recomputed via generate_quantities over those
# draws. test-golden-configs.R recomputes log_lik on the frozen draws + data and
# asserts bit-for-bit equality, catching any change to a configuration's Stan
# code path (decision D-IMPL-1).
#
# The DATA contract is reused from the existing fixtures and EXTENDED with the
# fields a later Stan revision added (only the parameter block changed when the
# Level-2 latent path moved to the centered parametrization, D-IMPL-9.1; the data
# contract gained the measurement-SD latency dial sigma_phi_meas_fixed /
# sigma_phi_meas_value, D-IMPL-9.4). The draws are regenerated fresh under the
# current model. `single` (n_levels = 1) is unchanged by these surgeries but is
# regenerated here too so all goldens come from one current script.
#
# A sixth golden, `fixed_meas`, exercises the K-deterministic fixed-SD code path
# (sigma_phi_meas_fixed = 1, no sigma_phi_meas parameter), cloned from the
# canonical 2-level data contract. This guards the new branch added in D-IMPL-9.4.
#
# Run with NOT_CRAN=true from any dir:
#   NOT_CRAN=true Rscript validacion/make_golden_configs.R
suppressWarnings(suppressMessages({
  GB <- "/mnt/kingston/Carpeta de Estudio/[Teoría Marxista]/6. [Mis Investigaciones]/LIBRERÍAS EN R/CARPETA bayesianOU/bayesianOU"
  pkgload::load_all(GB, quiet = TRUE)
  library(cmdstanr); library(posterior)
}))

cache <- tryCatch({
  d <- tools::R_user_dir("bayesianOU", "cache")
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  if (dir.exists(d)) d else tempdir()
}, error = function(e) tempdir())

mod <- cmdstan_model(.stan_file_path(), dir = cache,
                     cpp_options = list(stan_threads = TRUE))

fixdir <- file.path(GB, "tests/testthat/fixtures")
configs <- c("single", "canonical", "both_full", "both_lean", "n1_lean")

# Ensure the data-contract fields a later Stan revision added exist on a reused
# fixture: the measurement-SD latency dial (D-IMPL-9.4) and the Level-3 value
# anchor V_anchor_z (D-IMPL-10, length 0 unless n_levels == 3). Defaults to the
# estimated SD mode (flag 0) and no value anchor; fixed_meas / level3 override.
.with_new_fields <- function(d, fixed = 0L, value = 0.5, V = NULL) {
  d$sigma_phi_meas_fixed <- as.integer(fixed)
  d$sigma_phi_meas_value <- as.numeric(value)
  d$V_anchor_z <- if (is.null(V)) matrix(0.0, 0L, 0L) else as.matrix(V)
  d
}

fit_and_save <- function(cfg, data_cfg, fx) {
  cat(sprintf("[%s] fitting (n_levels=%s, sigma_phi_meas_fixed=%s) ...\n",
              cfg, data_cfg$n_levels, data_cfg$sigma_phi_meas_fixed))
  fit <- mod$sample(
    data = data_cfg, chains = 1L, iter_warmup = 400L, iter_sampling = 60L,
    threads_per_chain = 1L, refresh = 0, seed = 20260610L,
    show_messages = FALSE, show_exceptions = FALSE,
    init = if (data_cfg$n_levels >= 2L) 0.3 else 2)
  draws <- fit$draws()
  gq <- mod$generate_quantities(fitted_params = draws, data = data_cfg,
                                threads_per_chain = 1L)
  loglik <- gq$draws("log_lik", format = "matrix")
  saveRDS(list(draws = draws, data = data_cfg, loglik = loglik), fx)
  cat(sprintf("    saved %s  (draws nvar=%d, loglik %s)\n",
              basename(fx), dim(draws)[3], paste(dim(loglik), collapse = "x")))
}

canonical_data <- NULL
for (cfg in configs) {
  fx  <- file.path(fixdir, sprintf("golden_%s.rds", cfg))
  old <- readRDS(fx)
  data_cfg <- .with_new_fields(old$data)       # reuse the data contract + extend
  if (cfg == "canonical") canonical_data <- data_cfg
  fit_and_save(cfg, data_cfg, fx)
}

# Sixth golden: K-deterministic fixed-SD path (no sigma_phi_meas parameter),
# cloned from the canonical 2-level data with the dial fixed at 0.05 (D-IMPL-9.4).
data_fixed <- .with_new_fields(canonical_data, fixed = 1L, value = 0.05)
fit_and_save("fixed_meas", data_fixed,
             file.path(fixdir, "golden_fixed_meas.rds"))

# Seventh golden: Level-3 (values) path (D-IMPL-10). Clones the canonical 2-level
# data, sets n_levels = 3 and a synthetic standardized value anchor V (T x S), so
# the m_v * V coupling code path is exercised and bit-exact-guarded.
{
  Tn <- canonical_data$T; S <- canonical_data$S
  set.seed(20260610L)
  Vz <- matrix(stats::rnorm(Tn * S), Tn, S)
  Vz <- scale(Vz)                              # standardized, like the R side does
  data_l3 <- .with_new_fields(canonical_data, V = Vz)
  data_l3$n_levels <- 3L
  fit_and_save("level3", data_l3, file.path(fixdir, "golden_level3.rds"))
}

cat("LISTO make_golden_configs.\n")
