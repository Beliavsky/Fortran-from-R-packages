# Mock-based tests for the 2-level diagnostics: OOS recursion (market reverts to
# the latent Phi, never to a realized X), the mu_t trajectory extractor, and the
# S3 methods. No Stan backend required.

# ---- helpers ---------------------------------------------------------------

mock_fit_df <- function(cols) {
  # cols: named list of numeric vectors -> a CmdStanMCMC-like object whose
  # draws(p, format="df") returns the columns matching the parameter name.
  df <- as.data.frame(cols, check.names = FALSE)
  structure(list(
    draws = function(vars, format = "df") {
      keep <- unlist(lapply(vars, function(v)
        grep(sprintf("^%s(\\[|$)", v), names(df), value = TRUE)))
      df[, keep, drop = FALSE]
    }
  ), class = "CmdStanMCMC")
}

base_meds <- function(S) {
  list(kappa_tilde = rep(0, S), a3_s = rep(-0.05, S), beta1 = 0.2,
       gamma = 0.1, kappa_p = rep(0.2, S), mu_const = rep(0, S), m1 = 0.5)
}

# ---- evaluate_oos_nested ----------------------------------------------------

test_that("evaluate_oos_nested returns the documented structure", {
  set.seed(11)
  Tn <- 20L; S <- 3L; T_train <- 14L
  Yz <- matrix(stats::rnorm(Tn * S), Tn, S)
  Phi <- matrix(stats::rnorm(Tn * S), Tn, S)
  COM <- matrix(abs(stats::rnorm(Tn * S)) + 0.5, Tn, S)
  K   <- matrix(stats::runif(Tn * S, 1, 5), Tn, S)
  res <- evaluate_oos_nested(base_meds(S), Yz, Phi, stats::rnorm(Tn),
                             stats::rnorm(Tn), T_train, COM_ts = COM, K_ts = K,
                             kappa_cap = 2, com_in_mean = TRUE,
                             horizons = c(1L, 4L, 8L))
  expect_named(res, c("h1", "h4", "h8"))
  expect_true(all(vapply(res, function(o) all(c("h","RMSE","MAE","n_obs") %in% names(o)),
                         logical(1))))
})

test_that("evaluate_oos_nested guards horizons longer than the test window", {
  set.seed(12)
  Tn <- 10L; S <- 2L; T_train <- 7L
  Yz <- matrix(stats::rnorm(Tn * S), Tn, S)
  Phi <- matrix(stats::rnorm(Tn * S), Tn, S)
  COM <- matrix(abs(stats::rnorm(Tn * S)) + 0.5, Tn, S)
  K   <- matrix(1, Tn, S)
  res <- evaluate_oos_nested(base_meds(S), Yz, Phi, stats::rnorm(Tn),
                             stats::rnorm(Tn), T_train, COM_ts = COM, K_ts = K,
                             kappa_cap = 2, horizons = c(1L, 8L))
  expect_true(is.na(res$h8$RMSE)); expect_equal(res$h8$n_obs, 0L)
  expect_true(is.finite(res$h1$RMSE)); expect_gt(res$h1$n_obs, 0L)
})

test_that("zero market speed yields a persistence forecast (phi does not move)", {
  # kappa_tilde very negative -> kappa_m ~ 0 -> phi_pred stays at the origin value,
  # so the h-step error equals Yz[t0+h-1] - Yz[t0-1] regardless of Phi/G'.
  set.seed(13)
  Tn <- 12L; S <- 2L; T_train <- 8L
  Yz <- matrix(stats::rnorm(Tn * S), Tn, S)
  Phi <- matrix(stats::rnorm(Tn * S), Tn, S)   # arbitrary; must not matter here
  COM <- matrix(abs(stats::rnorm(Tn * S)) + 0.5, Tn, S)
  K   <- matrix(1, Tn, S)
  meds <- base_meds(S); meds$kappa_tilde <- rep(-50, S)
  res <- evaluate_oos_nested(meds, Yz, Phi, stats::rnorm(Tn), stats::rnorm(Tn),
                             T_train, COM_ts = COM, K_ts = K, kappa_cap = 2,
                             com_in_mean = FALSE, horizons = 1L)
  origins <- seq.int(T_train + 1, Tn)
  expected_errs <- as.vector(t(Yz[origins, ] - Yz[origins - 1, ]))
  expect_equal(res$h1$RMSE, sqrt(mean(expected_errs^2)), tolerance = 1e-10)
})

test_that("evaluate_oos_nested adds the Level-3 value term to the Phi attractor (D-IMPL-10)", {
  set.seed(31)
  Tn <- 20L; S <- 2L; T_train <- 12L
  Yz  <- matrix(stats::rnorm(Tn * S), Tn, S)
  Phi <- matrix(stats::rnorm(Tn * S), Tn, S)
  COM <- matrix(abs(stats::rnorm(Tn * S)) + 0.5, Tn, S)
  K   <- matrix(stats::runif(Tn * S, 1, 5), Tn, S)
  Gp  <- stats::rnorm(Tn); zt <- stats::rnorm(Tn)
  Vz  <- matrix(stats::rnorm(Tn * S), Tn, S)

  meds0 <- base_meds(S)                  # no value coupling
  meds1 <- c(meds0, list(m_v = 0.8))     # nonzero value coupling

  run <- function(meds, Vz_arg)
    evaluate_oos_nested(meds, Yz, Phi, Gp, zt, T_train, COM_ts = COM, K_ts = K,
                        kappa_cap = 2, com_in_mean = TRUE, horizons = c(2L, 4L),
                        V_z = Vz_arg)

  # (a) Gate on the anchor: m_v present but V_z = NULL reproduces the no-value run.
  expect_equal(run(meds1, NULL), run(meds0, NULL))
  # (b) Gate on the coupling: V_z present but no m_v reproduces the no-value run.
  expect_equal(run(meds0, Vz), run(meds0, NULL))
  # (c) Active when both present: the value term shifts Phi, so the >=2-step
  #     forecast (which consumes the updated Phi) must move.
  expect_false(isTRUE(all.equal(run(meds1, Vz)$h2$RMSE, run(meds1, NULL)$h2$RMSE)))
})

# ---- extract_mu_trajectory --------------------------------------------------

test_that("extract_mu_trajectory reconstructs m0 + m1 * G' with correct dims", {
  set.seed(14)
  S <- 2L; Tn <- 15L; nd <- 200L
  m0 <- c(0.3, -0.4); m1 <- 0.7
  fit <- mock_fit_df(list(
    "mu_const[1]" = rnorm(nd, m0[1], 1e-6),
    "mu_const[2]" = rnorm(nd, m0[2], 1e-6),
    "m1"          = rnorm(nd, m1,   1e-6)
  ))
  Gp <- stats::rnorm(Tn)
  mu <- extract_mu_trajectory(fit, Gp, sector_names = c("A", "B"))
  expect_equal(dim(mu$median), c(Tn, S))
  expect_equal(colnames(mu$median), c("A", "B"))
  expect_equal(mu$median[, 1], m0[1] + m1 * Gp, tolerance = 1e-3)
  expect_equal(mu$median[, 2], m0[2] + m1 * Gp, tolerance = 1e-3)
  expect_true(all(mu$q2.5 <= mu$median + 1e-6))
})

test_that("extract_mu_trajectory returns NULL without Level-2 parameters", {
  fit <- mock_fit_df(list("kappa_s[1]" = rnorm(10)))
  expect_null(extract_mu_trajectory(fit, stats::rnorm(5)))
})

# ---- S3 methods -------------------------------------------------------------

make_mock_2level <- function(Tn = 12L, S = 2L, richness = FALSE) {
  band <- function(v) {
    m <- cbind(q2.5 = v - 0.1, median = v, q97.5 = v + 0.1); m
  }
  mat <- function(v) matrix(v, Tn, S)
  # Level-2 richness extras (D-IMPL-10): present only in the "both_full"-like
  # mock; NULL otherwise, mirroring the canonical configuration.
  rich <- if (richness) list(
    a3_p      = band(c(-0.05, -0.06)),
    nu_p      = matrix(c(5, 7, 11), nrow = 1,
                       dimnames = list(NULL, c("q2.5", "median", "q97.5"))),
    sigma_p_t = list(median = mat(0.3), q2.5 = mat(0.2), q97.5 = mat(0.45))
  ) else list(a3_p = NULL, nu_p = NULL, sigma_p_t = NULL)
  structure(list(
    factor_ou = list(
      model = "ou_nested_2level", n_levels = 2L,
      level2 = c(list(
        kappa_p = band(c(0.20, 0.25)), kappa_m_base = band(c(0.60, 0.65)),
        mu_const = band(c(0.0, 0.1)), sigma_p = band(c(0.3, 0.3)),
        m1 = 0.5, sigma_phi_meas = 0.2), rich),
      factor_ou_info = list(T_train = 8L, fit_window = "train",
                            theta_separation = "soft", k_uncertainty = "meas",
                            kappa_cap = 2)),
    phi_latent = list(median = mat(0.1), q2.5 = mat(-0.2), q97.5 = mat(0.4)),
    mu_path = list(median = mat(0.05), q2.5 = mat(-0.1), q97.5 = mat(0.2),
                   Gprime_z = stats::rnorm(Tn)),
    separation = list(prob_sep_by_sector = c(0.9, 0.85), prob_sep_joint = 0.82),
    diagnostics = list(rhat_max = 1.01, rhat_share = 0.0, divergences = 0,
                       oos = list(h1 = list(h = 1, RMSE = 0.5, MAE = 0.4, n_obs = 6L)))
  ), class = c("ou_nested_2level", "list"))
}

test_that("ou_nested_2level print/summary work and dispatch on the top object", {
  x <- make_mock_2level()
  expect_output(print(x), "ou_nested_2level")
  s <- summary(x)
  expect_s3_class(s, "summary.ou_nested_2level")
  expect_output(print(s), "Level-2 parameters")
  expect_true(any(grepl("kappa_p", s$level2_table$parameter)))
})

test_that("ou_nested_2level plot draws for each type without error", {
  x <- make_mock_2level()
  pf <- tempfile(fileext = ".pdf"); grDevices::pdf(pf)
  on.exit({ grDevices::dev.off(); unlink(pf) }, add = TRUE)
  expect_silent(plot(x, type = "phi", sector = 1))
  expect_silent(plot(x, type = "mu", sector = 2))
  expect_silent(plot(x, type = "separation"))
})

# ---- D-IMPL-10: Level-2 richness exposure (nu_p, sigma_p(t)) ----------------

test_that("summary/plot expose Level-2 richness when present, and not otherwise", {
  # Canonical mock: no a3_p / nu_p rows, sv_p plot errors cleanly.
  x0 <- make_mock_2level(richness = FALSE)
  s0 <- summary(x0)
  expect_false(any(grepl("a3_p|nu_p", s0$level2_table$parameter)))
  expect_null(s0$sigma_p_t)
  pf <- tempfile(fileext = ".pdf"); grDevices::pdf(pf)
  on.exit({ grDevices::dev.off(); unlink(pf) }, add = TRUE)
  expect_error(plot(x0, type = "sv_p"), "Level-2 SV scale")

  # both_full-like mock: a3_p + nu_p rows present, sigma_p(t) printed and plotted.
  x1 <- make_mock_2level(richness = TRUE)
  s1 <- summary(x1)
  expect_true(any(grepl("a3_p", s1$level2_table$parameter)))
  expect_true(any(grepl("nu_p", s1$level2_table$parameter)))
  expect_false(is.null(s1$sigma_p_t))
  expect_output(print(s1), "sigma_p\\(t\\)")
  expect_silent(plot(x1, type = "sv_p", sector = 1))
})

test_that(".extract_level2_summary derives nu_p and sigma_p(t) from draws (D-IMPL-10)", {
  set.seed(21)
  S <- 2L; Tn <- 5L; nd <- 150L
  Hp <- matrix(c(-1.0, -0.5, 0.0, 0.5, 1.0,
                 -0.8, -0.3, 0.1, 0.4, 0.9), Tn, S)
  cols <- list()
  for (s in seq_len(S)) {
    cols[[sprintf("kappa_p[%d]", s)]]      <- rnorm(nd, 0.2, 1e-6)
    cols[[sprintf("kappa_m_base[%d]", s)]] <- rnorm(nd, 0.6, 1e-6)
    cols[[sprintf("mu_const[%d]", s)]]     <- rnorm(nd, 0.0, 1e-6)
    cols[[sprintf("sigma_p[%d]", s)]]      <- rnorm(nd, 0.3, 1e-6)
    for (t in seq_len(Tn))
      cols[[sprintf("h_p[%d,%d]", t, s)]]  <- rnorm(nd, Hp[t, s], 1e-6)
  }
  cols[["m1"]]             <- rnorm(nd, 0.5, 1e-6)
  cols[["sigma_phi_meas"]] <- rnorm(nd, 0.2, 1e-6)
  cols[["nu_p_tilde"]]     <- rnorm(nd, 6,   1e-6)   # nu_p = 2 + 6 = 8
  fit <- mock_fit_df(cols)

  l2 <- bayesianOU:::.extract_level2_summary(fit, S, Tn, c("A", "B"))
  expect_false(is.null(l2$nu_p))
  expect_equal(unname(l2$nu_p[, "median"]), 8, tolerance = 1e-3)
  expect_false(is.null(l2$sigma_p_t))
  expect_equal(dim(l2$sigma_p_t$median), c(Tn, S))
  expect_equal(colnames(l2$sigma_p_t$median), c("A", "B"))
  # median sigma_p(t) ~ exp(h_p / 2) (near-degenerate draws)
  expect_equal(l2$sigma_p_t$median[, 1], exp(0.5 * Hp[, 1]), tolerance = 1e-3)
  expect_equal(l2$sigma_p_t$median[, 2], exp(0.5 * Hp[, 2]), tolerance = 1e-3)
})

test_that(".extract_level2_summary keeps richness NULL in the canonical L2", {
  S <- 2L; Tn <- 5L; nd <- 80L
  cols <- list()
  for (s in seq_len(S)) {
    cols[[sprintf("kappa_p[%d]", s)]]      <- rnorm(nd, 0.2)
    cols[[sprintf("kappa_m_base[%d]", s)]] <- rnorm(nd, 0.6)
    cols[[sprintf("mu_const[%d]", s)]]     <- rnorm(nd, 0.0)
    cols[[sprintf("sigma_p[%d]", s)]]      <- rnorm(nd, 0.3)
  }
  cols[["m1"]] <- rnorm(nd, 0.5); cols[["sigma_phi_meas"]] <- rnorm(nd, 0.2)
  l2 <- bayesianOU:::.extract_level2_summary(mock_fit_df(cols), S, Tn, NULL)
  expect_null(l2$nu_p)
  expect_null(l2$sigma_p_t)
  expect_null(l2$a3_p)
})

test_that("fit_ou_nested_mi guards n_levels and requires V_value for Level 3 (D-IMPL-10)", {
  ph <- replicate(3, matrix(stats::rnorm(10 * 2), 10, 2), simplify = FALSE)
  X  <- matrix(stats::rnorm(10 * 2), 10, 2)
  args0 <- list(phi_draws = ph, X = X, TMG = stats::rnorm(10),
                COM = matrix(1, 10, 2), CAPITAL_TOTAL = matrix(1, 10, 2),
                Gprime = stats::rnorm(10), M = 2L)
  expect_error(do.call(fit_ou_nested_mi, c(args0, list(n_levels = 4L))),
               "1, 2 or 3")
  expect_error(do.call(fit_ou_nested_mi, c(args0, list(n_levels = 3L))),
               "V_value")
})

test_that("ou_nested_mi print/summary/plot work", {
  Tn <- 10L; S <- 2L; mat <- function(v) matrix(v, Tn, S)
  x <- structure(list(
    rubin = data.frame(parameter = c("kappa_p[1]", "m1"), estimate = c(0.2, 0.5),
                       total_sd = c(0.05, 0.1), q2.5 = c(0.1, 0.3),
                       q97.5 = c(0.3, 0.7), df = c(20, 18), fmi = c(0.3, 0.2)),
    phi_latent_pooled = list(mean = mat(0.1), sd = mat(0.2),
                             q2.5 = mat(-0.1), q97.5 = mat(0.3)),
    separation_pooled = list(prob_sep_by_sector = c(0.9, 0.8), prob_sep_joint = 0.78),
    per_imputation = list(list(rhat_max = 1.02), list(rhat_max = 1.01)),
    config = list(M = 2L, n_levels = 2L, n_available = 4L)
  ), class = c("ou_nested_mi", "list"))
  expect_output(print(x), "ou_nested_mi")
  s <- summary(x); expect_s3_class(s, "summary.ou_nested_mi")
  expect_output(print(s), "Rubin-pooled")
  pf <- tempfile(fileext = ".pdf"); grDevices::pdf(pf)
  on.exit({ grDevices::dev.off(); unlink(pf) }, add = TRUE)
  expect_silent(plot(x, type = "phi", sector = 1))
  expect_silent(plot(x, type = "fmi"))
})
