# ============================================================================
# S3 methods for the 2-level nested fit (ou_nested_2level) and the
# multiple-imputation driver (ou_nested_mi). Base-R plots only (no ggplot2
# dependency), matching the rest of the package's plotting helpers.
# ============================================================================

#' @keywords internal
#' @noRd
.fmt_band <- function(m) {
  if (is.null(m)) return("NA")
  sprintf("%.3f [%.3f, %.3f]", m[, "median"], m[, "q2.5"], m[, "q97.5"])
}

#' @keywords internal
#' @noRd
.range_str <- function(v) {
  v <- v[is.finite(v)]
  if (!length(v)) return("NA")
  sprintf("[%.3f, %.3f]", min(v), max(v))
}


#' Print a 2-level nested OU fit
#'
#' Compact one-screen overview: data dimensions, MCMC health, the time-scale
#' separation evidence \eqn{\kappa^m > \kappa^p}, and the median ranges of the
#' market and production reversion speeds.
#'
#' @param x An \code{ou_nested_2level} object from \code{\link{fit_ou_nested}}.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @export
print.ou_nested_2level <- function(x, ...) {
  fo  <- x$factor_ou; l2 <- fo$level2; dg <- x$diagnostics
  inf <- fo$factor_ou_info
  Tn  <- if (!is.null(x$phi_latent)) nrow(x$phi_latent$median) else NA_integer_
  S   <- if (!is.null(x$phi_latent)) ncol(x$phi_latent$median) else NA_integer_

  cat("<ou_nested_2level>  2-level nested gravitation OU\n")
  cat(sprintf("  market phi  ->  latent production Phi  ->  G'-driven mean mu\n"))
  cat(sprintf("  T = %s, S = %s | T_train = %s, fit_window = %s\n",
              Tn, S, inf$T_train %||% NA, inf$fit_window %||% NA))
  cat(sprintf("  theta_separation = %s | k_uncertainty = %s | kappa_cap implied\n",
              inf$theta_separation %||% NA, inf$k_uncertainty %||% NA))
  cat("  -- MCMC --\n")
  cat(sprintf("  max R-hat = %.4f | share(R-hat>1.01) = %.3f | divergences = %s\n",
              dg$rhat_max %||% NA_real_, dg$rhat_share %||% NA_real_,
              dg$divergences %||% NA))
  cat("  -- Time-scale separation kappa^m > kappa^p --\n")
  if (!is.null(x$separation)) {
    cat(sprintf("  P(separation, joint) = %.3f | by sector: %s\n",
                x$separation$prob_sep_joint,
                paste(sprintf("%.2f", x$separation$prob_sep_by_sector),
                      collapse = " ")))
  }
  cat("  -- Reversion speeds (median range over sectors) --\n")
  if (!is.null(l2)) {
    cat(sprintf("  kappa_m_base %s | kappa_p %s | m1 = %.3f\n",
                .range_str(l2$kappa_m_base[, "median"]),
                .range_str(l2$kappa_p[, "median"]),
                l2$m1 %||% NA_real_))
  }
  invisible(x)
}


#' Summarize a 2-level nested OU fit
#'
#' Assembles the Level-2 parameter table (production reversion speed, market
#' speed at \eqn{zTMG=0}, mean intercept, innovation scale, measurement scale),
#' the separation evidence, the MCMC diagnostics and the out-of-sample metrics.
#'
#' @param object An \code{ou_nested_2level} object from \code{\link{fit_ou_nested}}.
#' @param ... Ignored.
#' @return An object of class \code{summary.ou_nested_2level}.
#' @export
summary.ou_nested_2level <- function(object, ...) {
  fo <- object$factor_ou; l2 <- fo$level2; dg <- object$diagnostics
  band <- function(m, name) {
    if (is.null(m)) return(NULL)
    df <- as.data.frame(m); df$parameter <- sprintf("%s[%d]", name, seq_len(nrow(df)))
    df[, c("parameter", "q2.5", "median", "q97.5")]
  }
  # The Level-2 richness rows (a3_p cubic, nu_p Student-t df) are present only
  # when their switch is on; band(NULL) returns NULL and is dropped, so the
  # canonical table is unchanged (D-IMPL-10).
  level2_table <- do.call(rbind, Filter(Negate(is.null), list(
    band(l2$kappa_p,      "kappa_p"),
    band(l2$kappa_m_base, "kappa_m_base"),
    band(l2$mu_const,     "mu_const"),
    band(l2$sigma_p,      "sigma_p"),
    band(l2$a3_p,         "a3_p"),
    band(l2$nu_p,         "nu_p")
  )))
  out <- list(
    info         = fo$factor_ou_info,
    level2_table = level2_table,
    m1           = l2$m1,
    sigma_phi_meas = l2$sigma_phi_meas,
    sigma_p_t    = l2$sigma_p_t,          # time-varying L2 SV scale (NULL if l2_sv off)
    separation   = object$separation,
    diagnostics  = list(rhat_max = dg$rhat_max, rhat_share = dg$rhat_share,
                        divergences = dg$divergences,
                        loo_pareto_k = dg$loo_pareto_k),
    oos          = dg$oos
  )
  class(out) <- "summary.ou_nested_2level"
  out
}

#' @rdname summary.ou_nested_2level
#' @param x A \code{summary.ou_nested_2level} object.
#' @export
print.summary.ou_nested_2level <- function(x, ...) {
  cat("Summary: 2-level nested gravitation OU\n")
  cat(sprintf("  theta_separation = %s | k_uncertainty = %s\n\n",
              x$info$theta_separation %||% NA, x$info$k_uncertainty %||% NA))
  cat("-- Level-2 parameters (95% credible intervals) --\n")
  if (!is.null(x$level2_table)) {
    print(format(x$level2_table, digits = 3), row.names = FALSE)
  }
  cat(sprintf("\n  m1 (G'-driven mean slope) median = %.4f\n", x$m1 %||% NA_real_))
  cat(sprintf("  sigma_phi_meas (anchor measurement SD) median = %.4f\n",
              x$sigma_phi_meas %||% NA_real_))
  if (!is.null(x$sigma_p_t)) {
    cat(sprintf("  sigma_p(t) (Level-2 SV scale) median range = %s\n",
                .range_str(as.vector(x$sigma_p_t$median))))
  }
  if (!is.null(x$separation)) {
    cat(sprintf("\n-- Separation kappa^m > kappa^p --\n  P(joint) = %.3f\n",
                x$separation$prob_sep_joint))
  }
  cat(sprintf("\n-- MCMC --\n  max R-hat = %.4f | divergences = %s\n",
              x$diagnostics$rhat_max %||% NA_real_, x$diagnostics$divergences %||% NA))
  if (!is.null(x$oos)) {
    cat("\n-- Out-of-sample (standardized phi units) --\n")
    for (nm in names(x$oos)) {
      o <- x$oos[[nm]]
      cat(sprintf("  h=%-2d RMSE = %s  MAE = %s  (n = %s)\n", o$h,
                  if (is.na(o$RMSE)) "NA" else sprintf("%.4f", o$RMSE),
                  if (is.na(o$MAE)) "NA" else sprintf("%.4f", o$MAE), o$n_obs))
    }
  }
  invisible(x)
}


#' Plot a 2-level nested OU fit
#'
#' @param x An \code{ou_nested_2level} object from \code{\link{fit_ou_nested}}.
#' @param type Character. \code{"phi"} (latent production-price path with 95\%
#'   band), \code{"mu"} (the \eqn{G'}-driven Level-2 mean trajectory with band),
#'   \code{"sv_p"} (the time-varying Level-2 stochastic-volatility scale
#'   \eqn{\sigma_p(t)} with band; only available when the Level-2 SV switch is on)
#'   or \code{"separation"} (per-sector posterior probability of \eqn{\kappa^m >
#'   \kappa^p}). Default \code{"phi"}.
#' @param sector Integer. Sector index for \code{"phi"} / \code{"mu"} /
#'   \code{"sv_p"}. Default 1.
#' @param ... Passed to the underlying base graphics call.
#' @return \code{NULL} invisibly; draws as a side effect.
#' @export
plot.ou_nested_2level <- function(x, type = c("phi", "mu", "sv_p", "separation"),
                                  sector = 1, ...) {
  type <- match.arg(type)
  band_plot <- function(med, lo, hi, s, ylab, main) {
    tt <- seq_len(nrow(med))
    graphics::plot(tt, med[, s], type = "n",
                   ylim = range(c(lo[, s], hi[, s]), na.rm = TRUE),
                   xlab = "t", ylab = ylab, main = main, ...)
    graphics::polygon(c(tt, rev(tt)), c(lo[, s], rev(hi[, s])),
                      col = grDevices::adjustcolor("steelblue", 0.25), border = NA)
    graphics::lines(tt, med[, s], lwd = 2, col = "steelblue4")
    graphics::grid()
  }
  if (type == "phi") {
    p <- x$phi_latent
    if (is.null(p)) stop("No latent Phi trajectory in this object.", call. = FALSE)
    band_plot(p$median, p$q2.5, p$q97.5, sector,
              expression(Phi[t]),
              sprintf("Latent production price (sector %d)", sector))
  } else if (type == "mu") {
    m <- x$mu_path
    if (is.null(m)) stop("No mu trajectory in this object.", call. = FALSE)
    band_plot(m$median, m$q2.5, m$q97.5, sector,
              expression(mu[t]),
              sprintf("Level-2 mean mu_t = m0 + m1*G' (sector %d)", sector))
  } else if (type == "sv_p") {
    sp <- x$factor_ou$level2$sigma_p_t
    if (is.null(sp)) {
      stop("No Level-2 SV scale path in this object (fit with the Level-2 ",
           "stochastic-volatility switch on, e.g. ou_level_spec(\"both_full\")).",
           call. = FALSE)
    }
    band_plot(sp$median, sp$q2.5, sp$q97.5, sector,
              expression(sigma[p](t)),
              sprintf("Level-2 SV scale sigma_p(t) (sector %d)", sector))
  } else {
    sep <- x$separation
    if (is.null(sep)) stop("No separation evidence in this object.", call. = FALSE)
    bp <- graphics::barplot(sep$prob_sep_by_sector, ylim = c(0, 1),
                            names.arg = seq_along(sep$prob_sep_by_sector),
                            xlab = "sector",
                            ylab = expression(P(kappa^m > kappa^p)),
                            main = "Time-scale separation by sector", ...)
    graphics::abline(h = sep$prob_sep_joint, lty = 2, col = "firebrick")
    graphics::text(max(bp), sep$prob_sep_joint,
                   sprintf("joint = %.2f", sep$prob_sep_joint),
                   pos = 3, col = "firebrick", cex = 0.8)
  }
  invisible(NULL)
}


# --------------------------------------------------------------------------
# ou_nested_mi (multiple-imputation driver)
# --------------------------------------------------------------------------

#' Print a multiple-imputation nested OU fit
#'
#' @param x An \code{ou_nested_mi} object from \code{\link{fit_ou_nested_mi}}.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @export
print.ou_nested_mi <- function(x, ...) {
  cfg <- x$config
  cat("<ou_nested_mi>  nested OU under multiple imputation (Rubin's rule)\n")
  cat(sprintf("  M = %d imputations (of %d available) | n_levels = %d\n",
              cfg$M, cfg$n_available %||% cfg$M, cfg$n_levels))
  rh <- vapply(x$per_imputation,
               function(z) if (is.list(z)) z$rhat_max %||% NA_real_ else NA_real_,
               numeric(1))
  if (any(is.finite(rh))) {
    cat(sprintf("  per-imputation max R-hat: median = %.4f, worst = %.4f\n",
                stats::median(rh, na.rm = TRUE), max(rh, na.rm = TRUE)))
  }
  if (!is.null(x$separation_pooled)) {
    cat(sprintf("  P(kappa^m > kappa^p, mixture posterior) = %.3f\n",
                x$separation_pooled$prob_sep_joint))
  }
  cat(sprintf("  Rubin-pooled parameters: %d rows (see summary())\n",
              if (is.null(x$rubin)) 0L else nrow(x$rubin)))
  invisible(x)
}

#' Summarize a multiple-imputation nested OU fit
#'
#' @param object An \code{ou_nested_mi} object from \code{\link{fit_ou_nested_mi}}.
#' @param ... Ignored.
#' @return An object of class \code{summary.ou_nested_mi}.
#' @export
summary.ou_nested_mi <- function(object, ...) {
  out <- list(
    config = object$config,
    rubin  = object$rubin,
    separation_pooled = object$separation_pooled,
    fmi_summary = if (!is.null(object$rubin)) {
      stats::quantile(object$rubin$fmi, probs = c(0, 0.5, 1), na.rm = TRUE)
    } else NULL
  )
  class(out) <- "summary.ou_nested_mi"
  out
}

#' @rdname summary.ou_nested_mi
#' @param x A \code{summary.ou_nested_mi} object.
#' @export
print.summary.ou_nested_mi <- function(x, ...) {
  cat("Summary: nested OU under multiple imputation\n")
  cat(sprintf("  M = %d | n_levels = %d\n\n", x$config$M, x$config$n_levels))
  cat("-- Rubin-pooled estimates --\n")
  if (!is.null(x$rubin)) print(format(x$rubin, digits = 3), row.names = FALSE)
  if (!is.null(x$fmi_summary)) {
    cat(sprintf("\n  Fraction of missing information (min/median/max): %.3f / %.3f / %.3f\n",
                x$fmi_summary[1], x$fmi_summary[2], x$fmi_summary[3]))
  }
  if (!is.null(x$separation_pooled)) {
    cat(sprintf("\n  P(kappa^m > kappa^p, mixture) = %.3f\n",
                x$separation_pooled$prob_sep_joint))
  }
  invisible(x)
}

#' Plot a multiple-imputation nested OU fit
#'
#' @param x An \code{ou_nested_mi} object from \code{\link{fit_ou_nested_mi}}.
#' @param type Character. \code{"phi"} (Rubin-pooled latent production price with
#'   95\% band) or \code{"fmi"} (fraction of missing information per pooled
#'   parameter). Default \code{"phi"}.
#' @param sector Integer. Sector index for \code{"phi"}. Default 1.
#' @param ... Passed to the underlying base graphics call.
#' @return \code{NULL} invisibly; draws as a side effect.
#' @export
plot.ou_nested_mi <- function(x, type = c("phi", "fmi"), sector = 1, ...) {
  type <- match.arg(type)
  if (type == "phi") {
    p <- x$phi_latent_pooled
    if (is.null(p)) stop("No pooled latent Phi in this object.", call. = FALSE)
    tt <- seq_len(nrow(p$mean))
    graphics::plot(tt, p$mean[, sector], type = "n",
                   ylim = range(c(p$q2.5[, sector], p$q97.5[, sector]), na.rm = TRUE),
                   xlab = "t", ylab = expression(Phi[t]),
                   main = sprintf("Rubin-pooled latent Phi (sector %d)", sector), ...)
    graphics::polygon(c(tt, rev(tt)),
                      c(p$q2.5[, sector], rev(p$q97.5[, sector])),
                      col = grDevices::adjustcolor("darkorange", 0.25), border = NA)
    graphics::lines(tt, p$mean[, sector], lwd = 2, col = "darkorange3")
    graphics::grid()
  } else {
    if (is.null(x$rubin)) stop("No Rubin table in this object.", call. = FALSE)
    graphics::barplot(x$rubin$fmi, names.arg = x$rubin$parameter, las = 2,
                      ylab = "fraction of missing information",
                      main = "FMI by pooled parameter", ...)
  }
  invisible(NULL)
}
