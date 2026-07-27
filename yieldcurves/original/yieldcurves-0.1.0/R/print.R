#' Print Method for Yield Curve Objects
#'
#' @param x A `yc_curve` object.
#' @param ... Additional arguments (currently unused).
#'
#' @return The input object, invisibly.
#'
#' @export
print.yc_curve <- function(x, ...) {
  method_names <- c(
    observed      = "Observed",
    nelson_siegel = "Nelson-Siegel",
    svensson      = "Svensson",
    cubic_spline  = "Cubic Spline"
  )
  method_label <- method_names[x$method]

  cli_h1("Yield Curve ({method_label})")

  bullets <- c(
    "*" = "Type: {.val {x$type}}",
    "*" = "Maturities: {.val {x$n_obs}} ({min(x$maturities)}Y to {max(x$maturities)}Y)",
    "*" = "Rate range: {.val {round(min(x$rates) * 100, 2)}}% to {.val {round(max(x$rates) * 100, 2)}}%"
  )

  if (!is.null(x$date)) {
    bullets <- c(bullets, "*" = "Date: {.val {x$date}}")
  }

  if (!is.null(x$residuals)) {
    rmse <- sqrt(mean(x$residuals^2))
    bullets <- c(bullets,
      "*" = "RMSE: {.val {format(rmse * 10000, digits = 2)}} bps"
    )
  }

  cli_bullets(bullets)

  if (x$method == "nelson_siegel") {
    cli_bullets(c(
      "*" = "Parameters: beta0={.val {round(x$params$beta0, 5)}}, beta1={.val {round(x$params$beta1, 5)}}, beta2={.val {round(x$params$beta2, 5)}}, tau={.val {round(x$params$tau, 3)}}"
    ))
  } else if (x$method == "svensson") {
    cli_bullets(c(
      "*" = "Parameters: beta0={.val {round(x$params$beta0, 5)}}, beta1={.val {round(x$params$beta1, 5)}}, beta2={.val {round(x$params$beta2, 5)}}, beta3={.val {round(x$params$beta3, 5)}}, tau1={.val {round(x$params$tau1, 3)}}, tau2={.val {round(x$params$tau2, 3)}}"
    ))
  }

  invisible(x)
}

#' Summary Method for Yield Curve Objects
#'
#' @param object A `yc_curve` object.
#' @param ... Additional arguments (currently unused).
#'
#' @return The input object, invisibly.
#'
#' @export
summary.yc_curve <- function(object, ...) {
  print(object, ...)

  if (!is.null(object$residuals)) {
    cli_h1("Fit Quality")
    resid_bps <- object$residuals * 10000
    cli_bullets(c(
      "*" = "Max absolute error: {.val {round(max(abs(resid_bps)), 2)}} bps",
      "*" = "Mean absolute error: {.val {round(mean(abs(resid_bps)), 2)}} bps",
      "*" = "Residual range: [{.val {round(min(resid_bps), 2)}}, {.val {round(max(resid_bps), 2)}}] bps"
    ))
  }

  invisible(object)
}

#' Plot Method for Yield Curve Objects
#'
#' @param x A `yc_curve` object.
#' @param ... Additional arguments passed to [plot()].
#'
#' @return The input object, invisibly.
#'
#' @export
plot.yc_curve <- function(x, ...) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))

  r_pct <- x$rates * 100
  ylim <- range(r_pct)

  if (!is.null(x$fitted)) {
    fitted_pct <- x$fitted * 100
    ylim <- range(c(r_pct, fitted_pct))
  }

  # Expand ylim slightly
  ylim <- ylim + c(-0.05, 0.05) * diff(ylim)

  method_names <- c(
    observed      = "Observed",
    nelson_siegel = "Nelson-Siegel",
    svensson      = "Svensson",
    cubic_spline  = "Cubic Spline"
  )

  title <- paste0("Yield Curve (", method_names[x$method], ")")
  if (!is.null(x$date)) title <- paste0(title, " - ", x$date)

  plot(x$maturities, r_pct, type = "p", pch = 19,
       xlab = "Maturity (years)", ylab = "Yield (%)",
       main = title, ylim = ylim, ...)
  grid()

  if (!is.null(x$fitted)) {
    # Plot smooth fitted curve
    m_smooth <- seq(min(x$maturities), max(x$maturities), length.out = 200)
    r_smooth <- yc_predict(x, m_smooth)$rate * 100
    lines(m_smooth, r_smooth, col = "#d95f02", lwd = 2)
    legend("topright", legend = c("Observed", "Fitted"),
           col = c("black", "#d95f02"), pch = c(19, NA),
           lty = c(NA, 1), lwd = c(NA, 2), bty = "n")
  } else {
    lines(x$maturities, r_pct, col = "black", lwd = 1)
  }

  invisible(x)
}
