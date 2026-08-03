#' Plot Method for the Package 'quarks'
#'
#' This function regulates how objects created by the package quarks are
#' plotted.
#'
#' @param x an input object of class \code{quarks}.
#' @param ... additional arguments of the standard plot method.
#'
#' @export
#'
#' @return
#' None

plot.quarks <- function(x, ...) {
  dots <- list(...)

  if(attr(x, "function") == "rollcast") {
    dots[["x"]] <- seq_along(x[["xout"]])
    VaR.plot <- ifelse(-x[["xout"]] >= x[["VaR"]], x[["VaR"]], NA)
    ES.plot <- ifelse(-x[["xout"]] >= x[["ES"]], x[["ES"]], NA)

    dots[["y"]] <- cbind(-x[["xout"]], x[["VaR"]], x[["ES"]], VaR.plot, ES.plot)
    if (is.null(dots[["col"]])) {
      dots[["col"]] <- c("darkgrey", "red", "green", "blue", "black")
    }
    if(is.null(dots[["panel.first"]])) {
      dots[["panel.first"]] = quote(graphics::grid())
    }
    if (is.null(dots[["type"]])) {
      dots[["type"]] <- "hllpp"
    }
    if (is.null(dots[["lty"]])) {
    dots[["lty"]] <- c(1, 5, 5)
    }
    if (is.null(dots[["pch"]])) {
    dots[["pch"]] <- 1
    }
    if (is.null(dots[["xlab"]])) {
    dots[["xlab"]] <- "out-of-sample observations"
    }
    if (is.null(dots[["ylab"]])) {
    dots[["ylab"]] <- "Negative returns, VaR and ES"
    }
    if (is.null(dots[["main"]])) {
    dots[["main"]] <- "VaR (red) and ES (green) - one-step forecasts"
    }

    do.call(what = graphics::matplot, args = dots)
  }
}

