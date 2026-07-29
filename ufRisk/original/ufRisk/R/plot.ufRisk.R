#' Plot Method for the Package 'ufRisk'
#'
#'This function regulates how objects created by the package \code{ufRisk} are
#'plotted.
#'
#'@param x an input object of class \code{ufRisk}.
#'@param plot.option plot choice for an object of class \code{ufRisk}; viable
#'choices are:
#'\describe{
#'\item{1}{Plotting out-of-sample loss series}
#'\item{2}{Plotting out-of-sample losses, VaR.v & breaches}
#'\item{3}{Plotting out-of-sample losses, VaR.e, ES & breaches}
#'}
#'Please note if no value is passed to \code{plot.option} a selection menu is
#'prompted; is set to \code{NULL} by default.
#'@param ... additional arguments of the standard plot method.
#'
#'@export
#'
#'@importFrom graphics lines
#'
#'@return
#'None
#'
#'@author
#'\itemize{
#'\item Sebastian Letmathe (Scientific Employee) (Department of Economics,
#'Paderborn University) \cr
#'\item Dominik Schulz (Research Assistant) (Department of Economics, Paderborn
#'University), \cr
#'}
#'

plot.ufRisk <- function(x, plot.option = NULL, ...) {
  dots <- list(...)

  if (!is.null(plot.option) && (!is.numeric(plot.option) || plot.option < 1 ||
      plot.option > 3)) {
    stop("A single numeric value must be passed to 'plot.option'. ",
         "Valid choices are 1, 2 or 3.")
  }

  if (is.null(plot.option)) {
    cat("Plot choices for ufRisk object:", fill = TRUE)
    choices <- c(1, 2, 3)
    choice_names <- c("Out-of-sample loss series:",
                      "Out-of-sample losses, VaR.v & breaches:",
                      "Out-of-sample losses, VaR.e, ES & breaches:")
    choices_df <- data.frame(choices)
    colnames(choices_df) <- ""
    rownames(choices_df) <- choice_names
    print.data.frame(choices_df)
    plot_choice <- readline(prompt="Please enter the corresponding number: ")
    plot_choice <- as.numeric(plot_choice)
  }
  else {
    plot_choice <- plot.option
  }
  if(attr(x, "function") == "varcast") {

    dots[["x"]] <- seq_along(-x[["ret.out"]])

    if(is.null(dots[["panel.first"]])) {
      dots[["panel.first"]] <- quote(graphics::grid())
    }

    if(substr(plot_choice, 1, 1) == 1) {

      dots[["y"]] <- -x[["ret.out"]]

      if(is.null(dots[["col"]])) {
        dots[["col"]] <- "darkgray"
      }
      if(is.null(dots[["xlab"]])) {
        dots[["xlab"]] <- "Trading days"
      }
      if(is.null(dots[["ylab"]])) {
        dots[["ylab"]] <- "Negative returns"
      }
      if(is.null(dots[["main"]])) {
        dots[["main"]] <- "Out-of-sample losses"
      }
      if(is.null(dots[["type"]])) {
        dots[["type"]] <- "h"
      }
      do.call(what = graphics::matplot, args = dots)

    }

    if(substr(plot_choice, 1, 1) == 2) {

      VaR99.plot <- ifelse(-x[["ret.out"]] >= x[["VaR.v"]], x[["VaR.v"]], NA)

      dots[["y"]] <- cbind(-x[["ret.out"]], x[["VaR.v"]], VaR99.plot)

      if(is.null(dots[["lty"]])) {
        dots[["lty"]] <- c(1, 5)
      }
      if(is.null(dots[["col"]])) {
        dots[["col"]] <- c("darkgrey", "green", "black")
      }
      if(is.null(dots[["xlab"]])) {
        dots[["xlab"]] <- "Trading days"
      }
      if(is.null(dots[["ylab"]])) {
        dots[["ylab"]] <- paste0("Negative returns & VaR ",
                                (1 - x[["a.v"]]) * 100, "%")
      }
      if(is.null(dots[["main"]])) {
        dots[["main"]] <- paste0("Out-of-sample losses & VaR ",
                                (1 - x[["a.v"]]) * 100, "% (green)")
      }
      if(is.null(dots[["type"]])) {
        dots[["type"]] <- "hlp"
      }
      if(is.null(dots[["pch"]])) {
        dots[["pch"]] <- 1
      }

      do.call(what = graphics::matplot, args = dots)
      graphics::grid()

    }

    if(substr(plot_choice, 1, 1) == 3) {

      VaR975.plot <- ifelse(-x[["ret.out"]] >= x[["VaR.e"]], x[["VaR.e"]], NA)

      ES.plot <- ifelse(-x[["ret.out"]] >= x[["ES"]], x[["ES"]], NA)

      dots[["y"]] <- cbind(-x[["ret.out"]], x[["VaR.e"]], x[["ES"]],
                          VaR975.plot, ES.plot)

      if(is.null(dots[["lty"]])) {
        dots[["lty"]] <- c(1, 5, 5)
      }
      if(is.null(dots[["col"]])) {
        dots[["col"]] <- c("darkgrey", "red", "green", "blue", "black")
      }
      if(is.null(dots[["xlab"]])) {
        dots[["xlab"]] <- "Trading days"
      }
      if(is.null(dots[["ylab"]])) {
        dots[["ylab"]] <- paste0("Negative returns, VaR ",
                                (1 - x[["a.e"]]) * 100, "% & ES")
      }
      if(is.null(dots[["main"]])) {
        dots[["main"]] <- paste0("Out-of-sample losses, VaR ",
                                (1 - x[["a.e"]]) * 100, "% (red) & ES (green)")
      }
      if(is.null(dots[["type"]])) {
        dots[["type"]] <- "hllpp"
      }
      if(is.null(dots[["pch"]])) {
        dots[["pch"]] <- c(1, 1)
      }

      do.call(what = graphics::matplot, args = dots)
      graphics::grid()

    }
  }
}

