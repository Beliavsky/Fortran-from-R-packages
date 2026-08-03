#' Run a package app
#' 
#' @param which A character. The name of the app to run.
#' Currently available:
#' 
#' * `risk-adjusted-returns` - Plotting risk-adjusted returns for various allocations to the risky asset allows you to find the optimal allocation.
#' * `purchasing-power` - Plotting the effect of real interest rates 
#' (positive or negative) on the purchasing power of savings over time.
#' * `retirement-ruin` - Plotting the probability of retirement ruin.
#' 
#' @param res A numeric. The initial resolution of the plots.
#' 
#' @param shinylive A logical. Whether to use `shinylive` for the app.
#' 
#' @returns A [shiny::shinyApp()] object if `shinylive` is `TRUE`.
#' Runs the app if `shinylive` is `FALSE` with [shiny::runApp()].
#'
#' @examplesIf interactive()
#' run_app("risk-adjusted-returns")
#' run_app("purchasing-power")
#' run_app("retirement-ruin")
#' @export

run_app <- function(
  which = c(
    "risk-adjusted-returns",
    "purchasing-power",
    "retirement-ruin"
  ),
  res = 120,
  shinylive = FALSE
) {

  if (!requireNamespace("shiny", quietly = TRUE)) {
    cli::cli_abort("The 'shiny' package is required but not installed.")
  }

  which <- match.arg(which)

  old_options <- options(R4GPF.plot_res = res)
  on.exit(options(old_options))

  if (shinylive) {

    app <- apps[[which]]
    temp_dir <- tempdir()

    temp_path_ui <- file.path(temp_dir, "ui.R")
    writeLines(app$ui, temp_path_ui)

    temp_path_server <- file.path(temp_dir, "server.R")
    writeLines(app$server, temp_path_server)

    ui     <- source(temp_path_ui,     local = TRUE)$value
    server <- source(temp_path_server, local = TRUE)$value

    shiny::shinyApp(ui = ui, server = server)

  } else {

    shiny::runApp(
      system.file("apps", which, package = "R4GoodPersonalFinances"),
      display.mode = "normal"
    )
  }
}
