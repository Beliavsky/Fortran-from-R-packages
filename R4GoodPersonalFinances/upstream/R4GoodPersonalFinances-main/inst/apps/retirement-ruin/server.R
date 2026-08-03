server <- function(input, output, session) {

  shiny::observe({

    plot_res <- input$res
    output$retirement_ruin_plot <- shiny::renderPlot({

      shiny::req(input$portfolio_return_mean)
      shiny::req(input$portfolio_return_sd)
      shiny::req(input$age)
      shiny::req(input$gompertz_mode)
      shiny::req(input$gompertz_dispersion)
      shiny::req(input$portfolio_value)
      
      if (shiny::isTruthy(input$monthly_spendings)) {
        monthly_spendings <- input$monthly_spendings * 1000
      } else {
        monthly_spendings <- NULL
      }

      asNamespace("R4GoodPersonalFinances")$plot_retirement_ruin(
        portfolio_return_mean = input$portfolio_return_mean / 100,
        portfolio_return_sd   = input$portfolio_return_sd / 100,
        age                   = input$age,
        gompertz_mode         = input$gompertz_mode,
        gompertz_dispersion   = input$gompertz_dispersion,
        portfolio_value       = input$portfolio_value,
        monthly_spendings     = monthly_spendings
      )
    }, res = plot_res)
  }) |> shiny::bindEvent(input$res)
}
