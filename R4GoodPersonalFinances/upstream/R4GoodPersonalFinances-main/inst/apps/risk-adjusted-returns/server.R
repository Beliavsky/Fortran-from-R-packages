server <- function(input, output, session) {

  shiny::observeEvent(input$res, {

    plot_res <- input$res
    
    output$rar_plot <- shiny::renderPlot({
      
      asNamespace("R4GoodPersonalFinances")$plot_risk_adjusted_returns(
        current_risky_asset_allocation = 
          input$current_risky_asset_allocation / 100,
        safe_asset_return       = input$safe_asset_return / 100,
        risky_asset_return_mean = input$risky_asset_return_mean / 100,
        risky_asset_return_sd   = input$risky_asset_return_sd / 100,
        risk_aversion           = input$risk_aversion
      )
    }, res = plot_res)
  })
}
