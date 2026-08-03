server <- function(input, output, session) {

  shiny::observeEvent(input$res, {

    plot_res <- input$res
    output$purchasing_power_plot <- shiny::renderPlot({

      asNamespace("R4GoodPersonalFinances")$plot_purchasing_power(
        x = input$x ,
        real_interest_rate = input$real_interest_rate / 100,
        years = input$years
      )
    }, res = plot_res)
  })
}
