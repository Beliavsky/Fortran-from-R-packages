ui <- 
  bslib::page_sidebar(
    title = "Probability of Retirement Ruin",
    sidebar = bslib::sidebar(
      width = 250,

      shiny::numericInput(
        inputId = "age",
        label = "Current age",
        step = 1,
        value = 65
      ),
      shiny::numericInput(
        inputId = "portfolio_value",
        label = "Initial portfolio value",
        step = 1,
        value = 1000000
      ),
      shiny::numericInput(
        inputId = "monthly_spendings",
        label = "Monthly spendings (in thousands)",
        step = 0.1,
        value = NA      
      ),
      shiny::numericInput(
        inputId = "portfolio_return_mean",
        label = "Mean of portfolio returns (%)",
        step = 0.1,
        value = 3.4
      ),
      shiny::numericInput(
        inputId = "portfolio_return_sd",
        label = "Standard deviation of portfolio returns (%)",
        step = 0.1,
        value = 15
      ),
      shiny::numericInput(
        inputId = "gompertz_mode",
        label = "Gompertz mode",
        step = 1,
        value = 84
      ),
      shiny::numericInput(
        inputId = "gompertz_dispersion",
        label = "Gompertz dispersion",
        step = 0.1,
        value = 10
      ),
      asNamespace("R4GoodPersonalFinances")$sidebar_footer()
    ),
    
    bslib::card(
      max_height = "650px",
      full_screen = TRUE,
      shiny::plotOutput("retirement_ruin_plot"),
      bslib::card_footer(
        bslib::popover(
          bsicons::bs_icon("gear"),
          shiny::numericInput(
            inputId = "res", 
            value = getOption("R4GPF.plot_res", default = 120), 
            label = "Plot resolution"
          ),
          title = "Settings"
        )
      )
    )
  )
