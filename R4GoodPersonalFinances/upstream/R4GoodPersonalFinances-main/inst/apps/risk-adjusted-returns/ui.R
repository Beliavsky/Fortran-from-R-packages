ui <- 
  bslib::page_sidebar(
    title = "Risk-adjusted Returns & Optimal Risky Asset Allocation",
    sidebar = bslib::sidebar(
      width = 250,
      shiny::sliderInput(
        inputId = "current_risky_asset_allocation",
        label = "Current risky asset allocation",
        step = 1,
        animate = shiny::animationOptions(
          interval = 700,
          loop = FALSE,
          playButton = NULL,
          pauseButton = NULL
        ),
        post = "%",
        min = 0, 
        max = 100, 
        value = 0
      ),
      shiny::numericInput(
        inputId = "risky_asset_return_mean",
        label = "Mean of yearly real returns of risky asset (%)",
        step = 0.1,
        value = 4
      ),
      shiny::numericInput(
        inputId = "risky_asset_return_sd",
        label = "Standard deviation of yearly real returns of risky asset (%)",
        step = 0.1,
        value = 15
      ),
      shiny::numericInput(
        inputId = "safe_asset_return",
        label = "Yearly real return of safe asset (%)",
        step = 0.1,
        value = 2
      ),
      shiny::numericInput(
        inputId = "risk_aversion",
        label = "Risk aversion (default = 2)",
        step = 0.1,
        value = 2
      ),
      asNamespace("R4GoodPersonalFinances")$sidebar_footer()
    ),
    
    bslib::card(
      max_height = "650px",
      full_screen = TRUE,
      shiny::plotOutput("rar_plot"),
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
