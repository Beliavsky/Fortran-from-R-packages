#' Application for downloading data from Yahoo Finance
#'
#' @export
#' @import shiny
#' @importFrom shinyjs 'useShinyjs'
#' @import dygraphs
#' @import ggplot2
#'
#' @return
#' None

runFTSdata <- function() {
  options(show.error.messages = FALSE)
  app <-
    list(ui = fluidPage(useShinyjs(),
                              sidebarLayout(
                                sidebarPanel(
                                  selectizeInput("ticker",
                                                 label = "Enter a ticker symbol or name",
                                                 choices = NULL,
                                                 options = list(create = TRUE)),

                                  dateRangeInput("Date",
                                                 "Period",
                                                 start = Sys.Date()-750,
                                                 end = Sys.Date()),
                                  br(),
                                  br(),
                                  h4("Upper Graph", style = "text-align: left;"),
                                  fluidRow(
                                    column(6,
                                  selectInput("dropdown3",
                                              label = "Choose frequency",
                                              choices = c("daily", "weekly", "monthly",
                                                          "yearly"),
                                              selected = "daily"),
                                  textOutput("selected_option"),
                                  # p(actionButton("loadData", "Load", icon = icon("r-project")),
                                  #   br(),
                                  #   "Selected data is loaded into global R-environment.",
                                  #   style = "font-family: 'consolas'; font-si16pt"),
                                  #br(),
                                  #br(),
                                    ),
                                  column(6,
                                  selectInput("dropdown",
                                              label = "Choose FTS",
                                              choices = c("no FTS selected", "price_open", "price_close",
                                                          "price_adjusted",
                                                          "price_high", "price_low",
                                                          "volume", "ret_adjusted_prices",
                                                          "ret_closing_prices",
                                                          "cumret_adjusted_prices"),
                                              selected = "price_close"),
                                  textOutput("selected_option"),
                                  br(),
                                  br(),
                                    )
                                  ),
                                  h4("Lower Graph", style = "text-align: left;"),
                                  fluidRow(
                                    column(6,
                                  selectInput("dropdown4",
                                              label = "Choose frequency",
                                              choices = c("daily", "weekly", "monthly",
                                                          "yearly"),
                                              selected = "daily"),
                                  textOutput("selected_option"),
                                    ),
                                  column(6,
                                  selectInput("dropdown2",
                                              label = "Choose FTS",
                                              choices = c("no FTS selected", "price_open", "price_close",
                                                          "price_adjusted",
                                                          "price_high", "price_low",
                                                          "volume", "ret_adjusted_prices",
                                                          "ret_closing_prices",
                                                          "cumret_adjusted_prices"),
                                              selected = "ret_closing_prices"),
                                  textOutput("selected_option"),
                                  br(),
                                  br(),
                                    )
                                  ),
                                  p(downloadButton("downloadData", "Save", icon = icon("floppy-disk")),
                                    br(),
                                    "Save selected data on local hard drive",
                                    style = "font-family: 'consolas'; font-si16pt"),
                                  br(),
                                  br(),
                                  p(actionButton("closeApp", "Close", icon = icon("power-off")),
                                    br(),
                                    "Close Application",
                                    style = "font-family: 'consolas'; font-si16pt"),
                                ),
                                mainPanel(br(),
                                          br(),
                                          conditionalPanel(
                                            condition = "input.dropdown !== 'no FTS selected'",
                                          dygraphOutput("plot1"),
                                          ),
                                          br(),
                                          br(),
                                          conditionalPanel(
                                            condition = "input.dropdown2 !== 'no FTS selected'",
                                          dygraphOutput("plot2"),
                                          ),
                                          )
                              )
  ),
  server = function(input, output, session){
    updateSelectizeInput(session, "ticker",
                         choices = c(paste0(ticker_stocks$Ticker[-1], " (", ticker_stocks$Name[-1], ")"),
                                     paste0(ticker_indices$Ticker, " (", ticker_indices$Name, ")")),
                         selected = '^GDAXI (DAX)',
                         server = TRUE)

    dataInput <- reactive({
      req(input$ticker)
      req(input$Date[1])
      req(input$Date[2])
      req(input$dropdown3)
      ticker = strsplit(input$ticker, " ")[[1]][1]
      suppressWarnings(
        suppressMessages(
          yfR::yf_get(tickers = ticker, first_date = input$Date[1],
                      last_date = input$Date[2], freq_data = input$dropdown3,
                      cache_folder = file.path(tempdir(), 'BGS_Cache'))
        )
      )
    })

    dataInput2 <- reactive({
      req(input$ticker)
      req(input$Date[1])
      req(input$Date[2])
      req(input$dropdown4)
      ticker = strsplit(input$ticker, " ")[[1]][1]
      suppressWarnings(
        suppressMessages(
          yfR::yf_get(tickers = ticker, first_date = input$Date[1],
                      last_date = input$Date[2], freq_data = input$dropdown4,
                      cache_folder = file.path(tempdir(), 'BGS_Cache'))
        )
      )
    })

  output$plot1 <- renderDygraph({
    req(input$ticker)
    req(input$Date[1])
    req(input$Date[2])
    req(input$dropdown)
    req(input$dropdown3)
    req(dataInput())
    dat = dataInput()
    plot_data = xts::xts(dat[paste(input$dropdown)], dat$ref_date)
    dygraphs::dygraph(plot_data, main = paste(input$ticker, input$dropdown, input$dropdown3),
                      #ylab = "X<sub>t</sub>",
                      group = "FZS1") %>%
      dygraphs::dyRangeSelector(height = 20) %>%
      dygraphs::dyUnzoom()
  })


    output$plot2 <- renderDygraph({
      req(input$ticker)
      req(input$Date[1])
      req(input$Date[2])
      req(input$dropdown)
      req(input$dropdown4)
      req(dataInput2())

      dat2 = dataInput2()#[-1, ]
      plot_data2 = xts::xts(dat2[paste(input$dropdown2)], dat2$ref_date)
      dygraphs::dygraph(plot_data2, main = paste(input$ticker, input$dropdown2, input$dropdown4),
                        #ylab = "Y<sub>t</sub>",
                        #xlab = "Daily observations (t)",
                        group = "FZS2") %>%
        dygraphs::dyRangeSelector(height = 20) %>%
        dygraphs::dyUnzoom()
    })

    # observeEvent(input$loadData, {
    #   data_raw = dataInput()
    #   ticker = strsplit(input$ticker, " ")[[1]][1]
    #   if (length(strsplit(ticker, "^", fixed = TRUE)[[1]]) != 1) {
    #     ticker = strsplit(ticker, "^", fixed = TRUE)[[1]][2]
    #   }
    #   assign(ticker, data_raw)#, envir = globalenv())
    #   cat("", fill = TRUE)
    #   cat("Data for", input$ticker, "loaded as object", ticker, "in global environment.", fill = TRUE)
    #   cat("", fill = TRUE)
    # })

    output$downloadData <- downloadHandler(
      filename = function() {
        ticker = strsplit(input$ticker, " ")[[1]][1]
        if (length(strsplit(ticker, "^", fixed = TRUE)[[1]]) != 1) {
          ticker = strsplit(ticker, "^", fixed = TRUE)[[1]][2]
        }
        paste0(ticker, ".csv")
      },
      content = function(file) {
        utils::write.csv(dataInput(), file, row.names = FALSE)
      }
    )
    observeEvent(input$closeApp, {
      stopApp()
    })
    }
)
shiny::runApp(app, display.mode = 'normal')
}
