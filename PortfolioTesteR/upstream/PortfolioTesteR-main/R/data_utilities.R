#' Download S&P 500 Sector Mappings from Wikipedia
#'
#' @description
#' Scrapes current S&P 500 constituent list with sector classifications
#' from Wikipedia and returns as a data.table.
#'
#' @return Data.table with columns: Symbol, Security, Sector, SubIndustry, Industry
#' @export
#' @examplesIf (requireNamespace("rvest", quietly = TRUE))
#' \donttest{
#' sectors <- download_sp500_sectors()
#' head(sectors)
#' }
download_sp500_sectors <- function() {
  # Check if rvest is available
  if (!requireNamespace("rvest", quietly = TRUE)) {
    stop("Package 'rvest' needed. Install with: install.packages('rvest')")
  }

  url <- "https://en.wikipedia.org/wiki/List_of_S%26P_500_companies"
  page <- rvest::read_html(url)

  # Extract the main table - NO PIPES
  sp500_node <- rvest::html_node(page, "#constituents")
  sp500_table <- rvest::html_table(sp500_node)

  # Clean and structure
  sector_mapping <- data.table(
    Symbol = sp500_table$Symbol,
    Security = sp500_table$Security,
    Sector = sp500_table$`GICS Sector`,
    SubIndustry = sp500_table$`GICS Sub-Industry`
  )

  # Clean symbols
  sector_mapping[, Symbol := gsub("\\.", "-", Symbol)]
  sector_mapping[, Industry := gsub(" \\(.*\\)", "", SubIndustry)]

  return(sector_mapping)
}
#' Update VIX data in database
#'
#' @param db_path Path to SQLite database
#' @param from_date Start date for update (NULL = auto-detect)
#' @return Number of rows updated (invisible)
#' @export
update_vix_in_db <- function(db_path, from_date = NULL) {
  # Check requirements
  if (!requireNamespace("RSQLite", quietly = TRUE)) {
    stop("RSQLite required. Install with: install.packages('RSQLite')")
  }
  if (!requireNamespace("quantmod", quietly = TRUE)) {
    stop("quantmod required. Install with: install.packages('quantmod')")
  }

  # Connect to database
  db <- RSQLite::dbConnect(RSQLite::SQLite(), db_path)  # CHANGED
  on.exit(RSQLite::dbDisconnect(db))  # CHANGED

  # Get last VIX date if from_date not specified
  if (is.null(from_date)) {
    last_date <- RSQLite::dbGetQuery(db, "SELECT MAX(date) as last_date FROM stocks WHERE symbol = 'VIX'")$last_date  # CHANGED
    from_date <- if(is.na(last_date)) "2014-01-01" else as.Date(last_date) + 1
  }

  cat(sprintf("Updating VIX from %s\n", from_date))

  # Download VIX (using ^VIX)
  tryCatch({
    vix_data <- quantmod::getSymbols("^VIX",  # CHANGED
                                     src = 'yahoo',
                                     from = from_date,
                                     to = Sys.Date(),
                                     auto.assign = FALSE,
                                     warnings = FALSE)

    if (nrow(vix_data) == 0) {
      cat("No new data to update\n")
      return(invisible(0))
    }

    # Convert to database format
    vix_df <- data.table::data.table(  # CHANGED
      date = format(index(vix_data), "%Y-%m-%d"),
      symbol = "VIX",
      price = as.numeric(quantmod::Cl(vix_data)),  # CHANGED
      price_adjusted = as.numeric(quantmod::Ad(vix_data))  # CHANGED
    )

    # Insert into database
    for (i in 1:nrow(vix_df)) {
      RSQLite::dbExecute(db,  # CHANGED
                         "INSERT OR REPLACE INTO stocks (date, symbol, price, price_adjusted) VALUES (?, ?, ?, ?)",
                         params = list(vix_df$date[i], vix_df$symbol[i], vix_df$price[i], vix_df$price_adjusted[i])
      )
    }

    cat(sprintf("Updated %d days of VIX data\n", nrow(vix_df)))
    return(invisible(nrow(vix_df)))

  }, error = function(e) {
    cat("Error updating VIX:", e$message, "\n")
    return(invisible(-1))
  })
}
