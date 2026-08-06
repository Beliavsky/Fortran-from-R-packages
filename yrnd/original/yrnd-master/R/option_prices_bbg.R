#' option_prices_bbg
#'
#' @param fut_bbg_ticker a single Bloomberg ticker of a STIR or bond futures contract, in character format
#' @param date a date for the recovery of options' closing prices, in Date format
#'
#' @returns provided there are listed options on the futures, call and put options strike prices and attached call and put options closing market prices in numeric format, the maturity date of the options in Date format, the option's style in character format, the option's currency in character format, the listing place of the options in character format, the closing price of the futures contract in numeric format. When options of different maturities are listed on the same futures contract, this information concerns options whose maturity date is closest to the maturity date of the futures contract.
#' @export
#' @import dplyr
#' @import Rblpapi
#' @import tibble
#' @importFrom stats na.omit
#'
#' @examples
#' \dontrun{
#' option_prices_bbg("ERM6", as.Date("2026-03-20"))
#' }
#'

option_prices_bbg <- function(fut_bbg_ticker, date){

  blpConnect()

  if(length(fut_bbg_ticker) > 1){message("please enter one ticker at a time")
  } else {

    underlying <- bdp(paste0(fut_bbg_ticker, " Comdty"), "FUT_NOTL_BOND")

    if(underlying == ""){message("please enter a STIR futures or a bond futures Bloomberg ticker")
    } else{

      option_tickers <- lookupSecurity(paste0(fut_bbg_ticker, " Comdty"), yellowkey = "cmdt",
                                       language = "english", maxResults = 1000, verbose = F) %>%
        filter(grepl("OP", description)) %>% mutate_at("security", ~gsub("<cmdty>", " Comdty", .))

      if(nrow(option_tickers) == 0){message("no option listed on this asset")
      } else{

        extract_prices <- bdh(option_tickers$security, "PX_LAST", start.date = date, end.date = date)
        extract_prices <- do.call(rbind, extract_prices)

        if(nrow(extract_prices)  == 0){ message("please enter a trading day")
        } else if(ncol(extract_prices) == 1) { message("options prices not available anymore")
        } else {

          extract_prices <- extract_prices %>% rownames_to_column("security")

          opt_matu <- bdp(head(extract_prices$security, 1), "OPT_EXPIRE_DT") %>% unlist() %>% as.Date()
          opt_style <- bdp(head(extract_prices$security, 1), "OPT_EXER_TYP") %>% unlist()
          opt_curr <- bdp(head(extract_prices$security, 1), "CRNCY") %>% unlist()
          opt_exchange <- bdp(head(extract_prices$security, 1), "FUT_EXCH_NAME_LONG") %>% unlist()
          fut_price <- bdh(paste0(fut_bbg_ticker, " Comdty"), "PX_LAST", start.date = date,
                           end.date = date) %>% dplyr::select(PX_LAST) %>% unlist()
          names(opt_matu) <- names(opt_style) <- names(opt_curr) <- names(opt_exchange) <- names(fut_price) <- c()

          suppressWarnings({
            extract_prices <- extract_prices %>% mutate_at("security", ~gsub(" Comdty", "", .)) %>%
              mutate(option = gsub('[0-9].+', '', gsub(fut_bbg_ticker, "", security))) %>%
              mutate(strike_price = gsub(paste0(fut_bbg_ticker, option, collapse = "|"), "", security)) %>%
              mutate_at("option", ~gsub(" ", "", .)) %>% dplyr::select(-c(date, security)) %>%
              filter(!is.na(PX_LAST)) %>% mutate_at(c("PX_LAST", "strike_price"), as.numeric) %>%
              arrange(option, strike_price) %>% unique
            call <- extract_prices %>% filter(option == "C") %>% rename(call_price = PX_LAST) %>%
              dplyr::select(-option)
            put <- extract_prices %>% filter(option == "P") %>% rename(put_price = PX_LAST)  %>%
              dplyr::select(-option)
            options_prices <- left_join(call, put, by = "strike_price") %>% relocate(strike_price) %>%
              arrange(strike_price)
          })

          prices <- list(options_prices = options_prices, option_maturity_date = opt_matu,
                         option_exercice_type = opt_style, option_currency = opt_curr, option_exchange = opt_exchange,
                         future_price = fut_price)
          return(prices)}
      }
    }
  }
}
