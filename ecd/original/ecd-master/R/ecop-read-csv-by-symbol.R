#' Read option data csv
#'
#' Read option data csv into dataframe. The dataframe is enriched with
#' \code{Date, expiration_date, days}.
#'
#' @param symbol character, option data symbol
#' @param extdata_dir optionally specify user's own extdata folder
#' @param df dataframe, it is assumed to be in CBOE heading format
#'
#' @return dataframe
#'
#' @keywords data
#'
#' @author Stephen H-T. Lihn
#'
#' @export ecop.read_csv_by_symbol
#' @export ecop.enrich_option_df
#'
#' @examples
#'
#' df <- ecop.read_csv_by_symbol("spxoption2")
#'
### <======================================================================>
"ecop.read_csv_by_symbol" <- function(symbol, extdata_dir=NULL)
{
    df <- ecd.read_csv_by_symbol(symbol, extdata_dir=extdata_dir)
    ecop.enrich_option_df(df)
}
### <---------------------------------------------------------------------->
#' @rdname ecop.read_csv_by_symbol
"ecop.enrich_option_df" <- function(df)
{
    date_format = "%Y%m%d"
    df$Date <- as.Date(as.character(df[,"TRADE_DT"]), date_format)
    df$expiration_date <- as.Date(as.character(df[,"EXPR_DT"]), date_format)
    df$days <- as.numeric(difftime(df$expiration_date, df$Date, units="days"))
    df
}
### <---------------------------------------------------------------------->
