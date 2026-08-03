#' Build a data frame containing interest rates for CDS pricing
#' 
#' \code{build_rates} can create a data frame of interest rates from a start 
#' date to an end date, which are specified by the user. The interest rates will
#' be later used in CDS pricing. \code{build_rates} also builds the rates.RData 
#' in the package. 
#' 
#' @param start.date is the start date of the data frame; it is the earliest CDS
#'   pricing date that the user concerns.
#' @param end.date is the end date of the data frame; it is the latest CDS 
#'   pricing date that the user concerns.
#'   
#' @return a data frame that contains the CDS pricing date, the currency, the 
#'   interest rate expiry and the interest rate.
#'   
#' @details
#' Since the Markit website has the interest rates back to the 1990s, download_markit 
#' is responsible for building up all the USD interest rate data frame; for EUR 
#' and JPY, markit can only get from 2005-01-05 to now. But the biggest 
#' advantage of using download_markit is that, since the Markit website only lists 
#' the rate expiries that are actually used for CDS pricing, download_markit can
#' get the exact type of expiries of rates needed to price CDS. Also, it has 
#' expiry over 1Y to 30Y. Then we want to get data for EUR and JPY before 2005-01-04 (the limit of
#' markit). 
#' 
#' Another note is that, the rates on Markit have 
#' not been adjusted to the previous business day. In other words, the rates 
#' from both website is the exact rate on that day, rather than on the previous 
#' business day. But download_markit has adjusted the days 
#' already for the convenience of CDS pricing, so we don't have to worry here.
#' 
#' @seealso \link{download_markit} \link{rates}
#'   
#' @references \url{https://www.spglobal.com/en} 
#'

build_rates <- function(start.date, end.date){
  
  ## rid NOTE
  
  end <- NULL
  
  
  stopifnot(inherits(start.date, "Date") & inherits(end.date, "Date"))
  
  ## start date should be earlier than or at least same with end date
  
  stopifnot(end.date >= start.date)
  
  ## If start date is before and end date is after 2005-01-05,
  ## then we have to use markit for USD, JPY and EUR
  
  if(start.date < as.Date("2005-01-05") & end.date >= as.Date("2005-01-05")){
    
    ## First, markit gets USD rates.
    
    x.markit.USD <- download_markit(start = start.date, end = end.date, 
                                    currency = "USD")
    
    ## Then markit gets rate for EUR and JPY after 2005-01-05
    
    x.markit.EUR <- download_markit(start = as.Date("2005-01-05"), 
                                    end = end.date, currency = "EUR")
    x.markit.JPY <- download_markit(start = as.Date("2005-01-05"), 
                                    end = end.date, currency = "JPY")
    
    ## bind the data together
    
    x <- rbind(x.markit.USD, x.markit.EUR, x.markit.JPY)
  } 
  
  ## If start date  and end date are after 2005-01-05, then we only have to use
  ## markit
  
  else if(start.date >= as.Date("2005-01-05")){
    x.markit.USD <- download_markit(start = start.date, end = end.date, 
                                    currency = "USD")
    x.markit.EUR <- download_markit(start = start.date, end = end.date, 
                                    currency = "EUR")
    x.markit.JPY <- download_markit(start = start.date, end = end.date, 
                                    currency = "JPY")
    x <- rbind(x.markit.USD, x.markit.EUR, x.markit.JPY)
  } 
  
  ## If start date  and end date are before 2005-01-05, then we have to use
  ## markit for USD
  
  else if(end < as.Date("2005-01-05")){
    x.markit.USD <- download_markit(start = start.date, end = end.date, 
                                    currency = "USD")
    x <- rbind(x.markit.USD)
  } 
  
  ## sort the data as required: currently ascending date and ascending currency.
  ## Ought to sort by ascending expiry, but that is harder than it looks.
  
  x <- x[order(x$date, x$currency, decreasing = FALSE), ]
  
  return(x)
}