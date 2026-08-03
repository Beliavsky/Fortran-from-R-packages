#' Coupon date payments
#'
#' @description
#' Function to calculate the upcoming coupon payment dates of a given asset,
#' based on its payment frequency. The list of payment dates encompass
#' the time period between the analysis date and the maturity of the asset.
#'
#' @param maturity Last day of the contract: YYYY-MM-DD.
#' Alternatively, it can be a numeric value that represents the duration of the contract in years.
#' @param analysis.date Date in which the asset is valued. By default, the
#' current date.
#' @param freq Frequency of payments of a given asset in a year. For LIBOR and
#' IBR the default frequency is four (quarterly payments). TES has a default
#' frequency of one (annual payments).
#' @param asset.type String that determines the asset type to value. See also
#' 'Details'.
#' @param loc String related to the location of the asset. It is used to
#' calculate the effective dates, taking into account the business days
#' of the given location. See also 'Details'.
#' @param convention String that establishes if the effective dates are
#' calculated using Following, Modified Following, Backward or Backward Following.
#' See also 'Details'.
#' @param trade.date The date on which the transaction occurs. It is used to calculate
#' maturity as a date, when given in years. Also required for non-trivial cases such as
#' bonds with long first coupon.
#' @param coupon.schedule String that establishes if a bond first coupon period is a long
#' first coupon or a short first coupon. On the contrary, establishes if last coupon period
#' is long last coupon or a short last coupon. See also 'Details'.
#'
#' @note
#' If only maturity is given, function assumes that the coupon payments have
#' already started. If \code{maturity} and \code{trade.date} are included,
#' coupon dates are calculated from \code{trade.date} to \code{maturity}.
#' If by doing so, \code{trade.date} doesn't converge to \code{maturity}, month remainder
#' is adjusted according to \code{coupon.schedule}.
#' For LIBOR assets, function adds 2 business days to \code{trade.date}.
#'
#'
#' @details
#' \code{asset.type} makes reference to the following type of assets:
#' \itemize{
#'    \item "TES" for Colombian Treasury Bonds (default).
#'    \item "FixedIncome" for assets that are indexed to a fixed income with
#'    different frequency of payments.
#'    \item "IBR" for bonds and assets indexed to 3M IBR rate.
#'    \item "IBRSwaps" for swaps indexed to IBR rate.
#'    \item "LIBOR" for bonds and assets indexed to 3M LIBOR.
#'    \item "UVRSwaps" for cross-currency swaps indexed to UVR-IBR rate.
#'    \item "LIBORSwaps" for Interest Rate Swaps (IRS) indexed to 3M LIBOR.
#'      }
#'
#' @details
#' \code{loc} makes reference to the following locations:
#' \itemize{
#'    \item "BOG" for colombian issued assets and national business days
#'    (default).
#'    \item "LDN" for business days of London.
#'    \item "NY" for business days of New York
#'    \item "NYLDN" for the intersection of business days in New York and
#'    London.
#'    \item "BOGNY" for the intersection of business days in Bogota and
#'    New York.
#'      }
#'
#'@details
#' \code{convention} makes reference to the following type of business day conventions:
#' \itemize{
#'    \item "F" for Following business day convention.
#'    \item "MF" for Modified Following business day convention.
#'    \item "B" for Backward business day convention.
#'    \item "MB" for Modified Backward business day convention.
#'      }
#'
#'@details
#' \code{coupon.schedule} makes reference to the following type of coupon payment schedule
#' of a bond:
#' \itemize{
#'    \item "LF" for Long First coupon payment.
#'    \item "LL" for Long Last coupon payment.
#'    \item "SF" for Short First coupon payment.
#'    \item "SL" for Short Last coupon payment.
#'      }
#'
#' @return Upcoming coupon dates and dates of payment according to business day conventions.
#' @export
#'
#' @examples
#' coupon.dates("2028-04-03")
#' coupon.dates(maturity = 2, analysis.date = "2021-10-01")
#' coupon.dates(maturity = "2029-10-01", asset.type = "FixedIncome", freq = 2, convention = "MB")
#' coupon.dates(maturity = "2028-02-29", analysis.date = "2022-07-29", trade.date = "2022-07-29",
#'              asset.type = "TES", coupon.schedule = "SF")
#' coupon.dates(maturity = "2025-11-30", analysis.date = "2022-03-01", trade.date = "2021-05-31",
#'              asset.type = "IBR", loc = "NY", convention = "F")

coupon.dates <- function(maturity, analysis.date = Sys.Date(), asset.type = "TES",
                         freq = NULL, convention = "F", loc = "BOG",
                         trade.date = NULL, coupon.schedule = "SF") {

  # Parameter validation ----------------------------------------------------
  # Date format: Input maturity and date of analysis must have a proper date format.
  if (!lubridate::is.Date(analysis.date)) {
    try(
      analysis.date <- as.Date(analysis.date),
      stop(paste0(deparse(sys.call()), ":", analysis.date," is not valid as Date."),
           call. = FALSE)
    )
  }

  if (is.numeric(maturity)) {
    m.date <- FALSE
  } else if (is.character(maturity) && !is.na(as.Date(maturity, format = "%Y-%m-%d"))) {
    m.date <- TRUE
  } else if (lubridate::is.Date(maturity)) {
    m.date <- TRUE
  } else {
    stop(paste0(deparse(sys.call()), ":", maturity, " is not valid, has to be a Date or numeric value"),
         call. = FALSE)
  }


  if (m.date == TRUE) {

    # If trade.date input is introduced, then coupon dates are determined by starting from trade.date,
    # Otherwise, coupon dates are determined by starting from maturity and going backwards, therefore,
    # this differentiation has to be made with input, code will change afterwards depending on input value.

      if (!is.null(trade.date)) {
        trade.date <- as.Date(trade.date)
        input <- 2
      } else {
        trade.date <- "1970-01-01"
        input <- 1
      }


      # Establishes settlement date, two business days after trade date, only for LIBOR.
      if (asset.type %in% c("LIBOR", "LIBORSwaps")) {
      settlement.date <- quantdates::AddBusinessDays(
        date    = trade.date,
        numDate = 2,
        loc     = loc
      )
      } else {
        settlement.date <- trade.date
      }



  } else if (m.date == FALSE) {
    input <- 2
    if (!is.null(trade.date)) {
      trade.date <- as.Date(trade.date)
    } else {
      trade.date <- analysis.date
    }

    # Establishes settlement date, two business days after trade date.

    if (asset.type %in% c("LIBOR", "LIBORSwaps")) {
      settlement.date <- quantdates::AddBusinessDays(
        date    = trade.date,
        numDate = 2,
        loc     = loc
      )
    } else {
      settlement.date <- trade.date
    }


    # Establishes maturity as a date, from numeric value in years.
    maturity <- lubridate::`%m+%`(
      quantdates::AddBusinessDays(
        date    = settlement.date,
        numDate = 0,
        loc     = loc
      ),
      months(round(maturity * 12,0))
    )
  }

  # Function ----------------------------------------------------------------
  ## Sequence ---------------------------------------------------------------

  res <- list()

  if (length(freq) == 0) {
    if (asset.type == "TES") {
      freq <- 1
    } else if (asset.type %in% c("LIBOR", "IBR","LIBORSwaps","IBRSwaps", "UVRSwaps")) {
      freq <- 4
    }
  } else {
    if (freq == 0) {
      dates <- maturity
    } else if (freq != 1) {
      freq <- freq
    }
  }

  if (freq == 0) {
    freq <- freq
  } else if (!12 %% freq == 0) {
    stop("frequency of payments has to be a factor of 12")
  }

  # Creates a sequence of dates by taking the maturity and going backwards by
  # the desired annual frequency until reaching the analysis date

  if (input == 1) {
  total.months <- lubridate::interval(settlement.date, maturity) %/% months(1)
  remainder    <- ((total.months) %% (12 / freq))
  endpoint     <- as.numeric(total.months - remainder)
  int.seq      <- seq(0, endpoint, by = (12 / freq))

  dates <- rev(lubridate::`%m-%`(as.Date(maturity), months(int.seq)))

  # Creates a sequence of dates by taking the settlement date and going forward by
  # the desired annual frequency until reaching maturity.

  } else if (input == 2) {
    total.months <- lubridate::interval(settlement.date, maturity) %/% months(1)
    remainder    <- ((total.months) %% (12 / freq))

    if (remainder == 0) {
      coupon.schedule <- "LF"
    }

    if (coupon.schedule == "SF") {
      int.seq      <- seq(remainder, total.months,
                          by = (12 / freq))
    } else if (coupon.schedule == "LF") {
      int.seq      <- seq(remainder + 12 / freq, total.months,
                          by = (12 / freq))
    } else if (coupon.schedule == "SL") {
      int.seq      <- c(seq(12 / freq, total.months - remainder,
                            by = (12 / freq)), total.months)
    } else if (coupon.schedule == "LL") {
      int.seq      <- c(seq(12 / freq, total.months - remainder - 12 / freq,
                            by = (12 / freq)), total.months)
    }

    dates <- (lubridate::`%m+%`(as.Date(settlement.date), months(int.seq)))

  }
  dates <- dates[which(dates > analysis.date)]
  #Coupon dates must be at least, a day after settlement date.
  dates <- dates[which(dates > settlement.date )]


  #Settlement date cannot belong to coupon dates
  res$dates <- as.Date(dates)

  ## Effective payment dates ------------------------------------------------
  effective.dates <- c()

  if (convention == "F") {
    # Gives the actual business dates of the coupon dates, by returning the
    # next business date if the date is a holiday or weekend.
    for (i in seq_along(dates)) {
      effective.dates[i] <- quantdates::AddBusinessDays(date = dates[i], numDate = 0, loc = loc)
      effective.dates    <- as.Date(effective.dates, origin = "1970-01-01")
    }
  } else if (convention == "B") {
    # Gives the actual business dates of the coupon dates, by returning the
    # previous business date if the date is a holiday or weekend.
    for (i in seq_along(dates)) {
      business.date <- ifelse(dates[i] == quantdates::AddBusinessDays(date = dates[i], numDate = 0, loc = loc),
                              0,-1)
      effective.dates[i] <- quantdates::AddBusinessDays(date = dates[i], numDate = business.date, loc = loc)
      effective.dates    <- as.Date(effective.dates, origin = "1970-01-01")
    }
  } else if (convention == "MF") {

  for (i in seq_along(dates)) {
    effective.dates[i] <- quantdates::AddBusinessDays(date = dates[i], numDate = 0, loc = loc)
    effective.dates    <- as.Date(effective.dates, origin = "1970-01-01")
  }
  # Implements modified following. In the case of LIBOR assets
  # returns the last business date of the month.
  for (i in seq_along(effective.dates)) {
    effective.dates[i] <- if (effective.dates[i] <= quantdates::LastDayOfMonth(date = dates[i])) {
      effective.dates[i]
    } else {
      effective.dates[i] <- quantdates::AddBusinessDays(date = effective.dates[i], numDate = (-1), loc = loc)
    }
  }
  } else if (convention == "MB") {
    for (i in seq_along(dates)) {
      business.date <- ifelse(dates[i] == quantdates::AddBusinessDays(date = dates[i], numDate = 0, loc = loc),
                              0,-1)
      effective.dates[i] <- quantdates::AddBusinessDays(date = dates[i], numDate = business.date, loc = loc)
      effective.dates    <- as.Date(effective.dates, origin = "1970-01-01")
    }

    for (i in seq_along(effective.dates)) {
      effective.dates[i] <- if (effective.dates[i] >= lubridate::floor_date(dates[i])) {
        effective.dates[i]
      } else {
        effective.dates[i] <- quantdates::AddBusinessDays(date = effective.dates[i], numDate = (1), loc = loc)
      }
    }
  }

  res$effective.dates <- as.Date(effective.dates, origin = "1970-01-01")

  return(res)
}

#' Coupon payment calculation
#'
#' @description
#' Function that returns coupon values according to specified payment dates and
#' a day count convention. Yields the values of cash flows for the remaining
#' duration of assets, following a date payment structure, face value -or
#' principal- and a specified coupon rate.
#'
#' @param dates Coupon payment dates.
#' @param coupon.rate Coupon rate of the asset. Can be an unique numeric
#' value or a vector corresponding to each coupon payment date.
#' @param principal Notional amount for the asset.
#' @param daycount Day count convention. See also 'Details'.
#' @param maturity Only necessary in cases where coupon payment dates are not provided in the \code{dates} parameter.
#' Last day of the contract. Can be a numeric value that represents the duration of the contract in years.
#' @inheritParams coupon.dates
#'
#' @details
#' \code{asset.type} makes reference to the following type of assets:
#' \itemize{
#'    \item "TES" for Colombian Treasury Bonds (default).
#'    \item "FixedIncome" for assets that are indexed to a fixed income with
#'    different frequency of payments.
#'    \item "IBR" for bonds and assets indexed to 3M IBR rate.
#'    \item "IBRSwaps" for swaps indexed to IBR rate.
#'    \item "LIBOR" for bonds and assets indexed to 3M LIBOR.
#'    \item "UVRSwaps" for cross-currency swaps indexed to UVR-IBR rate.
#'    \item "LIBORSwaps" for Interest Rate Swaps (IRS) indexed to 3M LIBOR.
#'      }
#' @details
#' \code{daycount} convention accepts the following values:
#' \itemize{
#'     \item 30/360.
#'     \item ACT/365.
#'     \item ACT/360 (Default).
#'     \item ACT/365L.
#'     \item NL/365.
#'     \item ACT/ACT-ISDA
#'     \item ACT/ACT-AFB
#' }
#'
#'@details
#' \code{coupon.schedule} makes reference to the following type of coupon payment schedule
#' of a bond:
#' \itemize{
#'    \item "LF" for Long First coupon payment.
#'    \item "LL" for Long Last coupon payment.
#'    \item "SF" for Short First coupon payment.
#'    \item "SL" for Short Last coupon payment.
#'      }
#'
#' @return $Coupons or $Cash flows of the asset in analysis.
#' @export
#'
#' @examples
#' coupons(dates = c("2020-09-10", "2020-12-10", "2021-03-10"),
#'         coupon.rate = 0.06)
#' coupons(dates = c("2020-09-10", "2020-12-10", "2021-03-10"),
#'         coupon.rate = 0.08, principal = 1000,
#'         asset.type = "LIBOR", daycount = "ACT/365")
#' coupons(dates = c("2020-09-10","2020-12-10", "2021-03-10"),
#'         coupon.rate = 0.07, asset.type = "FixedIncome",
#'         freq = 4, daycount = "NL/365")
#' coupons(coupon.rate = c(0.04,0.04,0.42,0.045,0.05),
#'         maturity = "2024-01-05", analysis.date = "2023-01-03",
#'         asset.type = "IBR")
#' coupons(coupon.rate = 0.03, maturity = 1.08,
#'         analysis.date = "2020-02-29", trade.date = "2020-02-29",
#'         asset.type = "IBR", coupon.schedule = "LF")
#'
coupons <- function(dates = NULL, coupon.rate, principal = 1, asset.type = "TES",
                    freq = NULL, daycount = "ACT/360", loc = "BOG", maturity = NULL,
                    analysis.date = Sys.Date(), trade.date = NULL, coupon.schedule = "SF") {

  # Parameter validation ----------------------------------------------------
  if (is.null(dates)) {
    if (is.null(trade.date)) {
    dates <- coupon.dates(
      maturity = maturity, asset.type = asset.type, freq = freq,
      analysis.date = analysis.date, convention = "F",
      loc = loc
    )$dates
    } else if (!is.null(trade.date)) {
      dates <- coupon.dates(
        maturity = maturity, asset.type = asset.type, freq = freq,
        analysis.date = trade.date, convention = "F",
        loc = loc, trade.date = trade.date, coupon.schedule = coupon.schedule
      )$dates
    }
  } else {
    try(
      dates <- lubridate::ymd(dates),
      stop(paste0(deparse(sys.call()), ":", dates, " is not valid as Date."),
           call. = FALSE)
    )
  }

  if (asset.type %in% c("TES", "FixedIncome")) {
    if (length(coupon.rate) != 1) {
      stop("coupon.rate has to be a unique fixed rate")
    }
  }

  if (!length(coupon.rate) == 1) {
    if (!length(coupon.rate) == length(dates)) {
    stop("length of coupon.rate does not match with length of dates")
    }
  }

  if (length(freq) == 0) {
    if (asset.type == "TES") {
      freq <- 1
    } else if (asset.type %in% c("LIBOR", "IBR","LIBORSwaps","IBRSwaps", "UVRSwaps")) {
      freq <- 4
    }
  }

  if (freq == 0) {
    freq <- freq
  } else if (!12 %% freq == 0) {
    stop("frequency of payments has to be a factor of 12")
  }

  # Function ----------------------------------------------------------------
  if (asset.type == "TES") {
    # For TES bonds, makes a sequence of fixed coupon payments until maturity,
    # with last payment being face value plus a coupon.
    CashFlows <- rep(principal * coupon.rate, length(dates))
    CashFlows[length(dates)] <- CashFlows[length(dates)] + principal
  } else if (asset.type %in% c("LIBOR", "IBR", "LIBORSwaps", "UVRSwaps",
                               "FixedIncome", "IBRSwaps","CCS")) {
    if (is.null(trade.date)) {
      if (asset.type %in% c("LIBOR", "IBR", "LIBORSwaps",
                          "UVRSwaps", "IBRSwaps")) {
      # Return the most recent coupon payment. Important as analysis date could
      # be different from next coupon payment. Takes the first coupon payment
      # in the analysis time frame and finds the previous coupon payment.
      # The process is repeated 3 times, since there is the possibility that dates[1],
      # drops on 30 of April, outputs 30 Jan but actually previous date is 31 January,
      # but there is not 31 of April. This intricate cases are solved by using three dates,
      # when freq = 4 and for freq = x it is required a total of 4 backwards operations.
      # max is used to obtain 31 of Jan instead of 30 Jan.
        if (freq == 0) {
          analysis.date <- analysis.date
        } else {
      analysis.date  <- lubridate::`%m-%`(dates[1], months(3))
      analysis.date2 <- lubridate::`%m-%`(dates[2], months(6))
      analysis.date3 <- lubridate::`%m-%`(dates[3], months(9))
      analysis.date  <- max(analysis.date, analysis.date2, analysis.date3, na.rm = TRUE)
        }
    } else if (asset.type %in% c("FixedIncome","CCS")) {
      if (freq == 0) {
        analysis.date <- analysis.date
        } else {
      analysis.date  <- lubridate::`%m-%`(dates[1], months(12 / freq))
      analysis.date2 <- lubridate::`%m-%`(dates[2], months(24 / freq))
      analysis.date3 <- lubridate::`%m-%`(dates[3], months(36 / freq))
      analysis.date4 <- lubridate::`%m-%`(dates[4], months(48 / freq))
      analysis.date <- max(analysis.date, analysis.date2, analysis.date3, analysis.date4, na.rm = TRUE)
        }
      }
  } else if (!is.null(trade.date)) {
    dates <- c(trade.date, as.character(dates))
    analysis.date <- dates[(which(dates > analysis.date)) - 1][1]
    dates <- dates[which(dates > analysis.date)]
  }

    # Gives the coupon payment percentages according to the time between
    # payments.
    PaymentPercentages <- apply(
      X      = data.frame(c(analysis.date, as.character(dates[-length(dates)])), dates),
      MARGIN = 1,
      function(x) {
        quantdates::day_count(
          tinitial   = x[1],
          tfinal     = x[2],
          convention = daycount
        )
      }
    )

    # Finds the cash flows of the assets by taking the sequence of coupons or
    # payments given by: i. The rate of return, or ii. The coupon rate, or
    # iii. The market rate times the faces value of the asset. This is then
    # adjusted by the payment percentage given by the time difference.

    CashFlows <- coupon.rate * principal * PaymentPercentages
    CashFlows[length(CashFlows)] <- CashFlows[length(CashFlows)] + principal
}

  return(CashFlows)
}

#' Discount factors
#'
#' @description
#' Function that calculates discount factors given effective payment dates
#' and a discount rate. Optional parameters available to calculate discrete
#' or continuous discount factors.
#'
#' @param rates Discount rates given by the zero coupon rate curve.
#' Can also be a unique discount rate.
#' @param rate.type (1) for discrete compounded discount rates and (0) for continuosly
#' compounded discount rates. By default rates are assumed to be discrete.
#' @inheritParams coupon.dates
#' @inheritParams coupons
#'
#' @return Discount factors.
#' @export
#'
#' @examples
#' discount.factors(dates = c("2020-09-10", "2020-12-10", "2021-03-10"), rates = c(0.07, 0.075, 0.08),
#'                  analysis.date = "2010-09-01")
#' discount.factors(dates = c("2025-09-01", "2025-12-01", "2026-03-01", "2026-06-01"),
#'                  rates = c(0.01, 0.015, 0.017, 0.02), analysis.date = "2025-06-01",
#'                  rate.type = 1, freq = 4)
#'
discount.factors <- function(dates, rates, analysis.date = Sys.Date(),
                             rate.type = 1, freq = 1) {
  DiscountFactors <- c()

  # Parameter validation ----------------------------------------------------
  if (length(rates) == 1) {
    rates <- rep(rates, length(dates))
  } else if (length(dates) > length(rates)) {
    stop("Missing rates for all coupon payments")
  }

  if (!is.null(dates)) {
    dates <- lubridate::ymd(dates)
  }

  # Function ----------------------------------------------------------------
  for (i in seq_along(dates)) {
    if (rate.type == 1) {
      # Calculates the discount factor for each time frame given by the analysis
      # date and the input dates. Follows a present value formula in discrete
      # time.
      DiscountFactors[i] <- (1 / (1 + (as.numeric(rates[i])) / freq))^(discount.time(
        tinitial = analysis.date,
        tfinal   = lubridate::ymd(dates[i])
      ) * freq)
    } else if (rate.type == 0) {
      # Calculates the discount factor for each time frame given by the analysis
      # date and the input dates. Follows a present value formula in continuous
      # time.
      DiscountFactors[i] <- exp(-(as.numeric(rates[i]) * discount.time(
        tinitial = analysis.date,
        tfinal   = dates[i]
      )))
    }
  }
  return(DiscountFactors)
}

#' Quantil's discount time convention
#'
#' @description Function to count the number of years between dates according
#' to Quantil's discount convention. A year is defined as the difference in one year,
#' between two dates with the exact month and day. Meanwhile, partial years are defined as the
#' quotient between the number of elapsed days within a year and the total number of days that
#' make up that year. Total number of years between the two dates is then the sum between complete full
#' years and the partial portion.
#'
#' @param tinitial Initial date of analysis.
#' @param tfinal Final date of analysis.
#'
#' @details There is an exception. For example, for initial date 29-February, a year is defined as the
#' 28 of February of the next year. Meanwhile four years, is defined as 29 of February four years after.
#'
#'
#'
#' @return Number of years between the specified dates.
#' @export
#'
#' @examples
#' discount.time(tinitial = "2024-07-13", tfinal = "2025-03-01")
#' discount.time(tinitial = "2024-02-29", tfinal = "2025-02-28")
#' discount.time(tinitial = "2024-02-29", tfinal = "2028-02-29")
#'
discount.time <- function(tinitial, tfinal) {
  # Parameter validation ----------------------------------------------------
  # Date format: Input initial and final date must have a proper date format.
  if (!lubridate::is.Date(tinitial)) {
    try(
      tinitial <- as.Date(tinitial),
      stop(paste0(deparse(sys.call()), ":", tinitial, " is not valid as Date."),
           call. = FALSE)
    )
  }

  if (!lubridate::is.Date(tfinal)) {
    try(
      tfinal <- as.Date(tfinal),
      stop(paste0(deparse(sys.call()), ":", tfinal, " is not valid as Date."),
           call. = FALSE)
    )
  }

  # Function ----------------------------------------------------------------
  # Create a sequence of days between input days. A year is defined as two dates
  # with the same day, month but a difference of a year. Counts the integer years,
  # defined by this rule, and the remainder is calculated as the number of days between
  # last date counted, divided by the amount of days that make up a year, starting from
  # last date counted.
  if(tinitial > tfinal) {
    stop("The analysis date is after payments.")
  } else {

    fullyears <- length(seq(from = tinitial, to = tfinal, by = "year"))-1

    next.year <- lubridate::`%m+%`(tinitial,
                                   months(12*fullyears))

    final.year <- lubridate::`%m+%`(next.year,
                                    months(12))

    x.num <- seq(from = next.year, to = tfinal, by = "day")
    x.denom <- seq(from = next.year, to = final.year, by = "day")

    return (fullyears + (length(x.num)-1)/(length(x.denom)-1))

  }
}

#' Accrued interest
#'
#' @description Calculates the accumulated coupon or accrued interests of the
#' asset, from its last coupon or cash flow payment.
#'
#' @inheritParams coupon.dates
#' @inheritParams coupons
#'
#' @details
#' \code{asset.type} makes reference to the following type of assets:
#' \itemize{
#'    \item "TES" for Colombian Treasury Bonds (default).
#'    \item "FixedIncome" for assets that are indexed to a fixed income with
#'    different frequency of payments.
#'    \item "IBR" for bonds and assets indexed to 3M IBR rate.
#'    \item "IBRSwaps" for swaps indexed to IBR rate.
#'    \item "LIBOR" for bonds and assets indexed to 3M LIBOR.
#'    \item "UVRSwaps" for cross-currency swaps indexed to UVR-IBR rate.
#'    \item "LIBORSwaps" for Interest Rate Swaps (IRS) indexed to 3M LIBOR.
#'      }
#'
#' \code{daycount} convention accepts the following values:
#' \itemize{
#'     \item 30/360.
#'     \item ACT/365.
#'     \item ACT/360 (Default).
#'     \item ACT/365L.
#'     \item NL/365.
#'     \item ACT/ACT-ISDA
#'     \item ACT/ACT-AFB
#' }
#'
#' @return Accrued interest of the bond from the last coupon payment to the
#' liquidation (valuation date).
#' @export
#'
#' @examples
#' accrued.interests(coupon.rate = 0.04, maturity = '2029-08-10',
#'                   asset.type = 'LIBOR', daycount = "30/360")
#' accrued.interests(coupon.rate = 0.04, maturity = '2029-08-10',
#'                   daycount = "NL/365")
#' accrued.interests(coupon.rate = 0.04, maturity = '2029-08-10',
#'                   asset.type= 'IBR', daycount = "ACT/360")
#' accrued.interests(coupon.rate = 0.04, maturity = '2029-08-10', freq= 2,
#'                   asset.type= 'FixedIncome', daycount = "ACT/365")
#'
#'
accrued.interests <- function(maturity, analysis.date = Sys.Date(),
                              coupon.rate, principal = 1, asset.type = "TES",
                              freq = NULL, daycount =  "ACT/360") {
  # Parameter validation ----------------------------------------------------
  # The asset may not have expired already.
  if (analysis.date >= maturity) {
    return(0)
  }

  if (length(freq) == 0) {
    if (asset.type == "TES") {
      freq <- 1
    } else if (asset.type %in% c("LIBOR", "IBR")) {
      freq <- 4
    }
  } else {
    if (freq == 0) {
      dates <- maturity
    } else if (freq != 1) {
      freq <- freq
    }
  }

  # Calculate coupon dates of an asset.
  dates <- coupon.dates(
    maturity      = maturity,
    asset.type    = asset.type,
    freq          = freq,
    analysis.date = analysis.date
  )$dates

  # Calculate the time difference between the first coupon date and the
  # analysis date.
  Delta <- as.numeric(quantdates::day_count(
    tfinal     = dates[1],
    tinitial   = analysis.date,
    convention = daycount
  ))

  # Calculate the time difference between the previous coupon date and the
  # next coupon date.
  previous.date <- lubridate::`%m-%`(dates[1], months(12 / freq))
  Big.Delta <- as.numeric(quantdates::day_count(
    tfinal     = dates[1],
    tinitial   = previous.date,
    convention = daycount
  ))

  if (analysis.date < previous.date) {
    return(0)
  }

  # Function ----------------------------------------------------------------
  # Calculate accrued interests.
  if (analysis.date == dates[1]) {
    accrued <- 0
  } else {
    accrued <- (1-(Delta/Big.Delta)) * coupon.rate * principal*Big.Delta
    accrued <- accrued[1]
  }

  return(accrued)
}
