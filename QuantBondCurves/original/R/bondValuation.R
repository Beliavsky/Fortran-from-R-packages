#' Bond valuation
#'
#' @description Function that values various asset types with varying payment frequencies.
#' It covers fixed-coupon assets, spread income assets, floating notes and fixed legs of
#' interest rate swaps.
#'
#' @inheritParams coupon.dates
#' @inheritParams coupons
#' @inheritParams discount.factors
#' @param rate.type (1) for annual compounded discount rates and (0) for continuosly
#' compounded discount rates. By default, rates are assumed to be the former.
#' @param spread Decimal value of spread added to coupon payment rate. By
#' default, \code{0}.
#' @param dirty Numeric value to determine if the calculated price is dirty or
#' clean. To calculate dirty price, set \code{dirty = 1}. Otherwise,
#' \code{dirty = 0}.
#' @param spread.only Logical condition that establishes if the output should just include
#' the spread or the complete bond value. By default, \code{FALSE}, referring to the
#' output being bond value.
#'
#' @details
#' \code{asset.type} makes reference to the following type of assets:
#' \itemize{
#'    \item "TES" for Colombian Treasury Bonds (default).
#'    \item "FixedIncome" for assets that are indexed to a fixed income with
#'    different frequency of payments.
#'    \item "IBR" for bonds and assets indexed to 3M IBR rate.
#'    \item "LIBOR" for bonds and assets indexed to 3M LIBOR.
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
#'     }
#'
#' @details
#' \code{convention} makes reference to the following type of business day conventions:
#' \itemize{
#'    \item "F" for Following business day convention.
#'    \item "MF" for Modified Following business day convention.
#'    \item "B" for Backward business day convention.
#'    \item "MB" for Modified Backward business day convention.
#'      }
#'
#' @details
#' \code{coupon.schedule} makes reference to the following type of coupon payment schedule
#' of a bond:
#' \itemize{
#'    \item "LF" for Long First coupon payment.
#'    \item "LL" for Long Last coupon payment.
#'    \item "SF" for Short First coupon payment.
#'    \item "SL" for Short Last coupon payment.
#'      }
#'
#' @return Bond value.
#' @importFrom stats approx optim stepfun runif
#' @export
#'
#' @examples
#' valuation.bonds(maturity = "2026-06-01", coupon.rate = 0.06, rates = 0.08,
#'                 analysis.date = "2022-06-01")
#' valuation.bonds(maturity = "2026-06-01", coupon.rate = 0.06, rates = rep(0.08,4),
#'                 analysis.date = "2022-06-01", rate.type = 0)
#' valuation.bonds(maturity = "2026-06-01", analysis.date= "2025-02-27",
#'                 coupon.rate = c(0.06, 0.062, 0.063, 0.065, 0.066, 0.068),
#'                 rates = c(0.08, 0.082, 0.078, 0.09, 0.077, 0.085),
#'                 asset.type = "IBR")
#' valuation.bonds(maturity = "2026-06-01", coupon.rate = 0.06,
#'                 rates = 0.08, analysis.date = "2022-06-01", asset.type = "IBR",
#'                 freq = 4, spread = 0.03)
#' valuation.bonds(maturity = 4.58, coupon.rate = 0.1256, rates = seq(0.05, 0.14, by = 0.005),
#'                 analysis.date = "2019-07-14", asset.type = "FixedIncome", freq = 4,
#'                 principal = 567, daycount = "ACT/360", rate.type = 0, trade.date = "2019-07-14",
#'                 coupon.schedule = "LL")
valuation.bonds <- function(maturity, coupon.rate, rates, principal = 1,
                            analysis.date = Sys.Date(), asset.type = "TES",
                            freq = NULL, rate.type = 1, spread = 0,
                            daycount = "NL/365", dirty = 1, convention = "F", trade.date = NULL,
                            coupon.schedule = "SF", spread.only = FALSE) {


  # Parameter validation ----------------------------------------------------

  if (length(freq) == 0) {
    if (asset.type == "TES") {
      freq <- 1
    } else if (asset.type %in% c("LIBOR", "IBR","LIBORSwaps","IBRSwaps", "UVRSwaps")) {
      freq <- 4
    }
  }

  if (asset.type %in% c("TES", "FixedIncome")) {
    if (length(coupon.rate) != 1) {
      stop("coupon.rate has to be a unique fixed rate")
    }
  }

  # Function ----------------------------------------------------------------
  BondValue <- c()

  # Calculate coupon dates of an asset.
  Dates <- coupon.dates(
    maturity        = maturity,
    asset.type      =  asset.type,
    analysis.date   = analysis.date,
    freq            = freq,
    convention      = convention,
    loc             = "BOG",
    trade.date      = trade.date,
    coupon.schedule = coupon.schedule
  )

  # Calculate coupon cash flows of an asset.
  cashFlows <- coupons(
    dates           = Dates$dates,
    coupon.rate     = coupon.rate,
    principal       = principal,
    asset.type      = asset.type,
    daycount        = daycount,
    freq            = freq,
    analysis.date   = analysis.date,
    trade.date      = trade.date,
    coupon.schedule = coupon.schedule
  )

    # Returns the discount factors of the given coupon dates adjusted to
    # modified following conventions.
    DiscountFactors <- discount.factors(
      rates         = rates,
      dates         = Dates$effective.dates,
      analysis.date = analysis.date,
      rate.type     = rate.type,
      freq          = 1
    )

    BondValue <- sum(cashFlows * (DiscountFactors))

  if (dirty == 1) {
    BondValue <- BondValue
  } else if (dirty == 0) {
    # Calculate the accrued interest of the bond, given the last coupon payment.
    accrued.coupon <- accrued.interests(
      analysis.date = analysis.date,
      coupon.rate   = coupon.rate,
      maturity      = maturity,
      daycount      = daycount,
      freq          = freq,
      asset.type    = asset.type
    )

    BondValue <- BondValue - accrued.coupon

  }
  SpreadValue <- 0

  if (!is.null(spread) & spread != 0) {
    # Warning -----------------------------------------------------------------
    # It would be better to include it at coupon.rate, though output is the same both ways.
    if (asset.type %in% c("TES", "FixedIncome")) {
      warning("Spread should be included in coupon.rate with fixed rate assets")
    }

    # Calculate coupons of a bond with face value 0 and coupon payments equal
    # to the spread rate.
    CouponSpread <- coupons(
      dates           = Dates$dates,
      coupon.rate     = spread,
      principal       = principal,
      asset.type      = asset.type,
      analysis.date   = analysis.date,
      daycount        = daycount,
      freq            = freq,
      trade.date      = trade.date,
      coupon.schedule = coupon.schedule
    )

    CouponSpread[length(Dates$dates)] <- CouponSpread[length(Dates$dates)] - principal

    SpreadValue <- c()

     DiscountFactors <- discount.factors(
        rates         = rates,
        dates         = Dates$effective.dates,
        analysis.date = analysis.date,
        rate.type     = rate.type,
        freq          = 1
      )

      for (i in seq_along(Dates$dates)) {
        SpreadValue[i] <- (CouponSpread[i] * as.numeric(DiscountFactors[i]))
      }

      SpreadValue <- sum(as.numeric(SpreadValue))
    # Bond value is equivalent to a par bond plus spread bond with face value 0.
    BondValue <- BondValue + SpreadValue
  }
  # This if was added because knowing spreadValue of a bond, without taking into account principal,
  # is helpful for CrossCurrency swap valuation.
  if (spread.only == TRUE) {
    BondValue <- SpreadValue
  }
  return(BondValue)
}

#' Swap valuation
#'
#' @description Function that values Interest Rate Swaps (IRS) and Cross Currency
#' Swaps (CCS).
#'
#' @inheritParams coupon.dates
#' @inheritParams coupons
#' @inheritParams valuation.bonds
#' @param coupon.rate For (IRS), coupon rate of the fixed leg. For (CCS), coupon rate
#' of local fixed leg.
#' @param principal For (IRS), notional amount for both legs. For (CCS), notional
#' amount of local leg.
#' @param float.rate For (IRS), last observed floating rate, necessary for variable leg when swap valuation
#' date doesn't belong to a coupon date. For (CCS), last observed local floating rate. By default,
#'  \code{NULL}.
#' @param rates For (IRS) discount rates given by the zero coupon rate curve. For (CCS),
#' represents discount rates for local currency. Can be a vector that corresponds
#' to each coupon date or a curve with at least, nodes with 3 decimals.
#' @param rates2 Discount rates given by the foreign zero coupon rate curve.
#' Can be a vector that corresponds to each coupon date or a curve with at least,
#' nodes with 3 decimals.
#' @param basis.rates "Discount rates" given by the basis curve. Can be a vector
#' that corresponds to each coupon date or a curve with at least,
#' nodes with 3 decimals.
#' @param Legs For (CCS), string that establishes the type of both legs that makeup
#' the Cross Currency Swap. See also 'Details'.
#' @param coupon.rate2 Coupon rate of the foreign leg. By default, \code{NULL}.
#' @param ex.rate Exchange rate on analysis date. Format has to be local currency divided
#' by foreign currency.
#' @param principal2 Notional amount for the foreign leg. By
#' default, \code{NULL}.
#' @param spread2 Decimal value of spread added to foreign floating rate. By
#' default, \code{0}.
#' @param float.rate2 Last observed foreign floating rate. By
#' default, \code{NULL}.
#'
#' @details
#' \code{asset.type} makes reference to the following type of assets:
#' \itemize{
#'    \item "IBRSwaps" for swaps indexed to IBR rate.
#'    \item "LIBORSwaps" for Interest Rate Swaps (IRS) indexed to 3M LIBOR.
#'    \item "CCS" for cross currency swaps.
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
#'     }
#'
#' @details
#' \code{convention} makes reference to the following type of business day conventions:
#' \itemize{
#'    \item "F" for Following business day convention.
#'    \item "MF" for Modified Following business day convention.
#'    \item "B" for Backward business day convention.
#'    \item "MB" for Modified Backward business day convention.
#'      }
#'
#' @details
#' \code{coupon.schedule} makes reference to the following type of coupon payment schedule
#' of a bond:
#' \itemize{
#'    \item "LF" for Long First coupon payment.
#'    \item "LL" for Long Last coupon payment.
#'    \item "SF" for Short First coupon payment.
#'    \item "SL" for Short Last coupon payment.
#'      }
#'
#' @details
#' \code{Legs} makes reference to the following types of legs composition of the
#' cross currency swap
#' \itemize{
#'    \item "FF" for fixed leg in local currency and fixed leg in foreign currency.
#'    \item "FV" for fixed leg in local currency and variable leg in foreign currency.
#'    \item "VF" for variable leg in local currency and fixed leg in foreign currency.
#'    \item "VV" for variable leg in local currency and variable leg in foreign currency.
#'      }
#'
#' @return Swap value for Interest Rate Swap (IRS) or Cross Currency Swap (CCR).
#' @export
#'
#' @examples
#' # (IRS) ----------------------------------------------------------------------
#' valuation.swaps(maturity = "2030-06-01", analysis.date= Sys.Date(),
#'                 coupon.rate = 0.06, rates = 0.08, float.rate = 0.03)
#' valuation.swaps(maturity = "2021-03-09", analysis.date= "2018-03-09",
#'                 coupon.rate = 0.05, rates = 0.08, rate.type = 0)
#' valuation.swaps(maturity = "2026-07-01", analysis.date = "2023-01-02",
#'                 asset.type = "IBRSwaps", freq = 4,
#'                 coupon.rate = 0.04, rates = rep(0.05,14),
#'                 float.rate = 0.03, spread = 0)
#' # (CCS) ----------------------------------------------------------------------
#' # Curve Calibration for rates input
#' yield.curve <- c(0.103,0.1034,0.1092, 0.1161, 0.1233, 0.1280, 0.1310, 0.1320, 0.1325)
#' names(yield.curve) <- c(0,0.08,0.25,0.5,1,2,3,5,6)
#' nodes <- seq(0, 10, by = 0.001) # Our curve has nodes with three decimals.
#' rates <- curve.calibration (yield.curve = yield.curve, market.assets = NULL,
#'                             analysis.date = "2023-03-01", asset.type = "IBRSwaps",
#'                             freq = 4, rate.type = 0, daycount= "ACT/365",
#'                             fwd = 0, npieces = NULL, nodes = nodes, approximation = "constant")
#'
#' # Curve Calibration for basis.rates input
#' nodes  <- seq(0, 10, by = 0.001)
#' rates2 <- rates/4 # It is assumed foreign curve is proportional to local spot curve.
#' # Swaps input for calibration
#' ex.rate <- 4814
#' swaps <- rbind(c("2024-03-01", "FF", 0.07 , 0.0325, NA   , NA    , 2000 * ex.rate, 2000),
#' c("2025-03-01", "VV", NA   , NA    , 0.015, 0.0175, 2000 * ex.rate, 2000),
#' c("2026-03-01", "FF", 0.075, 0.03  , NA   ,  NA   , 500000, 5000000 / ex.rate),
#' c("2027-03-01", "VV", NA   , NA    , 0.01 , 0.015 , 5000000, 5000000 / ex.rate),
#' c("2028-03-01", "FF", 0.08 ,0.035  , NA   , NA    , 3000000, 3000000 / ex.rate),
#' c("2029-03-01", "VV", NA   , NA    , 0.01 , 0.0125, 3000000, 3000000 / ex.rate))
#' colnames(swaps) <- c("Mat"  ,"Legs", "C1" , "C2", "spread1", "spread2", "prin1", "prin2")
#' # Calibration
#' basis.rates <- basis.curve(swaps, ex.rate = 4814, analysis.date = "2023-03-01",
#'                            freq = c(2,2,2,2,1,1), rates = rates, rates2 = rates2,
#'                            rate.type = 1, npieces = NULL, obj = "Price",
#'                            Weights = NULL, nodes = nodes, approximation = "linear")
#'
#' # Valuation
#' valuation.swaps (maturity = "2024-03-01", analysis.date = "2023-03-01",
#'                  asset.type = "CCS", freq = 2, Legs = "FF", ex.rate = 4814,
#'                  coupon.rate = 0.07, coupon.rate2 = 0.0325,
#'                  rates = rates, rates2 = rates2, basis.rates = basis.rates,
#'                  float.rate = NULL, float.rate2 = NULL, spread = 0,
#'                  spread2 = 0, principal = 2000 * 4814, principal2 = 2000,
#'                  rate.type = 0, daycount = "ACT/365", loc = "BOG",
#'                  convention = "F")
#'

valuation.swaps <- function(maturity, analysis.date = Sys.Date(), asset.type = "IBRSwaps", freq = 4,
                            coupon.rate, rates,  float.rate = NULL, spread = 0, principal = 1,
                            Legs = "FF", ex.rate = NULL, basis.rates = NULL,
                            coupon.rate2 = NULL, rates2 = NULL, float.rate2 = NULL, spread2 = 0, principal2 = NULL,
                            rate.type = 1, daycount = "NL/365", loc = "BOG",
                            convention = "F", trade.date = NULL, coupon.schedule = "SF") {
  # Parameter validation ----------------------------------------------------
  if (!lubridate::is.Date(analysis.date)) {
    try(
      analysis.date <- as.Date(analysis.date),
      stop(paste0(deparse(sys.call()), ":", analysis.date, " is not valid as Date."),
           call. = FALSE)
    )
  }

  swapValue <- c()

  Dates <- coupon.dates(
    maturity        = maturity,
    asset.type      = asset.type,
    analysis.date   = analysis.date,
    freq            = freq,
    convention      = convention,
    trade.date      = trade.date,
    coupon.schedule = coupon.schedule
    )

  if (length(rates) > length(Dates$dates)) {
    Payments <- c()
      for (i in 1:(length(Dates$dates))) {
        Payments <- c(Payments, discount.time(
          tinitial   = analysis.date,
          tfinal     = Dates$dates[i]
        ))
      }
      Payments <- round(Payments,3)

      termpayment <- c()
      for (i in 1:(length(Dates$dates))) {
        termpayment <- c(termpayment, which(names(rates) == Payments[i]))
      }

      rates <- rates[termpayment]
  }

  if (length(rates2) > length(Dates$dates)) {
  rates2 <- rates2[termpayment]
  }

  if (length(basis.rates) > length(Dates$dates)) {
    basis.rates <- basis.rates[termpayment]
  }

  # Function ----------------------------------------------------------------
  if (asset.type %in% c("IBRSwaps", "UVRSwaps", "LiborSwaps")) {
  cashFlows <- coupons(
    dates           = Dates$dates,
    coupon.rate     = coupon.rate,
    principal       = principal,
    asset.type      = asset.type,
    daycount        = daycount,
    analysis.date   = analysis.date,
    trade.date      = trade.date,
    coupon.schedule = coupon.schedule
  )

    DiscountFactors <- discount.factors(
      rates         = rates,
      dates         = Dates$effective.dates,
      analysis.date = analysis.date,
      rate.type     = rate.type
    )

    for (i in seq_along(Dates$dates)) {
      swapValue[i] <- (cashFlows[i] * as.numeric(DiscountFactors[i]))
    }

    swapValue <- sum(as.numeric(swapValue))

  # Determines if analysis date is belongs to a coupon date, if yes, then next coupon is
  # valued almost par (depends if daycount of coupon and df match) by defining next coupon
  # rate equal to discount rate. The rest of the coupons, are always valued par on first (next)
  # coupon date and brought back to analysis date.

  previous.date <- lubridate::`%m-%`(Dates$dates[1], months(3))
  previous.date2 <- lubridate::`%m-%`(Dates$dates[2], months(6))
  previous.date3 <- lubridate::`%m-%`(Dates$dates[3], months(9))
  previous.date <- max(previous.date, previous.date2, previous.date3, na.rm = TRUE)


  if (analysis.date == previous.date) {
    float.rate <- rates[1]
  } else if (is.null(float.rate)) {
    stop("float rate is required")
  }

    cashFlows.next <- coupons(
      dates           = Dates$dates,
      coupon.rate     = float.rate,
      principal       = principal,
      asset.type      = asset.type,
      freq            = freq,
      daycount        = daycount,
      trade.date      = trade.date,
      coupon.schedule = coupon.schedule
    )[1]

    if (rate.type == 1) {
      floating.leg2 <- principal/(1 + as.numeric(rates[1]))^discount.time(
        tinitial = analysis.date,
        tfinal   = as.Date(Dates$effective.dates[1])
      )

      floating.leg <- (cashFlows.next * 1/(1 + as.numeric(rates[1]))^discount.time(
        tinitial = analysis.date,
        tfinal   = as.Date(Dates$effective.dates[1])
      )) + floating.leg2

    } else if (rate.type == 0) {
      floating.leg2 <- principal*exp(-(rates[1]*discount.time(
        tinitial = analysis.date,
        tfinal   = as.Date(Dates$effective.dates[1])
      )))
      floating.leg <- (cashFlows.next * exp(-(rates[1]*discount.time(
        tinitial = analysis.date,
        tfinal   = as.Date(Dates$effective.dates[1])
      )))) + floating.leg2
    }

    swapValue <- floating.leg - swapValue

    CouponSpread <- coupons(
      dates           = Dates$dates,
      coupon.rate     = spread,
      principal       = 1,
      freq            = freq,
      asset.type      = asset.type,
      daycount        = daycount,
      trade.date      = trade.date,
      coupon.schedule = coupon.schedule
    )

    # Calculate coupons of a bond with face value 0 and coupon payments
    # equal to the spread rate.
    CouponSpread[length(Dates$dates)] <- CouponSpread[length(Dates$dates)] - 1

    SpreadValue <- c()

      DiscountFactors <- discount.factors(
        rates         = rates,
        dates         = Dates$effective.dates,
        analysis.date = analysis.date,
        rate.type     = rate.type
      )

      for (i in seq_along(Dates$dates)) {
        SpreadValue[i] <- (CouponSpread[i] * as.numeric(DiscountFactors[i]))
      }

      SpreadValue <- sum(as.numeric(SpreadValue))

    # The swap value takes the previous swap value -which is equals to the par
    # value for the variable part, minus the fixed part value-, and adds a
    # spread bond with face value 0.
    swapValue   <- swapValue + SpreadValue

    # There are four possibilities of CCS, fixed and variable legs are valued the
    # same as IRS, using respective discount rates of each currency. The main diff
    # is that after valuation of variable leg in foreign currency, the equivalent
    # in value fixed leg is found by obtaining the TIR coupon rate. Once both legs
    # are in fixed terms, we use the basis rates to discount the fixed foreign leg.


  } else if (asset.type == "CCS") {
    if (Legs == "FF") {

      local.leg <- valuation.bonds(
        maturity        = maturity,
        coupon.rate     = coupon.rate,
        rates           = rates,
        analysis.date   = analysis.date,
        asset.type      = asset.type,
        freq            = freq,
        spread          = 0,
        rate.type       = rate.type,
        principal       = principal,
        trade.date      = trade.date,
        coupon.schedule = coupon.schedule
      )

      ex.leg <- valuation.bonds(
        maturity        = maturity,
        coupon.rate     = coupon.rate2,
        rates           = basis.rates,
        analysis.date   = analysis.date,
        asset.type      = asset.type,
        freq            = freq,
        spread          = 0,
        rate.type       = rate.type,
        principal       = principal2,
        trade.date      = trade.date,
        coupon.schedule = coupon.schedule
      )
    } else if (Legs == "VF") {
      #Value of variable local leg
      previous.date <- lubridate::`%m-%`(Dates$dates[1], months(12/freq))
      previous.date2 <- lubridate::`%m-%`(Dates$dates[2], months(24/freq))
      previous.date3 <- lubridate::`%m-%`(Dates$dates[3], months(36/freq))
      previous.date4 <- lubridate::`%m-%`(Dates$dates[4], months(48/freq))
      previous.date <- max(previous.date, previous.date2,
                           previous.date3, previous.date4,na.rm = TRUE)


      if (analysis.date %in% previous.date) {
        float.rate <- rates[1]
      }

      cashFlows.next <- coupons(
        dates           = Dates$dates,
        coupon.rate     = float.rate,
        principal       = principal,
        asset.type      = asset.type,
        daycount        = daycount,
        freq            = freq,
        trade.date      = trade.date,
        coupon.schedule = coupon.schedule
      )[1]

      DiscountFactors <- discount.factors(
        rates         = rates,
        dates         = Dates$effective.dates,
        analysis.date = analysis.date,
        rate.type     = rate.type
      )

      if (rate.type == 1) {
        floating.leg2 <- principal/(1 + as.numeric(rates[1]))^discount.time(
          tinitial = analysis.date,
          tfinal   = as.Date(Dates$effective.dates[1])
        )
      } else if (rate.type == 0) {
        floating.leg2<- principal*exp(-(rates[1]*discount.time(
          tinitial = analysis.date,
          tfinal   = as.Date(Dates$effective.dates[1])
        )))
      }

      floating.leg <- (cashFlows.next * DiscountFactors[1]) + floating.leg2


      local.leg <- floating.leg + valuation.bonds(
        maturity        = maturity,
        coupon.rate     = 0,
        rates           = rates,
        analysis.date   = analysis.date,
        asset.type      = asset.type,
        freq            = freq,
        spread          = spread,
        rate.type       = rate.type,
        principal       = principal,
        spread.only     = TRUE,
        trade.date      = trade.date,
        coupon.schedule = coupon.schedule
      )
      ex.leg <- valuation.bonds(
        maturity        = maturity,
        coupon.rate     = coupon.rate2,
        rates           = basis.rates,
        analysis.date   = analysis.date,
        asset.type      = asset.type,
        freq            = freq,
        spread          = 0,
        rate.type       = rate.type,
        principal       = principal2,
        trade.date      = trade.date,
        coupon.schedule = coupon.schedule
      )
    } else if (Legs == "FV" || Legs == "VV") {

      #Value of variable foreign leg
      previous.date <- lubridate::`%m-%`(Dates$dates[1], months(12/freq))
      previous.date2 <- lubridate::`%m-%`(Dates$dates[2], months(24/freq))
      previous.date3 <- lubridate::`%m-%`(Dates$dates[3], months(36/freq))
      previous.date4 <- lubridate::`%m-%`(Dates$dates[4], months(48/freq))
      previous.date <- max(previous.date, previous.date2,
                           previous.date3, previous.date4, na.rm = TRUE)


      if (analysis.date == previous.date) {
        float.rate2 <- rates2[1]
      }

      cashFlows.next <- coupons(
        dates           = Dates$dates,
        coupon.rate     = float.rate2,
        principal       = principal2,
        asset.type      = asset.type,
        daycount        = daycount,
        trade.date      = trade.date,
        freq            = freq,
        coupon.schedule = coupon.schedule
      )[1]

      DiscountFactors <- discount.factors(
        rates         = rates2,
        dates         = Dates$effective.dates,
        analysis.date = analysis.date,
        rate.type     = rate.type
      )

      if (rate.type == 1) {
        floating.leg2 <- principal2/(1 + as.numeric(rates2[1]))^discount.time(
          tinitial = analysis.date,
          tfinal   = as.Date(Dates$effective.dates[1])
        )
      } else if (rate.type == 0) {
        floating.leg2<- principal2*exp(-(rates2[1]*discount.time(
          tinitial = analysis.date,
          tfinal   = as.Date(Dates$effective.dates[1])
        )))
      }

      floating.leg <- (cashFlows.next * DiscountFactors[1]) + floating.leg2
      val.leg <- as.numeric(floating.leg) + valuation.bonds(
        maturity        = maturity,
        coupon.rate     = 0,
        rates           = rates2,
        analysis.date   = analysis.date,
        asset.type      = asset.type,
        freq            = freq,
        spread          = as.numeric(spread2),
        rate.type       = rate.type,
        principal       = as.numeric(principal2),
        spread.only     = TRUE,
        trade.date      = trade.date,
        coupon.schedule = coupon.schedule
      )

      Payments.ex <- valuation.bonds(
        maturity        = maturity,
        coupon.rate     = 0.04,
        rates           = rates2,
        analysis.date   = analysis.date,
        asset.type      = asset.type,
        freq            = freq,
        spread          = 0,
        rate.type       = rate.type,
        principal       = as.numeric(principal2),
        spread.only     = FALSE,
        trade.date      = trade.date,
        coupon.schedule = coupon.schedule
      )

      DiscountFactors <- discount.factors(
        rates         = rates2,
        dates         = Dates$effective.dates,
        analysis.date = analysis.date,
        rate.type     = rate.type,
        freq          = 1
      )

      #TIR represents the coupon rate that makes a fixed leg equivalent to variable foreign leg.

      TIR <- ((val.leg - (as.numeric(principal2) *
                            (DiscountFactors[length(DiscountFactors)]))) * 0.04) /
        (Payments.ex - (as.numeric(principal2) * (DiscountFactors[length(DiscountFactors)])) )

      coupon.rate2 <- TIR

      ex.leg <- valuation.bonds(
        maturity        = maturity,
        coupon.rate     = coupon.rate2,
        rates           = basis.rates,
        analysis.date   = analysis.date,
        asset.type      = asset.type,
        freq            = freq,
        spread          = 0,
        rate.type       = rate.type,
        principal       = principal2,
        trade.date      = trade.date,
        coupon.schedule = coupon.schedule
      )

        if (Legs == "FV") {
          local.leg <- valuation.bonds(
            maturity        = maturity,
            coupon.rate     = coupon.rate,
            rates           = rates,
            analysis.date   = analysis.date,
            asset.type      = asset.type,
            freq            = freq,
            spread          = 0,
            rate.type       = rate.type,
            principal       = principal,
            trade.date      = trade.date,
            coupon.schedule = coupon.schedule
          )
        } else if (Legs == "VV") {
          local.leg <- principal + valuation.bonds(
            maturity        = maturity,
            coupon.rate     = 0,
            rates           = rates,
            analysis.date   = analysis.date,
            asset.type      = asset.type,
            freq            = freq,
            spread          = spread,
            rate.type       = rate.type,
            principal       = principal,
            spread.only     = TRUE,
            trade.date      = trade.date,
            coupon.schedule = coupon.schedule
          )
        }
    }

swapValue <- local.leg - (ex.leg * ex.rate)
}
  return(swapValue)
}

#' From one price to another
#'
#' @description Converts bond prices from dirty to clean and viceversa.
#'
#' @inheritParams coupon.dates
#' @inheritParams coupons
#' @inheritParams discount.factors
#' @param dirty Numeric value. Determines if the input price corresponds to the
#' dirty price or the clean price. For dirty price, set \code{dirty = 1}.
#' Otherwise, \code{dirty = 0}
#' @param price Numeric value. Price of the bond to convert.
#'
#' @details
#' \code{asset.type} makes reference to the following type of assets:
#' \itemize{
#'    \item "TES" for Colombian Treasury Bonds (default).
#'    \item "FixedIncome" for assets that are indexed to a fixed income with
#'    different frequency of payments.
#'    \item "IBR" for bonds and assets indexed to 3M IBR rate.
#'    \item "LIBOR" for bonds and assets indexed to 3M LIBOR.
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
#' @return The dirty price or clean price of a bond.
#' @export
#'
#' @examples
#' price.dirty2clean(maturity = "2026-01-03", analysis.date = "2023-01-02",
#'                   price = 1, dirty = 1, coupon.rate = 0.04, principal = 1)
#' price.dirty2clean(maturity = "2026-01-03", analysis.date = "2023-01-02",
#'                   price = 0.9601096, dirty = 0, coupon.rate = 0.04, principal = 1)
#
price.dirty2clean <- function(maturity, analysis.date = Sys.Date(), price, dirty = 1,
                              coupon.rate, principal = 1, asset.type = "TES", freq = NULL,
                              daycount = "NL/365") {
  # Calculate the accrued interest from the last coupon payment until date
  # of analysis. Depending on the input price of the asset, returns the dirty
  # or clean price. If the input is a dirty price, it subtracts accrued
  # interest to the input price.
  if (analysis.date >= maturity) {
    stop("analysis date is after maturity")
  }

  return(price + ifelse(dirty == 1, -1, 1) * accrued.interests(
    analysis.date = analysis.date,
    coupon.rate   = coupon.rate,
    maturity      = maturity,
    asset.type    = asset.type,
    daycount      = daycount,
    freq          = freq,
    principal     = principal
  ))
}

#' From a price to a rate
#'
#' @description Calculates the Internal Rate of Return (IRR) of a given asset
#' taking into account the market price, maturity, face value, and analysis
#' date.
#'
#' @inheritParams coupon.dates
#' @inheritParams coupons
#' @inheritParams discount.factors
#' @inheritParams valuation.bonds
#' @inheritParams price.dirty2clean
#'
#' @details
#' \code{asset.type} makes reference to the following type of assets:
#' \itemize{
#'    \item "TES" for Colombian Treasury Bonds (default).
#'    \item "FixedIncome" for assets that are indexed to a fixed income with
#'    different frequency of payments.
#'    \item "IBR" for bonds and assets indexed to 3M IBR rate.
#'    \item "LIBOR" for bonds and assets indexed to 3M LIBOR.
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
#' @details
#' \code{convention} makes reference to the following type of business day conventions:
#' \itemize{
#'    \item "F" for Following business day convention.
#'    \item "MF" for Modified Following business day convention.
#'    \item "B" for Backward business day convention.
#'    \item "MB" for Modified Backward business day convention.
#'      }
#'
#' @details
#' \code{coupon.schedule} makes reference to the following type of coupon payment schedule
#' of a bond:
#' \itemize{
#'    \item "LF" for Long First coupon payment.
#'    \item "LL" for Long Last coupon payment.
#'    \item "SF" for Short First coupon payment.
#'    \item "SL" for Short Last coupon payment.
#'      }
#'
#' @return The Yield to Maturity or Internal Rate of Return of a given asset.
#'
#' @importFrom stats uniroot
#' @export
#'
#' @examples
#' bond.price2rate(maturity = "2023-01-03", analysis.date = "2021-01-03",
#'                 price = 1, coupon.rate = 0.04, principal = 1,
#'                 asset.type = "TES", freq = 1)
#'
bond.price2rate <- function(maturity, analysis.date = Sys.Date(), price,
                            coupon.rate, principal = 1, asset.type = "TES",
                            freq = NULL, rate.type = 1, spread = 0, dirty = 1,
                            daycount = "NL/365", convention = "F",
                            trade.date = NULL, coupon.schedule = "SF") {
  fun_dif <- function(rate) {
    valuation.bonds(
      coupon.rate     = coupon.rate,
      analysis.date   = analysis.date,
      maturity        = maturity,
      rates           = rate,
      principal       = principal,
      dirty           = dirty,
      rate.type       = rate.type,
      daycount        = daycount,
      freq            = freq,
      spread          = spread,
      asset.type      = asset.type,
      convention      = convention,
      coupon.schedule = coupon.schedule,
      trade.date      = trade.date
    ) - price
    }

  # Given an asset price, the function calculates the Internal Rate of Return
  # of the asset that makes the estimated price with the valuation formula
  # equal to the input/market price.
  return(uniroot(fun_dif, c(-0.99, 1000000), tol = 1e-07, maxiter = 200)$root)
}

#' Bond Sensitivity
#'
#' @description Calculates the sensitivity of a given bond by numerically averaging
#' the percentage change in bonds price when moving upwards and downwards, by 1 basic point,
#' the Yield to Maturity vector.
#'
#' @inheritParams coupon.dates
#' @inheritParams coupons
#' @inheritParams discount.factors
#' @inheritParams valuation.bonds
#' @inheritParams price.dirty2clean
#' @param input String that establishes if the \code{price} input corresponds to the
#' Internal Rate of Return (IRR) of the bond or the market price. Set
#' \code{"rate"} for the IRR. Otherwise, \code{"price"}.
#' @param price Numeric value of either market price or Internal Rate of Return of a
#' given bond. Instead of IRR, can also be a vector of multiple rates, one for every
#' coupon date.
#'
#' @details
#' \code{asset.type} makes reference to the following type of assets:
#' \itemize{
#'    \item "TES" for Colombian Treasury Bonds (default).
#'    \item "FixedIncome" for assets that are indexed to a fixed income with
#'    different frequency of payments.
#'    \item "IBR" for bonds and assets indexed to 3M IBR rate.
#'    \item "LIBOR" for bonds and assets indexed to 3M LIBOR.
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
#' @details
#' \code{convention} makes reference to the following type of business day conventions:
#' \itemize{
#'    \item "F" for Following business day convention.
#'    \item "MF" for Modified Following business day convention.
#'    \item "B" for Backward business day convention.
#'    \item "MB" for Modified Backward business day convention.
#'      }
#'
#' @details
#' \code{coupon.schedule} makes reference to the following type of coupon payment schedule
#' of a bond:
#' \itemize{
#'    \item "LF" for Long First coupon payment.
#'    \item "LL" for Long Last coupon payment.
#'    \item "SF" for Short First coupon payment.
#'    \item "SL" for Short Last coupon payment.
#'      }
#'
#' @return Bond sensitivity
#' @export
#'
#' @examples
#' sens.bonds(input = c("price"), price = 0.98, maturity = "2023-01-03",
#'            analysis.date = "2019-01-05", coupon.rate = 0.04,
#'            principal = 1, asset.type = "IBR", rate.type = 1)
#' sens.bonds(input = c("rate"), price = rep(0.08,8), maturity = "2023-01-03",
#'            analysis.date = "2015-02-03", coupon.rate = 0.04,
#'            principal = 1, asset.type = "FixedIncome", freq = 1,
#'            rate.type = 1, daycount = "ACT/365", dirty = 1,
#'            convention = "MB", trade.date = "2015-02-03",
#'            coupon.schedule = "LF")
#'
sens.bonds <- function(input, price, maturity, analysis.date = Sys.Date(),
                       coupon.rate, principal = 1, asset.type = "TES",
                       freq = 1, rate.type = 1,  spread = 0,
                       daycount = "ACT/365", dirty = 1, convention = "F",
                       trade.date = NULL, coupon.schedule = "SF") {

  # Transforms the input of price to rate if given input is the market price.
  rate <- ifelse(input == "price", bond.price2rate(
    price           = price,
    analysis.date   = analysis.date,
    coupon.rate     = coupon.rate,
    maturity        = maturity,
    dirty           = dirty,
    spread          = spread,
    rate.type       = rate.type,
    freq            = freq,
    daycount        = daycount,
    principal       = principal,
    asset.type      = asset.type,
    convention      = convention,
    coupon.schedule = coupon.schedule,
    trade.date      = trade.date
  ), price)

  if (length(price) == 1) {
  # Calculate the value of a bond with a positive movement in the Yield to
  # Maturity of one percentage point.
  priceF1 <- valuation.bonds(
    rates           = rate + 0.01,
    maturity        = maturity,
    principal       = principal,
    coupon.rate     = coupon.rate,
    rate.type       = rate.type,
    asset.type      = asset.type,
    daycount        = daycount,
    analysis.date   = analysis.date,
    freq            = freq,
    spread          = spread,
    dirty           = dirty,
    convention      = convention,
    coupon.schedule = coupon.schedule,
    trade.date      = trade.date

  )

  # Calculate the value of a bond with a negative movement in the Yield to
  # Maturity of one percentage point.
  priceF2 <- valuation.bonds(
    rates           = rate - 0.01,
    maturity        = maturity,
    principal       = principal,
    coupon.rate     = coupon.rate,
    rate.type       = rate.type,
    asset.type      = asset.type,
    daycount        = daycount,
    analysis.date   = analysis.date,
    freq            = freq,
    spread          = spread,
    dirty           = dirty,
    convention      = convention,
    coupon.schedule = coupon.schedule,
    trade.date      = trade.date
  )

  # Calculate the value of a bond without any movement in the Yield to
  # Maturity.
  priceF <- ifelse (input == "price", price, (valuation.bonds(
    rates           = rate,
    maturity        = maturity,
    principal       = principal,
    coupon.rate     = coupon.rate,
    rate.type       = rate.type,
    asset.type      = asset.type,
    daycount        = daycount,
    analysis.date   = analysis.date,
    freq            = freq,
    spread          = spread,
    dirty           = dirty,
    convention      = convention,
    coupon.schedule = coupon.schedule,
    trade.date      = trade.date
  )))
  } else {

    priceF1 <- valuation.bonds(
      rates           = price + 0.01,
      maturity        = maturity,
      principal       = principal,
      coupon.rate     = coupon.rate,
      rate.type       = rate.type,
      asset.type      = asset.type,
      daycount        = daycount,
      analysis.date   = analysis.date,
      freq            = freq,
      spread          = spread,
      dirty           = dirty,
      convention      = convention,
      coupon.schedule = coupon.schedule,
      trade.date      = trade.date
    )

    # Calculate the value of a bond with a negative movement in the Yield to
    # Maturity of one percentage point.
    priceF2 <- valuation.bonds(
      rates           = price - 0.01,
      maturity        = maturity,
      principal       = principal,
      coupon.rate     = coupon.rate,
      rate.type       = rate.type,
      asset.type      = asset.type,
      daycount        = daycount,
      analysis.date   = analysis.date,
      freq            = freq,
      spread          = spread,
      dirty           = dirty,
      convention      = convention,
      coupon.schedule = coupon.schedule,
      trade.date      = trade.date
    )

    # Calculate the value of a bond without any movement in the Yield to
    # Maturity.
    priceF <- ifelse (input == "price", price, (valuation.bonds(
      rates           = price,
      maturity        = maturity,
      principal       = principal,
      coupon.rate     = coupon.rate,
      rate.type       = rate.type,
      asset.type      = asset.type,
      daycount        = daycount,
      analysis.date   = analysis.date,
      freq            = freq,
      spread          = spread,
      dirty           = dirty,
      convention      = convention,
      coupon.schedule = coupon.schedule,
      trade.date      = trade.date
    )))
}
  # Returns the percentage sensitivity between an upward movement and a
  # downward movement of the Yield to Maturity. Divides by 0.01 to give output
  # in percentage.

  return((priceF2-priceF1) / (2*priceF*0.01))
}



#' Weighted Average Life
#'
#' @description Calculates the weighted average life of a given bond by dividing the weighted
#' total payments by the total payments.
#'
#' @inheritParams coupon.dates
#' @inheritParams coupons
#' @inheritParams discount.factors
#' @inheritParams valuation.bonds
#' @inheritParams price.dirty2clean
#' @param input String that establishes if the price input corresponds to the
#' Internal Rate of Return (IRR) of the bond or the market price. Set
#' \code{"rate"} for the IRR. Otherwise, \code{"price"}.
#' @param price Numeric value of either market price or Internal Rate of Return (IRR) of a
#' given bond. Instead of IRR, can also be a rates vector that corresponds to coupon dates.
#'
#' @details
#' \code{asset.type} makes reference to the following type of assets:
#' \itemize{
#'    \item "TES" for Colombian Treasury Bonds (default).
#'    \item "FixedIncome" for assets that are indexed to a fixed income with
#'    different frequency of payments.
#'    \item "IBR" for bonds and assets indexed to 3M IBR rate.
#'    \item "LIBOR" for bonds and assets indexed to 3M LIBOR.
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
#' @details
#' \code{convention} makes reference to the following type of business day conventions:
#' \itemize{
#'    \item "F" for Following business day convention.
#'    \item "MF" for Modified Following business day convention.
#'    \item "B" for Backward business day convention.
#'    \item "MB" for Modified Backward business day convention.
#'      }
#'
#' @details
#' \code{coupon.schedule} makes reference to the following type of coupon payment schedule
#' of a bond:
#' \itemize{
#'    \item "LF" for Long First coupon payment.
#'    \item "LL" for Long Last coupon payment.
#'    \item "SF" for Short First coupon payment.
#'    \item "SL" for Short Last coupon payment.
#'      }
#'
#' @return Weighted average life of given bond
#' @export
#'
#' @examples
#' average.life(input = c("rate"), price = 0.08, maturity = "2026-06-01",
#'              analysis.date = "2025-06-01", coupon.rate = 0.06, principal = 1000,
#'              asset.type = "IBR", freq = 4)
#' average.life(input = c("rate"), price = c(0.043,0.05), maturity = "2023-01-03",
#'              analysis.date = "2021-01-03", coupon.rate = 0.04, principal = 1,
#'              asset.type = "FixedIncome", freq = 1, rate.type = 0)
#'

average.life <- function(input, price, maturity, analysis.date = Sys.Date(),
                         coupon.rate, principal = 1, asset.type = "TES",
                         freq = 1, rate.type = 1,  spread = 0,
                         daycount = "ACT/365", dirty = 1, convention = "F",
                         trade.date = NULL, coupon.schedule = "SF") {


  # Transforms the input of price to rate if given input is the market price.
  if (input == "price") {
  price <- bond.price2rate(
    price           = price,
    analysis.date   = analysis.date,
    coupon.rate     = coupon.rate,
    maturity        = maturity,
    dirty           = dirty,
    spread          = spread,
    rate.type       = rate.type,
    freq            = freq,
    daycount        = daycount,
    principal       = principal,
    asset.type      = asset.type,
    convention      = convention,
    coupon.schedule = coupon.schedule,
    trade.date      = trade.date
  )
  } else if (input == "rate")  {
    price <- price
  }

# Calculate coupon dates of an asset.
Dates <- coupon.dates(
  maturity        = maturity,
  asset.type      = asset.type,
  analysis.date   = analysis.date,
  freq            = freq,
  convention      = convention,
  loc             = "BOG",
  coupon.schedule = coupon.schedule,
  trade.date      = trade.date
)

cashFlows <- coupons(
  dates           = Dates$dates,
  coupon.rate     = coupon.rate,
  principal       = principal,
  asset.type      = asset.type,
  daycount        = daycount,
  freq            = freq,
  coupon.schedule = coupon.schedule,
  trade.date      = trade.date
)

Discount.Factors <- c()
if (length(price) > 1) {
  Discount.Factors <- discount.factors(
  rates         = price,
  dates         = Dates$effective.dates,
  analysis.date = analysis.date,
  rate.type     = rate.type
  )
  } else {
  if (rate.type == 1) {
  for (i in seq_along(Dates$dates)) {
    Discount.Factors[i] <- 1 / (1 + as.numeric(price))^discount.time(
    tinitial = analysis.date,
    tfinal   = Dates$effective.dates[i]
  )
  }
  } else if (rate.type == 0) {

    for (i in seq_along(Dates$dates)) {
      Discount.Factors[i] <- exp((-(as.numeric(price) * discount.time(
        tinitial = analysis.date,
        tfinal   = Dates$effective.dates[i]
      ))))
    }
  }
}

times <- c()
for (i in seq_along(Dates$dates)) {
  times[i] <- quantdates::day_count(
  tinitial   = analysis.date,
  tfinal     = Dates$dates[i],
  convention = daycount
  )
  }

if (length(price) > 1) {
priceF <- valuation.bonds(
  rates           = price,
  maturity        = maturity,
  principal       = principal,
  coupon.rate     = coupon.rate,
  rate.type       = rate.type,
  asset.type      = asset.type,
  daycount        = daycount,
  analysis.date   = analysis.date,
  freq            = freq,
  spread          = spread,
  dirty           = dirty,
  convention      = convention,
  coupon.schedule = coupon.schedule,
  trade.date      = trade.date
)
} else {
  priceF <- valuation.bonds(
    rates           = price,
    maturity        = maturity,
    principal       = principal,
    coupon.rate     = coupon.rate,
    rate.type       = rate.type,
    asset.type      = asset.type,
    daycount        = daycount,
    analysis.date   = analysis.date,
    freq            = freq,
    spread          = spread,
    dirty           = dirty,
    convention      = convention,
    coupon.schedule = coupon.schedule,
    trade.date      = trade.date
  )
}

return ((sum(cashFlows*Discount.Factors*times))/(priceF))
}

#' Basis Curve
#'
#' @description Function that calibrates a "discount basis rate" curve according
#' to data of cross currency swaps. Available methods are bootstrapping or residual
#' sum of squares (RSS) between the value of a given or inferred fixed foreign leg
#' and the value of the local leg.
#'
#' @inheritParams coupon.dates
#' @inheritParams coupons
#' @inheritParams discount.factors
#' @inheritParams valuation.bonds
#' @inheritParams valuation.swaps
#' @inheritParams curve.calculation
#'
#' @param swaps Matrix containing relevant information of cross currency swaps where each
#' row represents a swap and each column represents the next attributes: maturity, legs,
#' coupon rate of local leg, coupon rate of foreign leg, spread of local leg, spread of
#' variable leg, principal of local and principal of variable leg. \code{colnames(swaps)}
#' must be defined as the following: \code{c("Mat", "Legs", "C1", "C2", "spread1", "spread2", "prin1", "prin2")}.
#' @param rates Discount rates given by the local zero coupon rate curve. The curve has to
#' have nodes with at least, with 3 decimals.
#' @param rates2 Discount rates given by the foreign zero coupon rate curve. The curve has to
#' have nodes with at least, with 3 decimals.
#' @param npieces Number of constant or linear segments for the curve to have. By
#' default \code{NULL}, and bootstrapping method is used, otherwise, minimization of
#' RSS is used.
#' @param Weights Vector of weights used to dot product with residual squares in order to
#' calculate residual sum of squares. By default, each residual is assigned the same weight.
#' @param nsimul Number of simulations for the terms of the pieces. The more simulations,
#' the more likely to find a better local solution. By default \code{1}, and terms are
#' defined in such way each piece occupies the same length in the abscissa axis.
#' @param piece.term Vector that establishes a unique term structure for optimization to take place.
#' Each piece or segment must have a unique maturity, as numeric value in years,
#' that signifies the end of the segment. Last segment maturity must not be introduced, it is assumed to be equivalent to
#' the last term introduced on analysis date. Therefore, the \code{piece.term} vector must always have a length equal to \code{npieces} - 1.
#' @param obj String related to the definition of the error in the RSS methodology.
#' Set \code{"Price"} for minimization of error by price or \code{"Rate"} for
#' minimization of error by rate.
#' @param approximation String that establish the approximation. Set \code{'linear'} for a
#' linear approximation, or \code{'constant'} for a piecewise constant function.
#' @param nodes Desired output nodes of the curve.
#'
#'
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
#'     }
#'
#' @details
#' \code{swaps["Legs"]} makes reference to the following types of legs composition of the
#' cross currency swaps.
#' \itemize{
#'    \item "FF" for fixed leg in local currency and fixed leg in foreign currency.
#'    \item "FV" for fixed leg in local currency and variable leg in foreign currency.
#'    \item "VF" for variable leg in local currency and fixed leg in foreign currency.
#'    \item "VV" for variable leg in local currency and variable leg in foreign currency.
#'      }
#'
#' @author Camilo Díaz
#' @return Constant or Linear piecewise basis curve.
#' @importFrom stats approx optim stepfun runif
#' @export
#'
#' @examples
#' # Inputs for calibration of spot curve
#' yield.curve  <- c(0.015,0.0175, 0.0225, 0.0275, 0.0325, 0.0375,0.04,0.0425,0.045,0.0475,0.05)
#' names(yield.curve) <- c(0.5,1,2,3,4,5,6,7,8,9,10)
#' nodes <- seq(0,10,0.001)
#' # Calibration of local spot curve
#' rates <- curve.calibration (yield.curve = yield.curve, market.assets = NULL,
#'                             analysis.date = "2019-01-03" , asset.type = "IBRSwaps",
#'                             freq = 4, rate.type = 0, fwd = 0, npieces = NULL,
#'                             obj = "Price", nodes = nodes, approximation = "linear")
#' # Input for Basis Curve
#' ex.rate <- 4814
#' swaps <- rbind(c("2024-03-01", "FF", 0.07 , 0.0325, NA   , NA    , 2000 * ex.rate, 2000),
#'                c("2025-03-01", "VV", NA   , NA    , 0.015, 0.0175, 2000 * ex.rate, 2000),
#'                c("2026-03-01", "FF", 0.075, 0.03  , NA   ,  NA   , 5000000, 5000000 / ex.rate),
#'                c("2027-03-01", "VV", NA   , NA    , 0.01 , 0.015 , 5000000, 5000000 / ex.rate),
#'                c("2028-03-01", "FF", 0.08 ,0.035  , NA   , NA    , 3000000, 3000000 / ex.rate),
#'                c("2029-03-01", "VV", NA   , NA    , 0.01 , 0.0125, 3000000, 3000000 / ex.rate))
#' colnames(swaps) <- c("Mat"  ,"Legs", "C1" , "C2", "spread1", "spread2", "prin1", "prin2")
#' # Function
#' basis.curve(swaps = swaps, ex.rate = 4814, analysis.date = "2023-03-01",
#'             rates = rates, rates2 = rates / 4, freq = c(2,2,2,2,1,1),
#'             rate.type = 1, npieces = 4, obj = "Price", Weights = NULL,
#'             nsimul = 1, nodes = nodes, approximation = "linear")


basis.curve <- function(swaps, ex.rate = NULL, analysis.date = Sys.Date(),
                        rates, rates2, freq = 1, rate.type = 1,
                        daycount = "ACT/365", npieces = NULL, obj = "Price",
                        Weights = NULL, nsimul = 1, piece.term = NULL,
                        nodes = seq(0,15,0.001), approximation = "constant") {

  asset.type = "CCS"
  colnames(swaps) <- c("Mat", "Legs", "C1", "C2", "spread1", "spread2", "prin1", "prin2")

  if (length(freq) == 1) {
    freq <- rep(freq, NROW(swaps))
  } else if (!length(freq) == 1) {
    if (!length(freq) == NROW(swaps)) {
      stop(paste0(deparse(sys.call()), ":","
                if freq as vector, then length has to be", ": ", NROW(swaps)),
           call. = FALSE)
    }
  }

  swaps.ccs <- seq(1,NROW(swaps),1)

  PaymentDates <- lapply(swaps.ccs, function(p) {
    Dates <- coupon.dates(
      maturity      = swaps[p,"Mat"] ,
      asset.type    = asset.type,
      analysis.date = analysis.date,
      freq          = freq[p]
    )$dates
  })

  Payments <- c()
  Payments <- lapply(seq_along(swaps.ccs), function(p) {

    for (i in 1:(length(PaymentDates[[p]]))) {
      Payments <- c(Payments, discount.time(
        tinitial   = analysis.date,
        tfinal     = PaymentDates[[p]][i]
      ))
    }
  Payments <- round(Payments,3)
    return(Payments)
  })


  rates <- lapply(seq_along(swaps.ccs), function(p) {
    termpayment <- c()
    for (i in 1:(length(Payments[[p]]))) {
      termpayment <- c(termpayment, which(names(rates) == Payments[[p]][i]))
    }

    rates <- rates[termpayment]
    return(rates)
  })

  rates2 <- lapply(seq_along(swaps.ccs), function(p) {
    termpayment <- c()
    for (i in 1:(length(Payments[[p]]))) {
      termpayment <- c(termpayment, which(names(rates2) == Payments[[p]][i]))
    }

    rates2 <- rates2[termpayment]
    return(rates2)
  })

  # Local leg is valued normally with the discount local rates. Variable legs are valued par + spread.

  local.leg <- lapply(swaps.ccs, function(p) {
    j <- which(swaps.ccs == p)
    swap.i <- swaps[p, ]
    mat    <- as.Date(swaps[p,"Mat"])

    if (swap.i["Legs"] == "FF") {
      local.leg <- valuation.bonds(
        maturity      = mat,
        coupon.rate   = as.numeric(swap.i["C1"]),
        rates         = rates[[j]],
        analysis.date = analysis.date,
        asset.type    = asset.type,
        freq          = freq[j],
        spread        = 0,
        rate.type     = rate.type,
        principal     = as.numeric(swap.i["prin1"])
      )
    } else if (swap.i["Legs"] == "VF") {
      local.leg <- as.numeric(swap.i["prin1"]) + valuation.bonds(
        maturity = mat,
        coupon.rate = 0,
        rates = rates[[j]],
        analysis.date = analysis.date,
        asset.type = asset.type,
        freq = freq[j],
        spread = as.numeric(swap.i["spread1"]),
        rate.type = rate.type,
        principal = as.numeric(swap.i["prin1"]),
        spread.only = TRUE
      )
    } else if (swap.i["Legs"] == "FV") {
      local.leg <- valuation.bonds(
        maturity = mat,
        coupon.rate = as.numeric(swap.i["C1"]),
        rates = rates[[j]],
        analysis.date = analysis.date,
        asset.type = asset.type,
        freq = freq[j],
        spread = 0,
        rate.type = rate.type,
        principal = as.numeric(swap.i["prin1"])
      )

    } else if (swap.i["Legs"] == "VV") {
      local.leg <- as.numeric(swap.i["prin1"]) + valuation.bonds(
        maturity = mat,
        coupon.rate = 0,
        rates = rates[[j]],
        analysis.date = analysis.date,
        asset.type = asset.type,
        freq = freq[j],
        spread = as.numeric(swap.i["spread1"]),
        rate.type = rate.type,
        principal = as.numeric(swap.i["prin1"]),
        spread.only = TRUE
      )
    }
  })

  for (p in swaps.ccs) {
    j <- which(swaps.ccs == p)
    swap.i <- swaps[p, ]
    mat    <- swaps[p,"Mat"]

    # For variable foreign legs, the TIR coupon rate that makes a fixed leg equivalent
    # in value to a variable leg is found. Therefore, all foreign legs are transformed,
    # if necessary,  to fixed legs.

    if (swap.i["Legs"] == "FV" || swap.i["Legs"] == "VV") {
      val.leg <- as.numeric(swap.i["prin2"]) + valuation.bonds(
        maturity      = mat,
        coupon.rate   = 0,
        rates         = rates2[[j]],
        analysis.date = analysis.date,
        asset.type    = asset.type,
        freq          = freq[j],
        spread        = as.numeric(swap.i["spread2"]),
        rate.type     = rate.type,
        principal     = as.numeric(swap.i["prin2"]),
        spread.only   = TRUE
      )

      Payments.ex <- valuation.bonds(
        maturity = mat,
        coupon.rate = 0.04,
        rates = rates2[[j]],
        analysis.date = analysis.date,
        asset.type = asset.type,
        freq = freq[j],
        spread = 0,
        rate.type = rate.type,
        principal = as.numeric(swap.i["prin2"]),
        spread.only = FALSE
      )

      Dates <- coupon.dates(
        maturity      = mat,
        asset.type    = asset.type,
        analysis.date = analysis.date,
        freq          = freq[j],
        loc           = "BOG"
      )

      DiscountFactors <- discount.factors(
        rates         = rates2[[j]],
        dates         = Dates$effective.dates,
        analysis.date = analysis.date,
        rate.type     = rate.type
      )



      TIR <- ((val.leg - (as.numeric(swap.i["prin2"]) *
                            (DiscountFactors[length(DiscountFactors)]))) * 0.04) /
        (Payments.ex - (as.numeric(swap.i["prin2"]) * (DiscountFactors[length(DiscountFactors)])) )

      swaps[p,"C2"] <- TIR

    }
  }

  # After we have all foreign legs as fixed legs, coupons are calculated and the discount
  # factors that make those sum(coupons*df) equivalent to local leg value are found. This is done
  # either by bootstrapping or by minimizing an error^2. Note that in this case Q0 * exp(-rt),
  # where Q0 is inital exchange rate, allows us to bring to net present value in COP,
  # future cashflows in dollars. Here we find optimal r of Q0 * exp(-rt) for any t.
  # Therefore, with this "basis discount rate", we can value in COP, any fixed foreign leg.

  cashFlows <- lapply(swaps.ccs, function(p) {
    j <- which(swaps.ccs == p)
    swap.i <- swaps[p, ]

    cashFlows <- coupons(
      dates       = PaymentDates[[j]],
      coupon.rate = as.numeric(swap.i["C2"]),
      asset.type  = asset.type,
      daycount    = daycount,
      freq        = freq[j],
      principal   = as.numeric(swap.i["prin2"])
    )
  })

      # Function that determines the approximation error between rates and
      # market bond prices. Builds, recursively, a curve of rates where
      # the current curve has the available zero coupon rates, and the final
      # node are the rates that are calibrated through minimization of
      # the error.
      RateFunction1 <- function(Payments_i, Dates, Today, VP,
                               CurrentCurve, FinalNode) {

        CurrentCurve[1] <- FinalNode
        Curva        <- c(CurrentCurve, FinalNode)
        names(Curva) <- c(names(CurrentCurve), (Dates[length(Dates)] - Today) / 365)

        PaymentRates <- approx(
          x      = round(as.numeric(names(Curva)), 6),
          y      = Curva,
          xout   = round(as.numeric((Dates - Today) / 365), 6),
          method = "linear",
          ties   = min
        )$y

        # Use an approx. to bootstrap the rates between the available
        # calibrated rates and the input rate (Final Node).
        return(abs(sum(Payments_i / (1 + PaymentRates)^as.numeric((Dates - Today) / 365)) * ex.rate - VP))
      }

      RateFunction2 <- function(Payments_i, Dates, Today, VP,
                               CurrentCurve, FinalNode) {
        CurrentCurve[1] <- FinalNode
        Curva        <- c(CurrentCurve, FinalNode)
        names(Curva) <- c(names(CurrentCurve), as.numeric((Dates[length(Dates)] - Today) / 365))
        PaymentRates <- stepfun(
          x     = as.numeric(names(Curva))[-1],
          y     = Curva,
          right = FALSE,
          f     = 1,
          ties   = min
        )

        # Use an StepFun to bootstrap the rates between the available
        # calibrated rates and the input rate (Final Node). With stepfun
        # we are building a constant piecewise curve.
        return(abs(sum(Payments_i / (1 + PaymentRates(as.numeric((Dates - Today) / 365)))^as.numeric((Dates - Today) / 365))
                   * ex.rate - VP))
      }

      RateFunction3 <- function(Payments_i, Dates, Today, VP,
                               CurrentCurve, FinalNode) {
        CurrentCurve[1] <- FinalNode
        Curva        <- c(CurrentCurve, FinalNode)
        names(Curva) <- c(names(CurrentCurve), (Dates[length(Dates)] - Today) / 365)

        PaymentRates <- approx(
          x      = round(as.numeric(names(Curva)), 6),
          y = Curva,
          xout   = round(as.numeric((Dates - Today) / 365), 6),
          method = "linear",
          ties   = min
        )$y

        return(abs(sum(Payments_i * exp((-(PaymentRates) * as.numeric((Dates - Today) / 365))))* ex.rate - VP))
      }

      RateFunction4 <- function(Payments_i, Dates, Today, VP,
                               CurrentCurve, FinalNode) {
        CurrentCurve[1] <- FinalNode
        Curva        <- c(CurrentCurve, FinalNode)
        names(Curva) <- c(names(CurrentCurve), (Dates[length(Dates)] - Today) / 365)
        PaymentRates <- stepfun(
          x     = as.numeric(names(Curva))[-1],
          y     = Curva,
          right = FALSE,
          f     = 1,
          ties   = min
        )

        return(abs(sum(Payments_i * exp((-(PaymentRates(as.numeric((Dates - Today) / 365))) *
                                           as.numeric((Dates - Today) / 365))))* ex.rate - VP))
      }

if (rate.type == 1) {
  if (approximation == "linear") {
    RateFunction <- RateFunction1
  } else if (approximation == "constant") {
    RateFunction <- RateFunction2
  }
} else if (rate.type == 0) {
  if (approximation == "linear") {
    RateFunction <- RateFunction3
  } else if (approximation == "constant") {
    RateFunction <- RateFunction4
  }
}

  CurrentCurve        <- as.numeric(0)
  names(CurrentCurve) <- "0"

  # Function that minimizes the absolute error between a calculated bond with
  # the input rate and the market value. The output FinalNode reflects a zero
  # coupon rate equivalent that is calibrated with the market conditions.
  for (i in seq_along(swaps.ccs)) {
    Dates      <- PaymentDates[[i]]
    Payments_i <- cashFlows[[i]]
    Today      <- as.Date(analysis.date)
    VP         <- local.leg[[i]]

    FinalNode <- optim(
      par          = as.numeric(rates[[i]][1]),
      fn           = RateFunction,
      method       = "Brent",
      lower        = -0.5,
      upper        = 1,
      Payments_i   = Payments_i,
      Dates        = Dates,
      Today        = Today,
      VP           = VP,
      CurrentCurve = CurrentCurve
    )

    CurrentCurve        <- c(CurrentCurve, FinalNode$par)
    names(CurrentCurve) <- c(names(CurrentCurve)[-length(CurrentCurve)],
                             (Dates[length(Dates)] - Today) / 365)
  }

  if (! is.null(npieces)) {
    RateFunction.for <- function(FinalNode,PaymentDates, Payments,
                             Today, local.leg) {
      # Constructs a time payment structure for calculating the integral of the
      # forward rates and time difference matrix for the integrals.
      Curva        <- as.numeric(FinalNode [1:npieces])
      if (npieces == 1) {
        names(Curva)[length(Curva)] <- end.date
      } else {
        if (! is.null(piece.term)) {
          names(Curva) <- as.numeric(piece.term)
        } else if (is.null(piece.term)) {
          names(Curva) <- as.numeric(FinalNode[(npieces + 1):(2 * npieces - 1)])
        }
        names(Curva)[length(Curva)] <- end.date
      }

      CurrentCurve        <- as.numeric(Curva[1])
      names(CurrentCurve) <- "0"

      if (! names(CurrentCurve[1]) == 0) {
        CurrentCurve <- c(CurrentCurve[1],CurrentCurve)
        names(CurrentCurve)[1] <- "0"
      }

      Curva <- c(CurrentCurve,Curva)

      df <- lapply(1:(length(PaymentDates)), function (p) {
        termpayment <- as.numeric(PaymentDates[[p]] - Today) / 365

        if (approximation == "constant") {
          Curva[1] <- Curva[2]
          PaymentRates <- stepfun(
            x     = as.numeric(names(Curva))[-1],
            y     = Curva,
            right = FALSE,
            f     = 1,
            ties   = min
          )

          if (rate.type == 1) {
            df <- (1 / ((1 + PaymentRates(termpayment)) ^ (termpayment))) * ex.rate
          } else if (rate.type == 0) {
            df <- exp(- PaymentRates(termpayment) * termpayment) * ex.rate
          }

        } else if  (approximation == "linear") {
          PaymentRates <- approx(
            x      = round(as.numeric(names(Curva)), 6),
            y      = Curva,
            xout   = round(termpayment, 6),
            method = "linear",
            ties   = min
          )$y

          if (rate.type == 1) {
            df <- (1 / ((1 + PaymentRates) ^ (termpayment))) * ex.rate
          } else if (rate.type == 0) {
            df <- exp(- PaymentRates * termpayment) * ex.rate
          }
        }

        return(df)
      })


      if (obj == "Rates") {

        Cupon <- c()
        TIR <- c()
        error <- c()

        for (j in 1:(length(PaymentDates))) {
          z <- swaps.ccs[j]
          Cupon <- c(Payments[[j]][-length(PaymentDates[[j]])],
                     Payments[[j]][length(PaymentDates[[j]])] - as.numeric(swaps[z,"prin2"])) / as.numeric(swaps[z,"C2"])
          TIR[j] <- (local.leg[[j]] - ((df[[j]][length(PaymentDates[[j]])])*(as.numeric(swaps[z,"prin2"])))) / crossprod(Cupon, df[[j]])
          error[j] <- TIR[j] - as.numeric(swaps[z,"C2"])
        }


      } else if (obj == "Price") {

        error <- lapply(1:(length(PaymentDates)), function (j) {
          (sum(Payments[[j]]*df[[j]]))- local.leg[[j]]
        })
      }

      if (is.null(Weights)) {
        Weights <- rep(1,length(error))
      } else {
        if (!length(Weights) == length(error)) {
          stop(paste0(deparse(sys.call()), ":"," Weights vector must have a length of", ": ", length(error)),
               call. = FALSE)
        }
      }

      fobj <- crossprod((unlist(error)^2),Weights)
      return (as.numeric(fobj))
    }

    end.date <- as.numeric((PaymentDates[[length(PaymentDates)]][length(PaymentDates[[length(PaymentDates)]])]
                            - Today ) / 365)
    Today      <- as.Date(analysis.date)

    FinalNode.1 <- c(rep(mean(CurrentCurve[(2):(length(CurrentCurve))]), npieces),
                     seq(from = 0,
                         to = end.date,
                         by = end.date/npieces)[-1])


    if (npieces == 1 || !is.null(piece.term)) {
      if (npieces == 1) {
        FinalNode <- FinalNode.1[1]
      } else if (!is.null(piece.term)) {
        FinalNode <- c(rep(mean(CurrentCurve[(2):(length(CurrentCurve))]), npieces), piece.term)
      }
      eqn1 <- function(FinalNode) {
        z1 <- c(FinalNode)
        return(z1)
      }

      LB <- c(rep(-0.5, npieces), rep(0, (npieces*2) - (npieces + 1)))


      CurrentCurve2 <- Rsolnp::solnp(pars = FinalNode,
                                     fun          = RateFunction.for,
                                     ineqLB       = eqn1,
                                     LB           = LB,
                                     Payments     = cashFlows,
                                     PaymentDates = PaymentDates,
                                     Today        = Today,
                                     local.leg    = local.leg
      )$pars

    } else {
      FinalNode <- matrix( data = NA, nrow = nsimul + 1, ncol = npieces*2 - 1)
      FinalNode[1,] <- FinalNode.1[-length(FinalNode.1)]
      simul.t <- matrix( data = NA, nrow = nsimul + 1, ncol = npieces - 1)
      values <- matrix( data = NA, nrow = nsimul, ncol = 1)
      optim.values <- matrix( data = NA, nrow = nsimul + 1, ncol = npieces*2 - 1)
      for (i in 1:(nsimul + 1)) {
        h <- 0
        for (j in 1:(npieces - 1)) {
          simul.t[i,j] <- runif(1,min = h, max = end.date)
          h <- simul.t[i,j]
        }
        if (i == 1) {
          FinalNode[1,] <- FinalNode.1[-length(FinalNode.1)]
        } else {
          FinalNode[i,] <- c(rep(mean(CurrentCurve[(2):(length(CurrentCurve))]), npieces),
                             simul.t[i,])
        }

        eqn1 <- function(FinalNode) {
          z1 <- c(FinalNode)
          return(z1)
        }

        LB <- c(rep(-0.5, npieces), rep(0, (npieces*2) - (npieces + 1)))


        local.optim <- Rsolnp::solnp(pars = FinalNode[i,],
                                     fun          = RateFunction.for,
                                     ineqLB       = eqn1,
                                     LB           = LB,
                                     Payments     = cashFlows,
                                     PaymentDates = PaymentDates,
                                     Today        = Today,
                                     local.leg    = local.leg
        )

        values[i] <- local.optim$values[length(local.optim$values)]
        optim.values[i,] <- local.optim$pars
      }
      CurrentCurve2 <- optim.values[which(values == min(values)),]
      if (!NCOL(CurrentCurve2) == 1) {
        CurrentCurve2 <- CurrentCurve2[1,]
      }
    }

    CurrentCurve3 <- as.numeric(CurrentCurve2 [1:npieces])
    if (npieces == 1) {
      names(CurrentCurve3)[length(CurrentCurve3)] <- end.date
    } else {
      names(CurrentCurve3) <- as.numeric(CurrentCurve2[(npieces + 1):(2 * npieces - 1)])
      names(CurrentCurve3)[length(CurrentCurve3)] <- end.date
    }
    CurrentCurve <- CurrentCurve3

  }

  if (!names(CurrentCurve)[1] == 0) {
    CurrentCurve <- c(CurrentCurve[1],CurrentCurve)
    names(CurrentCurve)[1] <- "0"
  } else if (names(CurrentCurve)[1] == 0) {
    CurrentCurve[1] <- CurrentCurve[2]
  }

  if (approximation == "linear") {
    CurrentCurve  <- c(CurrentCurve, `1000` = as.numeric(CurrentCurve[length(CurrentCurve)]))

    CurrentCurve <- c(CurrentCurve[1],CurrentCurve)
    names(CurrentCurve)[1] <- "0"
    # Returns a curve with the desired nodes by the user using a linear
    # approximation.
    ConstantCurve <- approx(
      x    = as.numeric(names(CurrentCurve)),
      y    = CurrentCurve,
      xout = nodes,
      ties = min
    )$y
    names(ConstantCurve) <- nodes
  } else if (approximation == "constant") {
    # Returns a curve with the desired nodes by the user using a constant
    # piecewise approximation.
    CurrentCurve <- c(CurrentCurve[1],CurrentCurve)
    names(CurrentCurve)[1] <- "0"

    piecewisefun <- stepfun(
      x     = as.numeric(names(CurrentCurve))[-1],
      y     = CurrentCurve,
      right = FALSE,
      f = 1,
      ties = min
    )
    ConstantCurve        <- piecewisefun(nodes)
    names(ConstantCurve) <- nodes

  }

  return(ConstantCurve)
}
