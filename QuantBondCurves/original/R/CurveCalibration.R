#' Curve calculation
#'
#' @description Function that calculates zero coupon or instantaneous
#' forward curves for multiple analysis dates, according to historical data of internal rates
#' of return (IRR) and coupon rates of assets.
#' Extends previous market rate curves by minimizing Mean Absolute Errors (MAE) following a
#' bootstrapping recursive method. Alternatively, methodology of residual sum of squares (RSS)
#' can be employed.
#'
#'
#' @inheritParams coupon.dates
#' @inheritParams coupons
#' @inheritParams valuation.bonds
#' @inheritParams basis.curve
#'
#' @param previous.curve A matrix that stores historical curve values up to the earliest calibration date
#' on \code{serie}. The row names correspond to the dates on which the curve was calculated,
#' while the columns represent maturities as years.
#' @param market.assets A matrix containing market assets data, where the first column represents
#' the coupon rates of the assets and the second column represents their corresponding maturities as dates.
#' This input is required only if the IRR's of assets differs from their coupon rates.
#' @param serie A time series matrix that encompasses a sequence of IRR's emitted on distinct dates.
#' The columns correspond to different maturity periods, expressed in years,
#' while the row names indicate the precise dates when the rates were emitted.
#' @param npieces Number of constant or linear segments for the curve to incorporate. By
#' default \code{NULL}, and bootstrapping method is used, otherwise, minimization of
#' RSS is used.
#' @param noSpots Number of spot interest rates introduced in the \code{serie} input per row.
#' Function assumes spots are the first entries on the \code{serie} vector for every row.
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
#' @param nodes Desired output nodes of the curve.
#' @param approximation String that establish the approximation. Set
#' \code{'linear'} for a piecewise linear approximation, or \code{'constant'} for a
#' piecewise constant curve.
#' @param fwd Numeric value that determines if the desired output curve is a
#' forward or a spot curve. Set \code{0} for spot curve (default), \code{1}
#' otherwise.
#'
#' @details
#' \code{asset.type} makes reference to the following type of assets:
#' \itemize{
#'    \item "TES" for Colombian Treasury Bonds (default).
#'    \item "IBRSwaps" for swaps indexed to IBR rate.
#'    \item "LIBORSwaps" for Interest Rate Swaps (IRS) indexed to 3M LIBOR.
#'    \item "FixedIncome" for assets that are indexed to a fixed income with
#'    different frequency of payments.
#'      }
#'
#' @details If \code{npieces = NULL} uses a recursive iteration process based in
#' bootstrapping where the curve is constructed through a minimization of the MAE
#' between the dirty price of historical assets and an approximation of the theoretical
#' price of assets of same maturity. Uses the "L-BFGS-B" optimization method
#' to minimize the expected MAE. Otherwise, curve is constructed through minimization
#' of RSS where the error can be defined via price or rate.
#'
#' @return Zero Coupon curves for the corresponding analysis dates.
#' If \code{fwd = 1}, returns forward curves.
#' @export
#'
#' @examples
#' # Previous curve input
#' previous.curve <- matrix(0.04,nrow = 2,ncol = 8)
#' rownames(previous.curve) <- c("2014-01-01","2015-01-01")
#' colnames(previous.curve) <- c(0, 0.25, 0.5, 1:5)
#' # IRR's input
#' serie <- matrix(NA,nrow = 4,ncol = 6)
#' rownames(serie) <- c("2014-01-01","2015-01-01","2016-01-01","2017-01-01")
#' colnames(serie) <- c(0, 0.08333, 0.25, 0.5, 1, 2)
#' serie[1,1] <- 0.040; serie[1,2] <- 0.050; serie[1,3] <- 0.060; serie[1,4] <- 0.065
#' serie[1,5] <- 0.070; serie[1,6] <- 0.075
#' serie[2,1] <- 0.030; serie[2,2] <- 0.040; serie[2,3] <- 0.050; serie[2,4] <- 0.063
#' serie[2,5] <- 0.074; serie[2,6] <- 0.080
#' serie[3,1] <- 0.060; serie[3,2] <- 0.065; serie[3,3] <- 0.070; serie[3,4] <- 0.080
#' serie[3,5] <- 0.084; serie[3,6] <- 0.090
#' serie[4,1] <- 0.020; serie[4,2] <- 0.030; serie[4,3] <- 0.040; serie[4,4] <- 0.042
#' serie[4,5] <- 0.045; serie[4,6] <- 0.050
#' # Market Assets input
#' market.assets <- matrix(NA,nrow = 10,ncol = 2)
#' market.assets[1,1]  <- 0.040 ; market.assets[2,1]  <- 0.05
#' market.assets[3,1]  <- 0.060 ; market.assets[4,1]  <- 0.07
#' market.assets[5,1]  <- 0.080 ; market.assets[6,1]  <- 0.09
#' market.assets[7,1]  <- 0.060 ; market.assets[8,1]  <- 0.07
#' market.assets[9,1]  <- 0.075 ; market.assets[10,1] <- 0.07
#' market.assets[1,2]  <- "2016-01-01" ; market.assets[2,2]  <- "2016-02-01"
#' market.assets[3,2]  <- "2016-04-01" ; market.assets[4,2]  <- "2016-07-01"
#' market.assets[5,2]  <- "2017-01-01" ; market.assets[6,2]  <- "2017-02-01"
#' market.assets[7,2]  <- "2017-04-01" ; market.assets[8,2]  <- "2017-07-01"
#' market.assets[9,2]  <- "2018-01-01" ; market.assets[10,2] <- "2019-01-01"
#' #Calculation
#' curve.calculation(serie = serie, market.assets = market.assets,
#'                   previous.curve = previous.curve, asset.type = "TES",
#'                   freq = 1, rate.type = 1, fwd = 0,
#'                   nodes = c(0, 0.25, 0.5, 1:5), approximation = "linear")
#'
#'
curve.calculation <- function(serie, market.assets = NULL, noSpots = 1, previous.curve = NULL,
                              asset.type = "TES", freq = 1, rate.type = 1, daycount = NULL,
                              fwd = 0, npieces = NULL, obj = "Price", Weights = NULL,
                              nsimul = 1, piece.term = NULL, nodes = seq(0,10,0.001),
                              approximation = "constant") {

  # Parameter validation ----------------------------------------------------
  if (!is.null(daycount)) {
    daycount <- daycount
  } else if (asset.type == "TES") {
    daycount <- "NL/365"
  } else if (asset.type == "IBRSwaps") {
    daycount <- "ACT/360"
  } else if (asset.type == "LIBOR") {
    daycount <- "30/360"
  } else if (asset.type == "LIBORSwaps") {
    daycount <- "30/360"
  } else if (asset.type == "UVRSwaps") {
    daycount <- "ACT/365"
  }

  # Select the dates where available information of rates exist. This dates are
  # greater than the dates in the previous input curve. Hence, filters and
  # selects dates to add new curves that have a date of analysis greater than
  # the last available rate curve.

  HistoricalRates <- serie

  if (!is.null(previous.curve)) {
    HistoricalRates <- HistoricalRates[as.Date(rownames(HistoricalRates)) >
     as.Date(max(rownames(previous.curve))), ]
  }

  maturity.numeric <- colnames(HistoricalRates)

  if (nrow(HistoricalRates)) {
    # Calculate the following parameters for each resulting date from
    # the previous filter.
    FullCurve_L <- lapply(seq_along(rownames(HistoricalRates)), function(j) {

      # We re-define some variables to be taken in each iteration.
      market.assets           <- market.assets
      analysis.date   <- as.Date(rownames(HistoricalRates)[j])

      for (i in seq_along(colnames(HistoricalRates))) {

      colnames(HistoricalRates)[i] <- as.character(lubridate::`%m+%`(
        analysis.date,
        months(round(as.numeric(colnames(HistoricalRates)[i])*12),2)
      ))
      }


      if (sum(!is.na(as.Date(colnames(HistoricalRates)))) > 0) {
          # Select rates that are greater than the analysis date.
          HistoricalRates <- HistoricalRates[, as.Date(colnames(HistoricalRates)) >= analysis.date]
        }


      # Select available rates for the analysis date.
      dailyrates <- t(as.matrix(HistoricalRates[j, which(!is.na(HistoricalRates[j, ]))]))
      colnames(dailyrates) <- colnames(HistoricalRates)

      # Select only the market assets which maturity matches the available market rates.
      if (asset.type %in% c( "IBRSwaps", "LIBORSwaps", "UVRSwaps")) {
        market.assets <- NULL
      } else if (asset.type == "TES") {
          market.assets <- market.assets[as.character(market.assets[, 2]) %in% colnames(dailyrates), ]
        }


      names(dailyrates) <- maturity.numeric

      if (fwd == 0) {
        curve.calibration(
          market.assets = market.assets,
          yield.curve   = dailyrates,
          asset.type    = asset.type,
          noSpots       = noSpots,
          daycount      = daycount,
          nodes         = nodes,
          fwd           = 0,
          rate.type     = rate.type,
          analysis.date = analysis.date,
          freq          = freq,
          approximation = approximation,
          npieces       = npieces,
          obj           = obj,
          nsimul        = nsimul,
          piece.term    = piece.term,
          Weights       = Weights
        )
      } else if (fwd == 1) {
        if (rate.type == 1) {
          stop("Function requires continuous compounded rates")
        }

        # Function that calibrates and constructs a curve based on the input
        # from the market rates and available matching market assets.
        curve.calibration (
          market.assets = market.assets,
          yield.curve   = dailyrates,
          asset.type    = asset.type,
          noSpots       = noSpots,
          daycount      = daycount,
          nodes         = nodes,
          fwd           = 1,
          analysis.date = analysis.date,
          freq          = freq,
          approximation = approximation,
          npieces       = npieces,
          obj           = obj,
          nsimul        = nsimul,
          piece.term    = piece.term,
          Weights       = Weights
        )
      }
    })

    # Construct a matrix that contains the new calculated curves.
    FullCurve <- do.call(rbind, FullCurve_L)

    rownames(FullCurve) <- rownames(HistoricalRates)
    colnames(FullCurve) <- nodes
  } else {
    FullCurve <- NULL
  }

  if (!is.null(previous.curve)) {
    # Appends the the new calculated curves of the new dates of analysis to the
    # previous available curve.
    FullCurve <- rbind(previous.curve, FullCurve)
  }

  return(FullCurve)
}


#' Curve calibration
#'
#' @description Function that calibrates and returns a Zero Coupon curve based
#' on the coupon rates and IRR's of the assets. Uses the bootstrap method to find, recursively, the
#' corresponding Zero Coupon rates given by the market data. This rates are
#' then optimized by the minimization of the MAE between bond values given by
#' the constructed rates and bond market value. Alternatively, uses minimization of
#' residual sum of squares (RSS), allowing user to optimize or define an specific
#' term structure of the segments.
#'
#' @inheritParams coupon.dates
#' @inheritParams coupons
#' @inheritParams valuation.bonds
#' @inheritParams curve.calculation
#' @param yield.curve Internal rates of return of the market assets.
#' Series of rates with different maturities that adjust to the
#' time structure of the curve to calibrate and construct. \code{names(yield.curve)} needs to
#' include the numeric maturities in years of each rate, unless the input \code{market.assets} includes the maturities
#' as dates.
#' @param market.assets Matrix containing the market assets. The first column
#' contains the coupon rates of the assets and the second column represents their
#' corresponding maturities as a date. This input is only required if IRR's of assets differ from
#' coupon rate.
#' @param analysis.date Date for which the curve is going to be calibrated and,
#' thus, constructed.
#' @param noSpots Number of spot interest rates introduced in the \code{yield.curve} input. Function assumes spots are the
#' first entries on the \code{yield.curve} vector.
#'
#'
#' @details
#' \code{asset.type} makes reference to the following type of assets:
#' \itemize{
#'    \item "TES" for Colombian Treasury Bonds (default).
#'    \item "IBRSwaps" for swaps indexed to IBR rate.
#'    \item "LIBORSwaps" for Interest Rate Swaps (IRS) indexed to 3M LIBOR.
#'    \item "FixedIncome" for assets that are indexed to a fixed income with
#'    different frequency of payments.
#'      }
#'
#' @details If \code{npieces = NULL} uses a recursive iteration process based in
#' bootstrapping where the curve is constructed through a minimization of the MAE
#' between the dirty price of historical market assets and an approximation of the theoretical
#' price of assets of same maturity. Uses the "L-BFGS-B" optimization method
#' to minimize the expected MAE. Otherwise, curve is constructed through minimization
#' of RSS where the error can be defined via price or rate.
#'
#' @author Andres Galeano & Camilo Díaz
#'
#' @return Zero Coupon curve for a specific date based on historical spot rates
#' and bond structures.
#'
#' @importFrom stats approx optim stepfun runif
#' @export
#'
#' @examples
#' # Create input
#' yield.curve <- c(0.103,0.1034,0.1092, 0.1161, 0.1233, 0.1280, 0.1310, 0.1320, 0.1325, 0.1320)
#' names(yield.curve) <- c(0,0.08,0.25,0.5,1,2,3,5,7,10)
#' nodes <- seq(0,10,0.001)
#'
#' market.assets <- matrix(NA,nrow = 10,ncol = 2)
#' market.assets[1,1] <- 0.1030 ; market.assets[2,1]  <- 0.1044
#' market.assets[3,1] <- 0.1083 ; market.assets[4,1]  <- 0.1010
#' market.assets[5,1] <- 0.1120 ; market.assets[6,1]  <- 0.1130
#' market.assets[7,1] <- 0.1150 ; market.assets[8,1]  <- 0.1160
#' market.assets[9,1] <- 0.1150 ; market.assets[10,1] <- 0.13
#' market.assets[1,2] <- "2019-01-03" ; market.assets[2,2]  <- "2019-02-03"
#' market.assets[3,2] <- "2019-04-03" ; market.assets[4,2]  <- "2019-07-03"
#' market.assets[5,2] <- "2020-01-03" ; market.assets[6,2]  <- "2021-01-03"
#' market.assets[7,2] <- "2022-01-03" ; market.assets[8,2]  <- "2024-07-03"
#' market.assets[9,2] <- "2026-01-03" ; market.assets[10,2] <- "2029-01-03"
#'
#' # Function
#' curve.calibration (yield.curve = yield.curve, market.assets = market.assets,
#'                   analysis.date = "2019-01-03" , asset.type = "IBRSwaps",
#'                   freq = 4, daycount = "ACT/365", fwd = 0, nodes = nodes,
#'                   approximation = "linear")
#'
#'
curve.calibration <- function(yield.curve, market.assets = NULL, noSpots = NULL,
                              analysis.date = Sys.Date(), asset.type = "IBRSwaps",
                              freq = 4, rate.type = 0, daycount = NULL, fwd = 0,
                              npieces = NULL, obj = "Price", Weights = NULL,
                              nsimul = 1, piece.term = NULL, nodes = seq(0,10,0.001),
                              approximation = "linear") {


  #Parameter Validation
  analysis.date <- as.Date(analysis.date)

  if (!is.null(daycount)) {
    daycount <- daycount
  } else if (asset.type == "TES") {
    daycount <- "NL/365"
  } else if (asset.type == "IBRSwaps") {
    daycount <- "ACT/360"
  } else if (asset.type == "LIBORSwaps") {
    daycount <- "30/360"
  } else if (asset.type == "UVRSwaps") {
    daycount <- "ACT/365"
  }

  # How many zero spot rates are there? In the case of LIBOR, the overnight
  # rate, 3-month rate, and 6-month rate are zero coupon bond rates.
  if (length(noSpots) == 0) {
    if (asset.type == "TES") {
      noSpots <- 1
    } else if (asset.type == "IBRSwaps") {
      noSpots <- 4
    } else if (asset.type == "LIBORSwaps") {
      noSpots <- 3
    }
  } else if (noSpots <= 0) {
    stop("Calibration has to have at least one spot rate")
  }

  if (length(freq) == 1) {
    if (fwd == 0) {
      freq <- rep(freq,length(yield.curve) - noSpots)
    } else if (fwd == 1) {
      freq <- rep(freq,length(yield.curve) - 1)
    }
  } else if (!length(freq) == 1) {
    if (fwd == 0) {
    if (!length(freq) == length(yield.curve) - noSpots) {
      stop(paste0(deparse(sys.call()), ":","
                if freq as vector, then length has to be", ": ", length(yield.curve) - noSpots),
           call. = FALSE)
    }
    } else if (fwd == 1) {
      if (!length(freq) == length(yield.curve) - 1) {
        stop(paste0(deparse(sys.call()), ":","
                if freq as vector, then length has to be", ": ", length(yield.curve) - 1),
             call. = FALSE)
      }
    }
  }

  if (!is.null(npieces)) {
    if (!is.null(piece.term)) {
  if (!length(piece.term) == (npieces - 1)) {
    stop(paste0(deparse(sys.call()), ":","
                piece.term vector must have a length of", ": ", npieces - 1),
         call. = FALSE)
  }
    }
    }

 # If maturities in numeric are not given in the yield.curve input, then maturities dates must be included
 # in market.assets, in order to find this numeric maturities.

  if (is.null(names(yield.curve))) {
    names <- c()
    for (i in seq_along(yield.curve)) {
      names[i] <- discount.time(
        tinitial   = analysis.date,
        tfinal     = market.assets[i,2]
      )
    }
    names(yield.curve) <- names
  } else {
    names(yield.curve) <- names(yield.curve)
  }

 # If market.assets is null, then TIR of swaps has to be exactly its coupon rate,
 # and maturities must be provided as numeric in the yield.curve input.
 # In conclusion, either maturities are given as numeric in yield.curve, or as date
 # in market.assets.

  if (is.null(market.assets)) {
    coupon.null <- 1
    market.assets <- matrix(data = NA, nrow = length(yield.curve), ncol = 2)
    for (i in 1:(length(yield.curve))) {
      market.assets[i,2] <- as.character(lubridate::`%m+%`( analysis.date,
                                       months(round(as.numeric(names(yield.curve[i])) * 12 ,1))))
    }
  } else {
    coupon.null <- 0
  }

  if (fwd == 0) {

    # Number of pieces shall not be bigger than amount of input swaps.

    if (!is.null(npieces)) {
      if (fwd == 0) {
    if (npieces > (length(yield.curve) - noSpots)) {
      stop(paste0(deparse(sys.call()), ":",
                  " max number of pieces is", ": ", (length(yield.curve) - noSpots)),
           call. = FALSE)
    }
    } else if (fwd == 1) {
      if (npieces > (length(yield.curve) - 1)) {
        stop(paste0(deparse(sys.call()), ":",
                    " max number of pieces is", ": ", (length(yield.curve) - noSpots)),
             call. = FALSE)
      }
    }
    }

  # Calculate coupon dates for the input market assets.
  PaymentDates <- lapply((noSpots + 1):(length((yield.curve))), function(j) {
    coupon.dates(
      maturity      = market.assets[j, 2],
      asset.type    = asset.type,
      analysis.date = analysis.date,
      freq          = freq[j - noSpots]
    )$dates
  })

  # As variable rate assets depend on the market rate, we establish the
  # coupon rate needed for the coupon calculation as the market rate of
  # this assets.
  if (asset.type %in% c( "IBRSwaps", "LIBORSwaps", "UVRSwaps") || coupon.null == 1) {
    CouponRate <- as.numeric(yield.curve)
  } else {
    CouponRate <- as.numeric(market.assets[, 1])
  }

  # Calculate cash flows for the input market assets, given the coupon rate of
  # the available data.
  Payments <- lapply((noSpots + 1):length(CouponRate), function(j) {
    coupons(
      dates = PaymentDates[[j - noSpots]], coupon.rate = as.numeric(CouponRate[[j]]),
      asset.type = asset.type, analysis.date = analysis.date, freq = freq[j - noSpots],
      daycount = daycount
    )
  })

  if (asset.type %in% c("IBRSwaps", "LIBORSwaps", "UVRSwaps") || coupon.null == 1) {
    BondValues <- c()

    # Takes bond value as par.
    for (i in (noSpots + 1):nrow(market.assets)) {
      BondValues[[i - noSpots]] <- 1
    }
  } else {
    # Calculate market asset value for each type of asset.
    BondValues <- sapply((noSpots + 1):nrow(market.assets), function(j) {
      valuation.bonds(
        maturity      = market.assets[j, 2],
        coupon.rate   = as.numeric(market.assets[j, 1]),
        rate.type     = rate.type,
        rates         = yield.curve[j],
        analysis.date = analysis.date,
        asset.type    = asset.type,
        daycount      = daycount,
        freq          = freq[j - noSpots],
        spread        = 0,
        dirty         = 1
      )
    })
  }

      # Function that determines the approximation error between rates and
      # market bond prices. Builds, recursively, a curve of rates where
      # the current curve has the available zero coupon rates, and the final
      # node are the rates that are calibrated through minimization of
      # the error.
      RateFunction1 <- function(Payments_i, Dates, Today, VP,
                               CurrentCurve, FinalNode) {

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
        return(abs(sum(Payments_i / (1 + PaymentRates)^as.numeric((Dates - Today) / 365)) - VP))
      }

      RateFunction2 <- function(Payments_i, Dates, Today, VP,
                               CurrentCurve, FinalNode) {
        Curva        <- c(CurrentCurve, FinalNode)
        names(Curva) <- c(names(CurrentCurve), as.numeric((Dates[length(Dates)] - Today) / 365))
        Curva[1] <- Curva[2]
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
        return(abs(sum(Payments_i / (1 + PaymentRates(as.numeric((Dates - Today) / 365)))^as.numeric((Dates - Today) / 365)) - VP))
      }

      RateFunction3 <- function(Payments_i, Dates, Today, VP,
                               CurrentCurve, FinalNode) {
        Curva        <- c(CurrentCurve, FinalNode)
        names(Curva) <- c(names(CurrentCurve), (Dates[length(Dates)] - Today) / 365)

        PaymentRates <- approx(
          x      = round(as.numeric(names(Curva)), 6),
          y = Curva,
          xout   = round(as.numeric((Dates - Today) / 365), 6),
          method = "linear",
          ties   = min
        )$y

        return(abs(sum(Payments_i * exp((-(PaymentRates) * as.numeric((Dates - Today) / 365)))) - VP))
      }

      RateFunction4 <- function(Payments_i, Dates, Today, VP,
                               CurrentCurve, FinalNode) {
        Curva        <- c(CurrentCurve, FinalNode)
        names(Curva) <- c(names(CurrentCurve), (Dates[length(Dates)] - Today) / 365)
        Curva[1] <- Curva[2]
        PaymentRates <- stepfun(
          x     = as.numeric(names(Curva))[-1],
          y     = Curva,
          right = FALSE,
          f     = 1,
          ties   = min
        )

        return(abs(sum(Payments_i * exp((-(PaymentRates(as.numeric((Dates - Today) / 365))) * as.numeric((Dates - Today) / 365)))) - VP))
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

  CurrentCurve        <- as.numeric(yield.curve[1:noSpots])
  names(CurrentCurve) <- names(yield.curve)[1:noSpots]

  if (! names(CurrentCurve[1]) == 0) {
    CurrentCurve <- c(CurrentCurve[1],CurrentCurve)
    names(CurrentCurve)[1] <- "0"
  }

  if (asset.type == "LIBORSwaps") {
    PaymentPorcentage <- apply(
      X      = data.frame(
        analysis.date,
        lubridate::`%m+%`(analysis.date, months(round(as.numeric(names(CurrentCurve)) * 12, 0)))
        ),
      MARGIN =  1,
      function(x) {
        quantdates::day_count(tinitial = x[1], tfinal = x[2], convention = daycount)
      }
    )

    # Corrects for the LIBOR case that contains 3 zero coupon rates: 0, 3-month, and 6-month.
    CurrentCurve <- (1 + CurrentCurve * PaymentPorcentage)^(1 / as.numeric(lubridate::`%m+%`(analysis.date, months(round(as.numeric(names(CurrentCurve)) * 12, 0))) - analysis.date) * 365) - 1
    names(CurrentCurve)[1] <- "0"
    CurrentCurve[1] <- yield.curve[1]
  }

  # Function that minimizes the absolute error between a calculated bond with
  # the input rate and the market value. The output FinalNode reflects a zero
  # coupon rate equivalent that is calibrated with the market conditions.
  for (i in (noSpots + 1):nrow(market.assets)) {
    Dates      <- PaymentDates[[i - noSpots]]
    Payments_i <- Payments[[i - noSpots]]
    Today      <- as.Date(analysis.date)
    VP         <- BondValues[[i - noSpots]]

    FinalNode <- optim(
      par          = as.numeric(yield.curve[i]),
      fn           = RateFunction,
      method       = "L-BFGS-B",
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
  RateFunction.spot <- function(FinalNode,PaymentDates, Payments,
                           Today, BondValues) {

    # Constructs a time payment structure for calculating the integral of the
    # forward rates and time difference matrix for the integrals.

    # When user defines a desired number of pieces, optimization includes rates and the T of each piece. Therefore,
    # FinalNode is a vector whose first half refers to rates and the second one to Ts.
    # The T for last piece is not included in FinalNode, since it has to be the furthest
    # away maturity.


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

    CurrentCurve        <- as.numeric(yield.curve[1:noSpots])
    names(CurrentCurve) <- "0"

    if (noSpots > 1) {
      names(CurrentCurve) <- names(yield.curve)[1:noSpots]
    }

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
          df <- 1 / ((1 + PaymentRates(termpayment)) ^ (termpayment))
        } else if (rate.type == 0) {
          df <- exp(- PaymentRates(termpayment) * termpayment)
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
          df <- 1 / ((1 + PaymentRates) ^ (termpayment))
        } else if (rate.type == 0) {
          df <- exp(- PaymentRates * termpayment)
        }
      }

      return(df)
    })

 # After obtaining discount factors, an error is calculated, can be either an error when valuating the
 # swap or an error of the TIR that makes the swap value par. Weights vector allows to prioritize one
 # swap contract over another.

    if (obj == "Rates") {

      Cupon <- c()
      TIR <- c()

      for (j in 1:(length(PaymentDates))) {
        Cupon <- c(Payments[[j]][-length(PaymentDates[[j]])],
                   Payments[[j]][length(PaymentDates[[j]])] - 1) / CouponRate[noSpots + j]
        TIR[j] <- (1 - df[[j]][length(PaymentDates[[j]])]) / crossprod(Cupon, df[[j]])
      }
      error <- TIR - yield.curve[(noSpots + 1):(length(yield.curve))]

    } else if (obj == "Price") {

      error <- lapply(1:(length(PaymentDates)), function (j) {
        sum(Payments[[j]]*df[[j]])-BondValues[[j]]
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

  # Ideally, optimization finds optimal rates and Ts for each piece. The problem is that this optimization
  # problem is non-differentiable, so solution is only local for introduced intial terms.
  # It is not trivial understanding the change in fobj, when a T changes.
  # With rsolnp, initial term guess if fixed and only rates are optimized, therefore, there is an option
  # to randomly simulate terms of pieces and optimize rates. FinalNode.1 is first guess of these Ts.
  # These guess is just separating equally each term of every piece. First guess of rates is average of
  # rates found in bootstrapping method.

  FinalNode.1 <- c(rep(mean(CurrentCurve[(2):(length(CurrentCurve))]), npieces),
                   seq(from = as.numeric(names(yield.curve)[2]),
                       to = end.date,
                       by = (end.date-as.numeric(names(yield.curve)[2]))/npieces)[-1])

 # There is also the possibility to avoid optimization issue by giving the term for each piece,
 # When npiece = 1 or when user defines the pieces terms, then optimization doesn't take into
 # account terms, they are fixed.

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
                                   fun          = RateFunction.spot,
                                   ineqLB       = eqn1,
                                   LB           = LB,
                                   Payments     = Payments,
                                   PaymentDates = PaymentDates,
                                   Today        = Today,
                                   BondValues   = BondValues
    )$pars

    # Here simulations of pieces terms are created, nsimul determines the amount of simulations that
    # user decides to do. The more simulations, the more likely to find a better local solution,
    # but computational cost is high.

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
                                   fun          = RateFunction.spot,
                                   ineqLB       = eqn1,
                                   LB           = LB,
                                   Payments     = Payments,
                                   PaymentDates = PaymentDates,
                                   Today        = Today,
                                   BondValues   = BondValues
      )

      values[i] <- local.optim$values[length(local.optim$values)]
      optim.values[i,] <- local.optim$pars
    }
    CurrentCurve2 <- optim.values[which(values == min(values)),]
    if (!NCOL(CurrentCurve2) == 1) {
      CurrentCurve2 <- CurrentCurve2[1,]
    }
  }

  # CurrentCurve2 possess al the relevant information, rates and Ts as one vector
  # CurrentCurve3 just transforms CC2 vector into the desired form vector: names
  # of vector as T and vector as rates.

  CurrentCurve3 <- as.numeric(CurrentCurve2 [1:npieces])
  if (npieces == 1) {
    names(CurrentCurve3)[length(CurrentCurve3)] <- end.date
  } else {
    names(CurrentCurve3) <- as.numeric(CurrentCurve2[(npieces + 1):(2 * npieces - 1)])
    names(CurrentCurve3)[length(CurrentCurve3)] <- end.date
  }
  CurrentCurve <- c(yield.curve[1:noSpots],CurrentCurve3)
  names(CurrentCurve)[1:noSpots] <- c(names(yield.curve)[1:noSpots])
}
  } else if (fwd == 1) {

    if (!is.null(npieces)) {
    if (npieces > (length(yield.curve) - 1)) {
      stop(paste0(deparse(sys.call()), ":",
                  " max number of pieces is", ": ", (length(yield.curve) - 1)),
           call. = FALSE)
    }
    }

    # Calculates the coupon payment dates for the input market.assets.
    PaymentDates <- lapply(2:nrow(market.assets), function(j) {
      coupon.dates(
        maturity      = market.assets[j, 2],
        asset.type    = asset.type,
        analysis.date = analysis.date,
        freq          = freq[j - 1]
      )$dates
    })

    if (asset.type %in% c( "IBRSwaps", "LIBORSwaps", "UVRSwaps") || !noSpots == 1) {
      for (i in 1:(noSpots - 1)) {
        PaymentDates[[i]] <- as.Date(market.assets[i + 1,2])
      }
    }

    if (asset.type %in% c( "IBRSwaps", "LIBORSwaps", "UVRSwaps") || coupon.null == 1) {
      # As variable rate assets depend on the market rate, we establish the coupon
      # rate needed for the coupon calculation as the market rate of this assets.
      CouponRate <- as.numeric(yield.curve)
    } else {
      CouponRate <- market.assets[, 1]
    }

    # Calculates cash flows for the input market assets given the coupon rate of
    # the available data.
    if (asset.type %in% c( "IBRSwaps", "LIBORSwaps", "UVRSwaps") || !noSpots == 1) {

      Payments.spot <- lapply(2:noSpots, function(j) {

        exp(CouponRate[j]*(
          quantdates::day_count(
            tinitial   = as.Date(analysis.date),
            tfinal     = PaymentDates[[j - 1]],
            convention = daycount
          )))
      })

      Payments <- lapply((noSpots+1):length(CouponRate), function(j) {
        coupons(
          dates         = PaymentDates[[j - 1]],
          coupon.rate   = as.numeric(CouponRate[[j]]),
          asset.type    = asset.type,
          analysis.date = analysis.date,
          daycount      = daycount,
          freq          = freq[j - 1]
        )
      })

      Payments <- c(Payments.spot, Payments)
    } else {

      Payments <- lapply(2:nrow(market.assets), function(j) {
        coupons(
          dates         = PaymentDates[[j - 1]],
          coupon.rate   = as.numeric(CouponRate[[j]]),
          asset.type    = asset.type,
          analysis.date = analysis.date,
          daycount      = daycount,
          freq          = freq[j - 1]
        )
      })
    }

    if (asset.type %in% c("IBRSwaps", "LIBORSwaps", "UVRSwaps") || coupon.null == 1) {
      # Takes bond value as par.
      BondValues <- c()
      for (i in 2:nrow(market.assets)) {
        BondValues[[i - 1]] <- 1
      }
    } else {
      # Calculates bond value for each type of asset.
      BondValues <- c()
      BondValues <- sapply(2:nrow(market.assets), function(j) {
        valuation.bonds(
          maturity      = market.assets[j, 2],
          coupon.rate   = as.numeric(market.assets[j, 1]),
          rate.type     = 0,
          rates         = yield.curve[j],
          analysis.date = analysis.date,
          asset.type    = asset.type,
          daycount      = daycount,
          freq          = freq[j - 1],
          spread        = 0,
          dirty         = 1
        )
      })
    }

    RateFunction.for <- function(Payments_i, Dates, Today, VP,
                             CurrentCurve, FinalNode) {

      # Constructs a time payment structure for calculating the integral of the
      # forward rates and time difference matrix for the integrals.

      # When user defines a desired number of pieces, optimization includes rates and the T
      # of each piece. Therefore, FinalNode is a vector whose first half refers to rates and
      # the second one to Ts.
      # The T for last piece is not included in FinalNode, since it has to be the furthest
      # away maturity.


      Curva        <- c(CurrentCurve, FinalNode)
      names(Curva) <- c(names(CurrentCurve), (Dates[length(Dates)] - Today) / 365)

      nodesi      <- as.numeric(names(Curva))[-1]
      nodesi      <- t(t(replicate(1, nodesi)))
      termpayment <- as.numeric((Dates - Today) / 365)

      plazos <- matrix(data = NA, nrow = length(nodesi), ncol = length(termpayment))
      for (i in 1:(nrow(nodesi))) {
        plazos[i,] <- termpayment
      }

      if (length(termpayment) == 1) {

        plazos <- t(t(apply(cbind(nodesi, plazos), 1, function(x) {
          pmin(x[-1], x[1])
        })))

      } else{
        plazos <- t(apply(cbind(nodesi, plazos), 1, function(x) {
          pmin(x[-1], x[1])
        }))
      }
      # Numerical approximation for the integral of forward rates.
      if (length(nodesi) == 1) {
        difft <- round(t((apply(rbind(0, plazos), 2, diff))),8)
      } else if (!length(nodesi) == 1) {
      difft <- round(t(t((apply(rbind(0, plazos), 2, diff)))),8)
      }

      if (approximation == "constant") {

        # Discount factors given by the exponential of the integration of
        # forward rates.
        df <- exp(-t(difft) %*% t(t(Curva[2:length(Curva)])))
      } else if (approximation == "linear") {

        #Vector Delta R supra C
        fwd.up <- matrix(data = NA, nrow = length(nodesi), ncol = ncol(plazos))
        fwd.dw <- matrix(data = NA, nrow = length(nodesi), ncol = ncol(plazos))
        for (i in 1:(ncol(plazos))) {
          fwd.up[,i] <- t(t(approx(
            x      = round(as.numeric(names(Curva)), 6),
            y      = Curva,
            xout   = round(plazos[,i], 6),
            method = "linear",
            ties   = min
          )$y))

          for (j in 1:(nrow(plazos))) {
            fwd.dw[j,i] <- Curva[max(which(as.numeric(names(Curva)) < plazos[j,i]))]
          }
        }
        fwd.diff <- fwd.up -fwd.dw


        fwd.linear <- t(t(Curva))

        area <- matrix(data = NA, nrow = nrow(difft), ncol = ncol(difft))
        df <- c()
        for (i in 1:(ncol(difft))) {
          for (j in 1:(nrow(difft))) {
            area[j,i] <- (fwd.linear[j] * difft[j,i] + (difft[j,i]*fwd.diff[j,i])/2)
          }
          df[i] <- exp(-sum(area[,i]))
        }
      }
      return(abs(sum(Payments_i * df) - VP))
    }

    # First market rate is taken as first forward rate.
    CurrentCurve        <- as.numeric(yield.curve[1])
    names(CurrentCurve) <- "0"
    for (i in (2):nrow(market.assets)) {
      Dates      <- PaymentDates[[i - 1]]
      Payments_i <- Payments[[i - 1]]
      Today      <- as.Date(analysis.date)
      VP         <- BondValues[[i - 1]]

      # Function that minimizes the absolute error between a calculated bond with
      # the input rate and the market value. The output FinalNode reflects a zero
      # coupon rate equivalent that is calibrated with the market conditions.
      FinalNode <- optim(
        par           = as.numeric(yield.curve[i]),
        fn            = RateFunction.for,
        method       = "L-BFGS-B",
        lower        = -0.5,
        upper        = 1,
        Payments_i    = Payments_i,
        Dates         = Dates,
        Today         = Today,
        VP            = VP,
        CurrentCurve  = CurrentCurve
      )

      CurrentCurve        <- c(CurrentCurve, FinalNode$par)
      names(CurrentCurve) <- c(names(CurrentCurve)[-length(CurrentCurve)],
                               (Dates[length(Dates)] - Today) / 365)
    }
    if (!is.null(npieces)) {
   RateFunction.for2 <- function(FinalNode,PaymentDates, Payments,
                            Today, BondValues) {
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

     nodesi      <- as.numeric(names(Curva))
     nodesi      <- t(t(replicate(1, nodesi)))
     df <- lapply(1:(length(PaymentDates)), function (p) {

       termpayment <- as.numeric(PaymentDates[[p]] - Today) / 365

       plazos <- matrix(data = NA, nrow = length(nodesi), ncol = length(termpayment))
       for (i in 1:(nrow(nodesi))) {
         plazos[i,] <- termpayment
       }

       if (length(termpayment) == 1) {

         plazos <- t(t(apply(cbind(nodesi, plazos), 1, function(x) {
           pmin(x[-1], x[1])
         })))

       } else{
         plazos <- t(apply(cbind(nodesi, plazos), 1, function(x) {
           pmin(x[-1], x[1])
         }))
       }

       # Numerical approximation for the integral of forward rates.
       difft <- (apply(rbind(0, plazos), 2, diff))
       if (npieces == 1) {
         difft <- t(difft)
       }

       # Discount factors given by the exponential of the integration of
       # forward rates.
       if (approximation == "constant") {

         df <- exp(-t(difft) %*% t(t(Curva[1:length(Curva)])))

       } else if (approximation == "linear") {
         Curva.0 <- c(yield.curve[1], Curva)
         names(Curva.0)[1] <- c("0")

         #Vector Delta R supra C
         fwd.up <- matrix(data = NA, nrow = length(nodesi), ncol = ncol(plazos))
         fwd.dw <- matrix(data = NA, nrow = length(nodesi), ncol = ncol(plazos))
         for (i in 1:(ncol(plazos))) {
           fwd.up[,i] <- t(t(approx(
             x      = round(as.numeric(names(Curva.0)), 8),
             y      = Curva.0,
             xout   = round(plazos[,i], 8),
             method = "linear",
             ties   = min
           )$y))

           for (j in 1:(nrow(plazos))) {
             fwd.dw[j,i] <- Curva.0[max(which(as.numeric(names(Curva.0)) < plazos[j,i]))]
           }
         }
         fwd.diff <- fwd.up -fwd.dw


         fwd.linear <- t(t(Curva.0))

         area <- matrix(data = NA, nrow = nrow(difft), ncol = ncol(difft))
         df <- c()
         for (i in 1:(ncol(difft))) {
           for (j in 1:(nrow(difft))) {
             area[j,i] <- (fwd.linear[j] * difft[j,i] + (difft[j,i]*fwd.diff[j,i])/2)
           }
           df[i] <- exp(-sum(area[,i]))
         }
       }
       return (df)
     })

     # After obtaining discount factors, an error is calculated, can be either an error when valuating the
     # swap or an error of the TIR that makes the swap value par. Weights vector allows to prioritize one
     # swap contract over another.


     if (obj == "Rates") {

       Cupon <- c()
       TIR <- c()
       for (i in 1:(noSpots - 1)) {
         TIR[i] <- (log(1/ df[[i]]))/(quantdates::day_count(
           tinitial   = as.Date(analysis.date),
           tfinal     = PaymentDates[[i]],
           convention = daycount
         ))
       }

       for (j in (noSpots):length(PaymentDates)) {
         Cupon <- c(Payments[[j]][-length(PaymentDates[[j]])],
                    Payments[[j]][length(PaymentDates[[j]])] - 1) / CouponRate[j + 1]
         TIR[j] <- (1 - df[[j]][length(PaymentDates[[j]])]) / crossprod(Cupon, df[[j]])
       }
       error <- TIR - yield.curve[-1]

     } else if (obj == "Price") {

       error <- lapply(1:(length(PaymentDates)), function (j) {
         sum(Payments[[j]]*df[[j]])-BondValues[[j]]
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

   # Ideally, optimization finds optimal rates and Ts for each piece. The problem is that this optimization
   # problem is non-differentiable, so solution is only local for introduced intial terms.
   # It is not trivial understanding the change in fobj, when a T changes.
   # With rsolnp, initial term guess if fixed and only rates are optimized, therefore, there is an option
   # to randomly simulate terms of pieces and optimize rates. FinalNode.1 is first guess of these Ts.
   # These guess is just separating equally each term of every piece. First guess of rates is average of
   # rates found in bootstrapping method.


   FinalNode.1 <- c(rep(mean(CurrentCurve[(2):(length(CurrentCurve))]), npieces),
                    seq(from = as.numeric(names(yield.curve)[2]),
                        to = end.date,
                        by = (end.date-as.numeric(names(yield.curve)[2]))/npieces)[-1])

   # There is also the possibility to avoid optimization issue by giving the term for each piece,
   # When npiece = 1 or when user defines the pieces terms, then optimization doesn't take into
   # account terms, they are fixed.

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
                                  fun          = RateFunction.for2,
                                  ineqLB       = eqn1,
                                  LB           = LB,
                                  Payments     = Payments,
                                  PaymentDates = PaymentDates,
                                  Today        = Today,
                                  BondValues   = BondValues
     )$pars

     # Here simulations of pieces terms are created, nsimul determines the amount of simulations that
     # user decides to do. The more simulations, the more likely to find a better local solution,
     # but computational cost is high.

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
     FinalNode[i,] <- c(rep(mean(CurrentCurve[(2):(length(CurrentCurve))]), npieces), simul.t[i,])
     }

     eqn1 <- function(FinalNode) {
       z1 <- c(FinalNode)
       return(z1)
     }

     LB <- c(rep(-0.5, npieces), rep(0, (npieces*2) - (npieces + 1)))


     local.optim <- Rsolnp::solnp(pars = FinalNode[i,],
                                  fun          = RateFunction.for2,
                                  ineqLB       = eqn1,
                                  LB           = LB,
                                  Payments     = Payments,
                                  PaymentDates = PaymentDates,
                                  Today        = Today,
                                  BondValues   = BondValues
     )

     values[i] <- local.optim$values[length(local.optim$values)]
     optim.values[i,] <- local.optim$pars
   }
   CurrentCurve2 <- optim.values[which(values == min(values)),]
   if (!NCOL(CurrentCurve2) == 1) {
     CurrentCurve2 <- CurrentCurve2[1,]
   }
   }

   # CurrentCurve2 possess al the relevant information, rates and Ts as one vector
   # CurrentCurve3 just transforms CC2 vector into the desired form vector: names
   # of vector as T and vector as rates.

   CurrentCurve3 <- as.numeric(CurrentCurve2 [1:npieces])
   if (npieces == 1) {
     names(CurrentCurve3)[length(CurrentCurve3)] <- end.date
   } else {
     names(CurrentCurve3) <- as.numeric(CurrentCurve2[(npieces + 1):(2 * npieces - 1)])
     names(CurrentCurve3)[length(CurrentCurve3)] <- end.date
   }
   CurrentCurve <- CurrentCurve3
    }
  }

  if (!names(CurrentCurve)[1] == 0) {
    CurrentCurve <- c(yield.curve[1],CurrentCurve)
    names(CurrentCurve)[1] <- "0"
  }

    if (approximation == "linear") {
      CurrentCurve  <- c(CurrentCurve, `1000` = as.numeric(CurrentCurve[length(CurrentCurve)]))

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
      CurrentCurve[1] <- CurrentCurve[2]
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

  if (names(ConstantCurve)[1] == 0) {
    ConstantCurve[1] <- yield.curve[1]
  }

  return(ConstantCurve)
}
