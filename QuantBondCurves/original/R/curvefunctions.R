#' Spot curve conversion
#'
#' @description Uses a recursive method to calculate the instantaneous forward
#' rates of a given spot curve. Calculations and formulas based on the
#' definition of forward rates where \eqn{\exp{-rT} = \exp{-\int{f(t)dt}}}.
#'
#' @param dates Term structure of rates.
#' @param spot Vector of spot rates to be converted.
#' @inheritParams curve.calculation
#'
#' @details Requires continuous rates. Recommended that the input spot curve
#' starts with maturity 0, if not, input function will approximate zero node
#' as equal to node 1 (first term structure). The time partition and available
#' data affects calculation and precision of resulting forward curve. Output
#' forward curve slightly differs from empirical curve as it calculates an
#' implied instantaneous forward curve.
#'
#' @return Instantaneous forward curve based on the input spot and the input term
#' structure.
#' @export
#'
#' @examples
#' # Inputs for calibration of spot curve
#' yield.curve <- c(0.015,0.0175, 0.0225, 0.0275, 0.0325, 0.0375,0.04,0.0425,0.045,0.0475,0.05)
#' names(yield.curve) <- c(0.5,1,2,3,4,5,6,7,8,9,10)
#' nodes <- seq(0,10,0.001)
#' # Calibration
#' spot <- curve.calibration (yield.curve = yield.curve, market.assets = NULL,
#'                            analysis.date = "2019-01-03", asset.type = "IBRSwaps",
#'                            freq = 4, rate.type = 0, fwd = 0, npieces = NULL,
#'                            obj = "Price", nodes = nodes, approximation = "linear")
#' # Spot to Forward
#' dates <- names(spot)
#' spot2forward(dates, spot, approximation = "linear")
#'

spot2forward <- function(dates, spot, approximation = "constant") {

  if(!length(dates) == length(spot)) {
    stop("Both inputs must have the same length")
  }

  dates <- as.numeric(dates)

  if (dates[1] != 0) {
    # Makes correction to rate to include a 0 rate. Takes the first rate as 0
    # rate if not included.
    dates <- c(0, as.numeric(dates))
    spot  <- c(spot[1], spot)
  }

  # Calculates the difference in time between rates.
  timediff <- diff(as.numeric(dates))


  if (approximation == "constant") {

    warning("a Zero spot constant piecewise function
        produces an identical instantaneous forward curve")

    fwd <- spot
    names(fwd) <- as.numeric(dates)

  } else if (approximation == "linear") {

    # Takes the first available rate as the forward between time 0 and time
    # of first rate.
    fwd <- c()
    fwd[1] <- spot[1]

    # - The for calculates the integral recursively from rate 1 to the last
    #   rate.
    # - Updates information based on the previous calculated forward rate.
    # - As the second rate of the spot is taken as the first forward, times
    #   and rates should be adjusted making spot to start in i + 1.
    # - This function calculates the areas or integral of the linear curve to
    #   calculate a constant piecewise forward curve.

    if (length(spot) > 1) {
      for (i in 2:(length(spot))) {
        fwd[i] <- ((spot[i]-spot[i-1])/(timediff[i-1]))*(as.numeric(dates[i-1]))+spot[i]
      }
    }

    names(fwd) <- as.numeric(dates)

  }

  return(fwd)
}

#' Forward curve conversion
#'
#' @description Uses a recursive method to calculate the implicit spot rates of
#' a given forward curve. Calculations and formulas based on the definition of
#' forward rates where \eqn{\exp{-rT} = \exp{-\int{f(t)dt}}}.
#'
#' @inheritParams curve.calculation
#' @inheritParams spot2forward
#' @param fwd Numeric vector of forward rates to be converted.
#'
#' @details Requires continuous rates. Recommended that the input forward curve
#' starts with maturity 0, if not, function will approximate zero node
#' as equal to node 1 (first term structure). Output forward
#' curve slightly differs from empirical curve as it calculates an implicit
#' forward curve.
#'
#' @return Implicit spot curve based on the input forward rates and input term
#' structure.
#' @export
#'
#' @examples
#' # Inputs for calibration of forward curve
#' yield.curve <- c(0.015,0.0175, 0.0225, 0.0275, 0.0325, 0.0375,0.04,0.0425,0.045,0.0475,0.05)
#' names(yield.curve) <- c(0.5,1,2,3,4,5,6,7,8,9,10)
#' nodes <- seq(0,10,0.5)
#' # Calibration
#' fwd <- curve.calibration (yield.curve = yield.curve, market.assets = NULL,
#'                           analysis.date = "2019-01-03", asset.type = "IBRSwaps",
#'                           freq = 4, rate.type = 0, daycount = "ACT/365",  fwd = 1,
#'                           npieces = NULL, nodes = nodes, approximation = "constant")
#' # Forward to Spot
#' dates <- names(fwd)
#' fwd2spot(dates, fwd, approximation = "constant")
#'
fwd2spot <- function(dates, fwd, approximation = "constant") {

  if(!length(dates) == length(fwd)) {
    stop("Both inputs must have the same length")
  }

  dates <- as.numeric(dates)

  if (dates[1] != 0) {
    # Makes correction to rate to include a 0 rate. Takes the first rate as
    # 0 rate if not included.
    dates <- c(0, as.numeric(dates))
    fwd   <- c(fwd[1], fwd)
  }

  timediff <- diff(as.numeric(dates))
  # Calculates the difference in time between rates.

  if (approximation == "constant") {
    # Takes the first available rate as the forward between time 0 and time
    # of first rate.
    spot    <- c()
    spot[1] <- fwd[1]

    if (length(fwd) > 1) {
      # - The for calculates the integral recursively from rate 1 to the last
      #   rate.
      # - Updates information based on the previous calculated forward rate.
      # - As the second rate of the spot is taken as the first forward, times
      #   and rates should be adjusted making spot to start in i + 1.
      # - This function calculates the area or integral of the linear curve to
      #   calculate a constant piecewise forward curve.
      for (i in 2:(length(fwd))) {
        spot[i] <- (sum(fwd[2:i] * timediff[1:(i-1)])) * (1 / dates[(i)])
      }
    }

    names(spot) <- as.numeric(dates)
  } else if (approximation == "linear") {
    fwddiff <- diff(fwd)

    spot    <- c()
    spot[1] <- fwd[1]
    spot[2] <- sum(timediff[1]*fwd[1], (fwddiff[1]*timediff[1])/2)*(1/dates[2])

    if (length(fwd) > 2) {
      # The for calculates a linear forward curve by finding recursively the
      # areas of squares and trapezoids form between spot rates.
      # - fwddiff[1:i-2] * diff [2:i-1]/2 represents the  triangles between
      #   rates following the integral formula:
      # e^(-r_2 *T_2 )=e^(-f_1 *T_1 )*e^((-f_1*(T_2-T_1 )-((f_2-f_1 )*(T_2-T_1 ))/2) )

      for (i in 3:(length(fwd))) {

        spot[i] <- sum(spot[i-1]*dates[(i-1)], timediff[i-1]*fwd[i-1],
                       ((fwddiff[i-1]*timediff[i-1])/2))*(1/dates[i])
      }

      names(spot) <- as.numeric(dates)
    }
  }

  return(spot)
}
