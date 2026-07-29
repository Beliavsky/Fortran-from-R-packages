#'Calculates Macaulay Duration of a traditional Fixed-Rate Bond.
#'@details
#'According to Adams and Smith (2019),the duration of a bond measures the sensitivity of the bond`'`s full price (including accrued interest) to changes in bond`'`s yield-to-maturity or, more generally, to changes in benchmark interest rates. Duration estimates changes in the bond price assuming that variables other than the yield-to-maturity or benchmark rates are held constant. Most importantly, the time-to-maturity is unchanged. Therefore, duration measures the instantaneous (or, at least, same-day) change in the bond price. The accrued interest is the same, so it is the flat price that goes up or down when the full price changes. Duration is a useful measure because it represents the approximate amount of time a bond would have to be held for the market discount rate at purchase to be realized if/when there is a single change in interest rate. If the bond is held for the duration period, an increase from reinvesting coupons is offset by a decrease in price if interest rates increase and a decrease from reinvesting coupons is offset by an increase in price if interest rates decrease.
#'Macaulay duration is named after Frederick Macaulay, the Canadian economist who first wrote about this measure in a book published in 1938. This method calculates the Macaulay duration of a traditional fixed-rate bond (Adams & Smith, 2019).
#'Here, \code{n} is number of periods, \code{ytm} is yield-to-maturity, \code{coupon} is dollar value of the coupon payment,\code{maturityVal} is maturity Value,\code{daysCpnToSettle} is the number of days from the last coupon payment to the settlement date, and \code{daysCouponPeriod} is the number of days in the coupon period.
#'@param n A number.
#'@param ytm A number.
#'@param coupon A number.
#'@param maturityVal A number.
#'@param daysCpnToSettle A number.
#'@param daysCouponPeriod A number
#'@return Input values to six arguments  \code{n} , \code{ytm}, \code{coupon}, \code{maturityVal}, \code{daysCpnToSettle} and \code{daysCouponPeriod}.
#'@importFrom Rdpack reprompt
#'@author  MaheshP Kumar, Clare Matuka
#'@references
#'Adams,J.F. & Smith,D.J.(2019). Understanding Fixed‑Income Risk and Return. In CFA Program Curriculum 2020 Level I Volumes 1-6. (Vol. 5, pp. 237-299). Wiley Professional Development (P&T). ISBN 9781119593577, <https://bookshelf.vitalsource.com/books/9781119593577>
#'@examples
#'macDuration(n=10,ytm=0.104,coupon=8,maturityVal=100,daysCpnToSettle=0, daysCouponPeriod=0)
#'macDuration(n=8*2,ytm=0.06/2,coupon=3,maturityVal=100,daysCpnToSettle=57,daysCouponPeriod=180)
#'@export
macDuration <- function(n, ytm, coupon, maturityVal, daysCpnToSettle,daysCouponPeriod){
    tT=daysCpnToSettle/daysCouponPeriod
    all = list(n, ytm, coupon, maturityVal, daysCpnToSettle,daysCouponPeriod)
      if (ytm < 0)
      stop("Note: ytm cannot be negative.")
      if (daysCpnToSettle == 0 | daysCouponPeriod == 0)
      tT=0
  num <- ((n-tT)*(coupon + maturityVal))/((1+ytm)^(n-tT))
  for(i in 1:(n-1)){
    num <- num + ((i-tT)*coupon)/((1+ytm)^(i-tT))
  }
  den <- (coupon + maturityVal)/((1 + ytm)^(n-tT))
  for(i in 1:(n-1)){
    den <- den + (coupon/((1 + ytm)^(i-tT)))
  }
  macdur1 <- num/den
  (macdur1 = round(macdur1, digits=4))
  }


