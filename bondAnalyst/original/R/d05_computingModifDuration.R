#'Calculates Modified Duration statistic of a traditional Fixed-Rate Bond.
#'@details
#'According to information provided by Adams and Smith (2019), the method \code{modifDuration()} is developed to calculate Modified Duration statistic of a traditional Fixed-Rate Bond.
#'Here, \code{n} is number of periods, \code{ytm} is yield-to-maturity, \code{coupon} is dollar value of the coupon payment, \code{maturityVal} is maturity Value, \code{daysCpnToSettle} is the number of days from the last coupon payment to the settlement date, and \code{daysCouponPeriod} is the number of days in the coupon period.
#'@param n A number.
#'@param ytm A number.
#'@param coupon A number.
#'@param maturityVal A number.
#'@param daysCpnToSettle A number.
#'@param daysCouponPeriod A number
#'@return Input values to six arguments  \code{n} , \code{ytm}, \code{coupon}, \code{maturityVal}, \code{daysCpnToSettle} and \code{daysCouponPeriod}.
#'@importFrom Rdpack reprompt
#'@author MaheshP Kumar, Clare Matuka
#'@references
#'Adams,J.F. & Smith,D.J.(2019). Understanding Fixed‑Income Risk and Return. In CFA Program Curriculum 2020 Level I Volumes 1-6. (Vol. 5, pp. 237-299). Wiley Professional Development (P&T). ISBN 9781119593577, <https://bookshelf.vitalsource.com/books/9781119593577>
#'@examples
#'modifDuration(n=10, ytm=0.104, coupon=8, maturityVal=100, daysCpnToSettle=0, daysCouponPeriod=0)
#'modifDuration(n=8*2,ytm=0.06/2,coupon=3,maturityVal=100,daysCpnToSettle=57,daysCouponPeriod=180)
#'@export
modifDuration <- function(n, ytm, coupon, maturityVal, daysCpnToSettle,daysCouponPeriod){
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
  modif_duration= macdur1/(1+ytm)
  (modif_duration = round(modif_duration, digits=4))
}


