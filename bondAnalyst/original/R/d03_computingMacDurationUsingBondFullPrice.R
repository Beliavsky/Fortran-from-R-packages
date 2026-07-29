#'Calculates Macaulay Duration using the Full Price of the Bond and Yield-To-Maturity.
#'@details
#' According to information provided by Adams and Smith (2019), the method \code{macDurationOnFP()} is developed to calculate Macaulay Duration using the Full Price of the Bond and Yield-To-Maturity.
#' Here, \code{fp} is Full Price of the bond, \code{n} is number of periods, \code{ytm} is yield-to-maturity, \code{coupon} is dollar value of the coupon payment, \code{maturityVal} is maturity Value, \code{daysCpnToSettle} is the number of days from the last coupon payment to the settlement date, and \code{daysCouponPeriod} is the number of days in the coupon period.
#'@param fp A number.
#'@param n A number.
#'@param ytm A number.
#'@param cpn A number.
#'@param mv A number.
#'@param daysCpnToSettle A number.
#'@param daysCouponPeriod A number
#'@return Input values to seven arguments  \code{fp} , \code{n} , \code{ytm}, \code{cpn}, \code{mv}, \code{daysCpnToSettle} and \code{daysCouponPeriod}.
#'@importFrom Rdpack reprompt
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'Adams,J.F. & Smith,D.J.(2019). Understanding Fixed‑Income Risk and Return. In CFA Program Curriculum 2020 Level I Volumes 1-6. (Vol. 5, pp. 237-299). Wiley Professional Development (P&T). ISBN 9781119593577, <https://bookshelf.vitalsource.com/books/9781119593577>
#'@examples
#'macDurationOnFP(fp=100.9404,n=8*2,ytm=0.06/2,cpn=3,mv=100,daysCpnToSettle=57,daysCouponPeriod=180)
#'macDurationOnFP(fp=85.5031,n=10, ytm=0.104, cpn=8, mv=100,daysCpnToSettle=0,daysCouponPeriod=0)
#'@export
macDurationOnFP <- function(fp,n, ytm, cpn, mv, daysCpnToSettle,daysCouponPeriod){
  tT=daysCpnToSettle/daysCouponPeriod
  all = list(n, ytm, cpn, mv, daysCpnToSettle,daysCouponPeriod)
  if (ytm < 0)
    stop("Note: ytm cannot be negative.")
  if (daysCpnToSettle == 0 | daysCouponPeriod == 0)
    tT=0
  num <- ((n-tT)*(cpn + mv))/((1+ytm)^(n-tT))
  for(i in 1:(n-1)){
    num <- num + ((i-tT)*cpn)/((1+ytm)^(i-tT))
  }
  macDuration <- num/fp
  (macDuration = round(macDuration, digits=4))
}

