#'#'Calculates Macaulay Duration using the Coupon Rate and Yield-To-Maturity.
#'@details
#'According to information provided by Adams and Smith (2019), the method \code{macDurationOnCouponRate()} is developed to calculate Macaulay Duration using the Coupon Rate and Yield-To-Maturity.
#'Here, \code{couponRate} is Coupon Rate (for example, 0.03 means a coupon rate of 3 percent), \code{n} is number of periods, \code{ytm} is yield-to-maturity, and \code{tT} is fraction that has the number of days from the last coupon payment to the settlement date in the numerator and the number of days in the coupon period in the denominator.
#'@param couponRate A number.
#'@param n A number.
#'@param ytm A number.
#'@param tT A number.
#'@return Input values to four arguments  \code{couponRate} , \code{n}, \code{ytm} and \code{tT}.
#'@importFrom Rdpack reprompt
#'@author MaheshP Kumar, Clare Matuka
#'@references
#'Adams,J.F. & Smith,D.J.(2019). Understanding Fixed‑Income Risk and Return. In CFA Program Curriculum 2020 Level I Volumes 1-6. (Vol. 5, pp. 237-299). Wiley Professional Development (P&T). ISBN 9781119593577, <https://bookshelf.vitalsource.com/books/9781119593577>
#'@examples
#'macDurationOnCouponRate(couponRate=0.03,n=8*2,ytm=0.06/2, tT=57/180)
#'macDurationOnCouponRate(couponRate=0.08,n=10, ytm=0.104, tT=0)
#'@export
macDurationOnCouponRate<- function(couponRate,n, ytm,  tT){
  leftside <- (1+ytm)/ytm
  rightside_num <- 1+ytm+(n*(couponRate-ytm))
  rightside_den <- (couponRate*(((1+ytm)^n)-1))+ytm
  macdur_val <- (leftside - (rightside_num/rightside_den))-tT
  (macdur_val = round(macdur_val, digits=4))
}


