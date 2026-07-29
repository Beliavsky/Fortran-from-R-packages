#'Calculates the Price of Bond making Deficient Coupon Payments.
#'@details
#'When the coupon rate is less than the market discount rate, the bond is priced at a discount below par value. This is so because the price of a fixed-rate bond, relative to par value, depends on the relationship of the coupon rate to the market discount rate. For this reason, when the coupon rate is less than the market discount rate, the bond is priced at a discount below par value (Adams & Smith, 2019). If a bond pays coupons at 3 percent, whereas the the market discount rate is 5 percent, then this bond has coupon deficiency of 2 percent. So for computing the price of a bond with Coupon Deficiency, present Value of the deficient Coupon Payments is computed and then is deducted from the Par Value of the bond. Of course such bond gets priced at a discount and trades below the par value for paying coupon rate less than that of market discount rate. Suppose present value of Deficient Coupon Payments comes out to be -8.425 then of the price this bond will be 91.575 (= 100 – 8.425).
#'In view of this, the method  \code{bondPriceDefCoupon()} is developed to compute the Price of Bond making Deficient Coupon Payments. So,  \code{bondPriceDefCoupon()} gives the Price of Bond making Deficient Coupon Payments for values passed to its four arguments. Here \code{parValue} represent the Maturity Value, \code{couponDeficiency} is the vector that has dollar values of Deficient Coupon Payments for periods till maturity, \code{times} is a vector of number of years ranging from 1 to any specified number of years till maturity, and \code{r} is Market Discount Rate or Required Rate of return. The output is rounded off to three decimal places.
#'@param parValue A number.
#'@param couponDeficiency A vector.
#'@param times A vector.
#'@param r A number.
#'@return Input values to four arguments \code{parValue},  \code{couponDeficiency} , \code{times} and \code{r}.
#'@importFrom Rdpack reprompt
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'Adams,J.F. & Smith,D.J.(2019). Introduction to fixed-income valuation. In CFA Program Curriculum 2020 Level I Volumes 1-6. (Vol. 5, pp. 107-151). Wiley Professional Development (P&T). ISBN 9781119593577, <https://bookshelf.vitalsource.com/books/9781119593577>
#'@examples
#'bondPriceDefCoupon(parValue=100,couponDeficiency=c(-2,-2,-2,-2,-2), times=c(1,2,3,4,5), r= 0.06)
#'bondPriceDefCoupon(parValue=100,couponDeficiency=c(-2,-2,-2,-2,-2), times=c(1:5), r= 0.06)
#'bondPriceDefCoupon(parValue=100,couponDeficiency=c(rep(-2,5)), times=c(1:5), r= 0.06)
#'bondPriceDefCoupon(100,c(rep(-2,5)), c(1:5), 0.06)
#'@export
bondPriceDefCoupon<-function (parValue, couponDeficiency, times, r)
{
  all = list(couponDeficiency, times, r)
  if (!is.vector(couponDeficiency) | !is.numeric(couponDeficiency))
    stop("Note:coupon payments must be a numeric values.")
  if (!is.vector(times) | !is.numeric(times))
    stop("Note: times represents number of coupon payments must be 1 to number of years as a numeric vector.")
  if (length(couponDeficiency) != length(times))
    stop("Note: Number of Coupon Payments not equal to times.")
  pv_couponDef <- sum(couponDeficiency/(1 + r)^times)
  bonPriceDefCoupon <- (parValue-abs(pv_couponDef))
  (bonPriceDefCoupon = round(bonPriceDefCoupon, digits=3))
}


