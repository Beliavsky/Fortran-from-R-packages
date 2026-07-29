#'Calculates the Present Value of the Deficiency as result of lower Coupon Payments as compared that of the Market.
#'@details
#'Price of a fixed-rate bond, relative to par value, depends on the relationship of the coupon rate to the market discount rate. When the coupon rate is less than the market discount rate, the bond is priced at a discount below par value. Coupon rate indicates the amount the issuer promises to pay the bondholders each year in interest. The market discount rate reflects the amount investors need to receive in interest each year in order to pay full par value for the bond (Adams & Smith, 2019). If a bond pays coupons at 4 percent, whereas the the market discount rate is 6 percent then this bond has coupon deficiency of 2 percent. So a deficiency of 2 percent also means a deficient coupon payment of 2 dollars (entered -2) for a par value of 100 dollars.
#'In light of this, the method  \code{pvCouponDeficiency()} is developed to compute the present value of the Coupon Deficiency. So,  \code{pvCouponDeficiency()} gives the discounted value of Deficient Coupon Payments for values passed to its three arguments. Here, \code{couponDeficiency} is the vector that has dollar values of Deficient Coupon Payments for periods till maturity, \code{times} is a vector of number of years ranging from 1 to any specified number of years till maturity, and \code{r} is Market Discount Rate or Required Rate of return. The output is rounded off to three decimal places.
#'@param couponDeficiency A vector.
#'@param times A vector.
#'@param r A number.
#'@return Input values to three arguments  \code{couponDeficiency} , \code{times} and \code{r}.
#'@importFrom Rdpack reprompt
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'Adams,J.F. & Smith,D.J.(2019). Introduction to fixed-income valuation. In CFA Program Curriculum 2020 Level I Volumes 1-6. (Vol. 5, pp. 107-151). Wiley Professional Development (P&T). ISBN 9781119593577, <https://bookshelf.vitalsource.com/books/9781119593577>
#'@examples
#'pvCouponDeficiency(couponDeficiency=c(-2,-2,-2,-2,-2), times=c(1,2,3,4,5), r= 0.06)
#'pvCouponDeficiency(couponDeficiency=c(-2,-2,-2,-2,-2), times=c(1:5), r= 0.06)
#'pvCouponDeficiency(couponDeficiency=c(rep(-2,5)), times=c(1:5), r= 0.06)
#'pvCouponDeficiency(c(rep(-2,5)), c(1:5), 0.06)
#'@export
pvCouponDeficiency <-function (couponDeficiency, times, r)
{
  all = list(couponDeficiency, times, r)
  if (!is.vector(couponDeficiency) | !is.numeric(couponDeficiency))
    stop("Note:coupon payments must be a numeric values.")
  if (!is.vector(times) | !is.numeric(times))
    stop("Note: times represents number of coupon payments must be 1 to number of years as a numeric vector.")
  if (length(couponDeficiency) != length(times))
    stop("Note: Number of Coupon Payments not equal to times.")
  pv_couponDef <- sum(couponDeficiency/(1 + r)^times)
  (pv_couponDef = round(pv_couponDef, digits=3))
}


