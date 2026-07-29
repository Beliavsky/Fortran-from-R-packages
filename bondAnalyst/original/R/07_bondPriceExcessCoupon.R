#'Calculates the Price of Bond making Excess Coupon Payments.
#'@details
#'For Example, if the present value of excess is +7.260, discounted using the required yield per period, then the price of the bond is 107.260 (= 100 + 7.260)(Adams & Smith, 2019).
#'Based on the information presented, the method \code{bondPriceExcessCoupon()} is developed to compute the Price of Bond making excess Coupon Payments. So,  \code{bondPriceExcessCoupon()} gives the Price of Bond making excess Coupon Payments for values passed to its three arguments. Here, \code{couponExcess} represents the dollar value of excess coupon payments, \code{times} is a vector of number of years ranging from 1 to any specified number of years till maturity, and \code{r} is Market Discount Rate or Required Rate of return. The output is rounded off to three decimal places. The given examples show various ways in which the arguments can be passed to \code{bondPriceExcessCoupon()}.
#'@param couponExcess A vector.
#'@param times A vector.
#'@param r A number.
#'@return Input values to three arguments  \code{couponExcess} , \code{times} and \code{r}.
#'@importFrom Rdpack reprompt
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'Adams,J.F. & Smith,D.J.(2019). Introduction to fixed-income valuation. In CFA Program Curriculum 2020 Level I Volumes 1-6. (Vol. 5, pp. 107-151). Wiley Professional Development (P&T). ISBN 9781119593577, <https://bookshelf.vitalsource.com/books/9781119593577>
#'bondPriceExcessCoupon(couponExcess=c(2,2,2,2), times=c(1,2,3,4), r= 0.04)
#'bondPriceExcessCoupon(couponExcess=c(2,2,2,2), times=c(1:4), r= 0.04)
#'bondPriceExcessCoupon(couponExcess=c(rep(2,4)), times=c(1:4), r= 0.04)
#'bondPriceExcessCoupon(c(rep(2,4)), c(1:4), 0.04)
#'@export
bondPriceExcessCoupon <-function (couponExcess, times, r)
{
  all = list(couponExcess, times, r)
  if (!is.vector(couponExcess) | !is.numeric(couponExcess))
    stop("Note:coupon payments must be a numeric values.")
  if (!is.vector(times) | !is.numeric(times))
    stop("Note: times represents number of coupon payments must be 1 to number of years as a numeric vector.")
  if (length(couponExcess) != length(times))
    stop("Note: Number of Coupon Payments not equal to times.")
  pv_couponExcss <- sum(couponExcess/(1 + r)^times)
  bond_PriceExCoupon <- (100+pv_couponExcss)
  (bond_PriceExCoupon = round(  bond_PriceExCoupon, digits=3))
}

