#'Calculates the Present Value of the Excess Coupon Payment resulting due to higher Coupon Rate as compared the Market Discount Rate.
#'@details
#'For Example, a Bond is trading at a premium because the coupon rate per period (6 percent) is greater than the market discount rate per period (4 percent). The excess per period is the coupon rate minus market discount rate, times the par value: (0.06 – 0.04) × 100 = +2. The present value of excess is +7.260, discounted using the required yield per period (Adams & Smith, 2019).
#'Based on this, the method  \code{pvExcessCoupon()} is developed to compute the present value of excess Coupon Payments. So,  \code{pvExcessCoupon()} gives the present value of excess Coupon Payments for values passed to its three arguments. Here, \code{couponExcess} represent the dollar value of excess coupon payment, \code{times} is a vector of number of years ranging from 1 to any specified number of years till maturity, and \code{r} is Market Discount Rate or Required Rate of return. The output is rounded off to three decimal places. The given examples show various ways in which the arguments can be passed to \code{pvExcessCoupon()}.
#'@param couponExcess A vector.
#'@param times A vector.
#'@param r A number.
#'@return Input values to three arguments  \code{couponExcess} , \code{times} and \code{r}.
#'@importFrom Rdpack reprompt
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'Adams,J.F. & Smith,D.J.(2019). Introduction to fixed-income valuation. In CFA Program Curriculum 2020 Level I Volumes 1-6. (Vol. 5, pp. 107-151). Wiley Professional Development (P&T). ISBN 9781119593577, <https://bookshelf.vitalsource.com/books/9781119593577>
#'@examples
#'pvExcessCoupon(couponExcess=c(2,2,2,2), times=c(1,2,3,4), r= 0.04)
#'pvExcessCoupon(couponExcess=c(2,2,2,2), times=c(1:4), r= 0.04)
#'pvExcessCoupon(couponExcess=c(rep(2,4)), times=c(1:4), r= 0.04)
#'pvExcessCoupon(c(rep(2,4)), c(1:4), 0.04)
#'@export
pvExcessCoupon <-function (couponExcess, times, r)
{
  all = list(couponExcess, times, r)
  if (!is.vector(couponExcess) | !is.numeric(couponExcess))
    stop("Note:coupon payments must be a numeric values.")
  if (!is.vector(times) | !is.numeric(times))
    stop("Note: times represents number of coupon payments must be 1 to number of years as a numeric vector.")
  if (length(couponExcess) != length(times))
    stop("Note: Number of Coupon Payments not equal to times.")
  pv_couponExcss <- sum(couponExcess/(1 + r)^times)
  (pv_couponExcss = round(pv_couponExcss, digits=3))
}

