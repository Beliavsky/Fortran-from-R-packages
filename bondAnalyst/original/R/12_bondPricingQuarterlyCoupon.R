#'Calculates Present Value or the Price of the Bond paying Quarterly Coupons.
#'@details
#'The method  \code{pricingQtrlyCpnBond()} is developed to compute the Price of Bond making quarterly Coupon Payments. So, \code{pricingQtrlyCpnBond()} gives the Price of Bond making quarterly Coupon Payments for values passed to its five arguments. Here, \code{qCoupons} represents the dollar values of all the quarterly coupon payments till maturity, \code{times} is a vector of number of years ranging from 1 to any specified number of quarterly periods till maturity, that is from 1 to times(n * 4), \code{mv} is Maturity Value, \code{n} is number of years till maturity, and \code{r} is annual Market Discount Rate or Required Rate of return. The output is rounded off to two decimal places. The given examples show various ways in which the arguments can be passed to \code{pricingQtrlyCpnBond()}.
#'@param qCoupons A vector.
#'@param times A vector.
#'@param mv A number.
#'@param n A number.
#'@param r A number.
#'@return Input values to five arguments  \code{qCoupons} , \code{times}, \code{mv},\code{n} and \code{r}.
#'@importFrom Rdpack reprompt
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'Adams,J.F. & Smith,D.J.(2019). Introduction to fixed-income valuation. In CFA Program Curriculum 2020 Level I Volumes 1-6. (Vol. 5, pp. 107-151). Wiley Professional Development (P&T). ISBN 9781119593577, <https://bookshelf.vitalsource.com/books/9781119593577>
#'@examples
#'pricingQtrlyCpnBond(qCoupons=c(2,2,2,2,2,2,2,2), times=c(1,2,3,4,5,6,7,8) ,mv=100,n=2,r=0.06)
#'pricingQtrlyCpnBond(qCoupons=c(2,2,2,2,2,2,2,2), times=c(1:8) ,mv=100,n=2,r=0.06)
#'pricingQtrlyCpnBond(qCoupons=c(rep(2,8)), times=c(1:8) ,mv=100,n=2,r=0.06)
#'pricingQtrlyCpnBond(c(rep(2,8)), c(1:8) ,100,2,0.06)
#'@export
pricingQtrlyCpnBond <-function (qCoupons, times,mv,n, r)
{
  all = list(qCoupons, times, r)
  if (!is.vector(qCoupons) | !is.numeric(qCoupons))
    stop("Note: Quarterly coupon payments Amount must be a numeric vector.")
  if (!is.vector(times) | !is.numeric(times))
    stop("Note: times(1 : n*4) represents number of Quarterly coupon payments as a numeric vector.")
  if (any(times <= 0))
    stop("Note: Occurance of Quarterly oupon cannot have negative values in times.")
  if (r < 0)
    stop("Note: r is Annual market Discount rate and cannot be negative.")
  if (length(qCoupons) == 0 | length(times) == 0)
    stop("Note: Insufficient information.")
  if (length(qCoupons) != length(times))
    stop("Note: Number of Coupon Payments not equal to times(n*4).")
  pv_coupon <- sum(qCoupons/(1 + r/4)^times)
  pv_mv <- mv/(1 + r/4)^(n*4)
  bond_price <- pv_coupon + pv_mv
 (bond_price= round(bond_price, digits=2))
}


