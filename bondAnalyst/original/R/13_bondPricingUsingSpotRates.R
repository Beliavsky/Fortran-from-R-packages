#'Calculate Present Value or the Price of the Bond using Spot Rates.
#'@details
#'When a fixed-rate bond is priced using the market discount rate, the same discount rate is used for each cash flow. A more fundamental approach to calculate the price of a bond is to use a sequence of market discount rates that correspond to the cash flow dates. These market discount rates are called spot rates. Spot rates are yields-to-maturity on zero-coupon bonds maturing at the date of each cash flow. Sometimes, these are called zero rates (Adams & Smith, 2019).
#'Based on the information provided, the method  \code{pricingWithSpots()} is developed to compute the Price of Bond using the Spot Rates. So,  \code{pricingWithSpots()} gives the Price of Bond using Spots for values passed to its five arguments. Here, \code{coupons} represent the dollar values of all the coupon payments, \code{spots} is a vector of spot rates, \code{times} is a vector of number of years ranging from 1 to any specified number of years till maturity, \code{mv} is Maturity Value, and \code{n} is number of years till maturity. The output is rounded off to two decimal places. The given examples show various ways in which the arguments can be passed to \code{pricingWithSpots()}.
#'Suppose that the one-year spot rate is 2 percent, the two-year spot rate is 3 percent, and the three-year spot rate is 4 percent. Then, the price of a three-year bond that makes a 5 percent annual coupon payment is 102.960.
#'@param coupons A vector.
#'@param spots A vector.
#'@param times A vector.
#'@param mv A number.
#'@param n A number.
#'@return Input values to five arguments  \code{coupons} , \code{spots},\code{times}, \code{mv}, and \code{n}.
#'@importFrom Rdpack reprompt
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'Adams,J.F. & Smith,D.J.(2019). Introduction to fixed-income valuation. In CFA Program Curriculum 2020 Level I Volumes 1-6. (Vol. 5, pp. 107-151). Wiley Professional Development (P&T). ISBN 9781119593577, <https://bookshelf.vitalsource.com/books/9781119593577>
#'@examples
#'pricingWithSpots(coupons=c(5,5,5), spots=c(0.02,0.03,0.04), times=c(1,2,3),mv=100,n=3)
#'pricingWithSpots(coupons=c(5,5,5), spots=c(0.02,0.03,0.04), times=c(1:3),mv=100,n=3)
#'pricingWithSpots(coupons=c(rep(5,3)), spots=c(0.02,0.03,0.04), times=c(1:3),mv=100,n=3)
#'pricingWithSpots(c(rep(5,3)), c(0.02,0.03,0.04), c(1:3),100,3)
#'pricingWithSpots(coupons=c(rep(10,2)), spots=c(0.01,0.02), times=c(1:2),mv=100,n=2)
#'@export
pricingWithSpots<-function (coupons, spots, times,mv,n )
{
  all = list(coupons,spots, times,mv,n )
  if (!is.vector(coupons) | !is.numeric(coupons))
    stop("coupon payments must be a numeric vector.")
  if (!is.vector(times) | !is.numeric(times))
    stop("times represents number of coupon payments must be 1 to number of years as a numeric vector.")
  if (any(times <= 0))
    stop("occurance of coupons cannot have negative values in times.")
  if (length(coupons) == 0 | length(times) == 0)
    stop("Insufficient information.")
  if (length(coupons) != length(times))
    stop("Number of Coupon Payments not equal to times.")
  pv_coupon <- sum(coupons/(1 + spots)^times)
  pv_mv <- mv/(1 + (spots[n]))^n
  bond_price <- pv_coupon + pv_mv
 ( bond_price= round(bond_price, digits=2))
}


