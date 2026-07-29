#'Calculates the Price of a Zero-Coupon Bond.
#'@details
#'The method  \code{pricingZeroCouponBond()} is developed to compute the Price of a Zero-Coupon Bond. So, \code{pricingZeroCouponBond()} gives the Price of a Zero-Coupon Bond for values passed to its three arguments. Here, \code{maturityVal} represents the Maturity Value of the Bond, \code{n} is number of years till maturity, and \code{r} is Market Discount Rate or Required Rate of return. The output is rounded off to three decimal places. The given examples show various ways in which the arguments can be passed to  \code{pricingZeroCouponBond()} for two different bonds.
#'@param maturityVal A number.
#'@param n A number.
#'@param r A number.
#'@return Input values to three arguments  \code{maturityVal} , \code{n} and \code{r}.
#'@importFrom Rdpack reprompt
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'Adams,J.F. & Smith,D.J.(2019). Introduction to fixed-income valuation. In CFA Program Curriculum 2020 Level I Volumes 1-6. (Vol. 5, pp. 107-151). Wiley Professional Development (P&T). ISBN 9781119593577, <https://bookshelf.vitalsource.com/books/9781119593577>
#'@examples
#'pricingZeroCouponBond(maturityVal=100, n=10, r=0.02)
#'pricingZeroCouponBond(100, 10, 0.02)
#'pricingZeroCouponBond(100, 60, 0.02527)
#'pricingZeroCouponBond(maturityVal=100, n=60, r=0.02527)
#'@export
pricingZeroCouponBond <- function (maturityVal,n, r){
  price_zcb <- (maturityVal/(1+r)^n)
  (price_zcb = (round(price_zcb, digits=3)))
}

