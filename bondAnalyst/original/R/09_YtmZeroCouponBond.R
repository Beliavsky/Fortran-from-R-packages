#'Calculates the Yield-To-Maturity(YTM) of a Zero-Coupon Bond.
#'@details
#'The method  \code{ytmZeroCouponBond()} is developed to compute the Yield-To-Maturity a Zero-Coupon Bond. So, \code{ytmZeroCouponBond()} gives the Price of a Zero-Coupon Bond for values passed to its three arguments. Here, \code{maturityVal} represents the Maturity Value of the Bond, \code{n} is number of years till maturity, and \code{price} is Market Price of Zero-Coupon Bond. The output is rounded off to three decimal places. The given examples show various ways in which the arguments can be passed to  \code{ytmZeroCouponBond()}.
#'@param maturityVal A number.
#'@param n A number.
#'@param price A number.
#'@return Input values to three arguments  \code{maturityVal} , \code{n} and \code{price}.
#'@importFrom Rdpack reprompt
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'Adams,J.F. & Smith,D.J.(2019). Introduction to fixed-income valuation. In CFA Program Curriculum 2020 Level I Volumes 1-6. (Vol. 5, pp. 107-151). Wiley Professional Development (P&T). ISBN 9781119593577, <https://bookshelf.vitalsource.com/books/9781119593577>
#'@examples
#'ytmZeroCouponBond(maturityVal=100, n=60, price=22.375)
#'ytmZeroCouponBond(100, 60, 22.375)
#'@export
ytmZeroCouponBond <- function (maturityVal,n, price){
   ytm_zcb <- ((maturityVal/price) ^ (1/n)-1)
  (ytm_zcb=(round(ytm_zcb, digits=5)))
  }


