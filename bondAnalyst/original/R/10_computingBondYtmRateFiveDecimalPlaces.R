#'Calculates the Yield-To-Maturity (value up to five decimal places) of the Bond paying Annual Coupons.
#'@details
#'According to Adams and Smith (2019), the yield-to-maturity is the rate of return on the bond to an investor given three critical assumptions: 1.	The investor holds the bond to maturity.2.	The issuer makes all the coupon and principal payments in the full amount on the scheduled dates. Therefore, the yield-to-maturity is the promised yield—the yield assuming the issuer does not default on any of the payments.3.	The investor can reinvest coupon payments at that same yield. This is a characteristic of an internal rate of return.
#'In view of this, the method  \code{computingBondYtmRateFiveDecimalPlaces()} computes the YTM of a Bond. So, \code{computingBondYtmRateFiveDecimalPlaces()} gives the YTM for values passed to its four arguments. Here, \code{couponPmt} represents the dollar value of Coupon Payment, \code{mv} is Maturity Value, \code{bondPv} is present value of the bond, and \code{period} is number of years. The output is rounded off to five decimal places.
#'@param couponPmt A number.
#'@param mv A number.
#'@param bondPv A number.
#'@param period A number.
#'@return Input values to four arguments  \code{couponPmt} , \code{mv}, \code{bondPv}, and \code{period}.
#'@importFrom Rdpack reprompt
#'@author  MaheshP Kumar, Clare Matuka
#'@references
#'Adams,J.F. & Smith,D.J.(2019). Introduction to fixed-income valuation. In CFA Program Curriculum 2020 Level I Volumes 1-6. (Vol. 5, pp. 107-151). Wiley Professional Development (P&T). ISBN 9781119593577, <https://bookshelf.vitalsource.com/books/9781119593577>
#'@examples
#'computingBondYtmRateFiveDecimalPlaces (couponPmt=5, mv=100, bondPv=105, period =4)
#'computingBondYtmRateFiveDecimalPlaces (5,100,105,4)
#'computingBondYtmRateFiveDecimalPlaces (1.25,100,98.175677,6)
#'computingBondYtmRateFiveDecimalPlaces (3.5,100,103.75,4)
#'computingBondYtmRateFiveDecimalPlaces (2.25,100,96.50,6)
#'computingBondYtmRateFiveDecimalPlaces (5,100,102.96,3)
#'computingBondYtmRateFiveDecimalPlaces (8,100,90.393,3)
#'computingBondYtmRateFiveDecimalPlaces (couponPmt=5.5, mv=100, bondPv=107.50, period=3)
#'computingBondYtmRateFiveDecimalPlaces (couponPmt=4.5, mv=100, bondPv=104.75, period=5)
#'computingBondYtmRateFiveDecimalPlaces (couponPmt=6, mv=100, bondPv=100.125, period=2)
#'computingBondYtmRateFiveDecimalPlaces (couponPmt=4, mv=100, bondPv=100.750, period=2)
#'@export
computingBondYtmRateFiveDecimalPlaces <- function(couponPmt, mv, bondPv, period)
{
  cfv <- mv + couponPmt
  vec <- c(-bondPv, rep(couponPmt,period-1), cfv)
  x <- polyroot(vec)
  r <- 1/x - 1
  roots <- Re(r[abs(Im(r))<1e-8])
  roots <- roots[roots>0]
  (roots=(round(roots, digits=5)))
}
