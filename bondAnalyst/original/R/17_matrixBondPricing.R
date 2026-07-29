#'Calculate Present Value or the Price of illiquid Bond using Matrix Method.
#'@details
#'Some fixed-rate bonds are not actively traded. Therefore, there is no market price available to calculate the rate of return required by investors. The same problem occurs for bonds that are not yet issued. In these situations, it is common to estimate the market discount rate and price based on the quoted or flat prices of more frequently traded comparable bonds. These comparable bonds have similar times-to-maturity, coupon rates, and credit quality. This estimation process is called matrix pricing. Matrix pricing also is used in underwriting new bonds to get an estimate of the required yield spread over the benchmark rate. The benchmark rate typically is the yield-to-maturity on a government bond having the same, or close to the same, time-to-maturity (Adams & Smith, 2019).
#'Based on the information provided, \code{matrixMethod()} is developed to compute the Price of Bond using Matrix Method. So,  \code{matrixMethod()} gives the Price of Bond using returns of two frequently traded comparable bonds for values passed to its six arguments. Here, \code{couponPmt} represents the dollar values of all the coupon payments, \code{times} is a vector of number of years ranging from 1 to any specified number of years till maturity, \code{maturityVal} is Maturity Value, \code{n} is number of years till maturity, \code{r1} return on first comparable bond, and \code{r2} is return on second comparable bond. The function computes simple average of the returns of two comparable bonds with \code{r=(r1+r2)/2}, that is being used for Bond pricing. The output is rounded off to three decimal places.
#'@param couponPmt A vector.
#'@param times A vector.
#'@param maturityVal A number.
#'@param n A number.
#'@param r1 A number.
#'@param r2 A number
#'@return Input values to six arguments  \code{couponPmt} ,\code{times}, \code{maturityVal},\code{n}, \code{r1} and \code{r2}.
#'@importFrom Rdpack reprompt
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'Adams,J.F. & Smith,D.J.(2019). Introduction to fixed-income valuation. In CFA Program Curriculum 2020 Level I Volumes 1-6. (Vol. 5, pp. 107-151). Wiley Professional Development (P&T). ISBN 9781119593577, <https://bookshelf.vitalsource.com/books/9781119593577>
#'@examples
#'matrixMethod(couponPmt=c(4.5,4.5,4.5,4.5),times=c(1:4),maturityVal=100,n=4,r1=0.02856,r2=0.03449)
#'@export
matrixMethod<- function(couponPmt, times, maturityVal,n,r1,r2)
{
  r=(r1+r2)/2
  pv_coupon <- sum(couponPmt/(1 + r)^times)
  pv_maturityVal <- maturityVal/(1 + r)^n
  bond_price <- pv_coupon + pv_maturityVal
 (bond_price= round(bond_price, digits=3))
  }



