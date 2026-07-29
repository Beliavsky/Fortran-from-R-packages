#'Calculates the Discounted Value of the the Par Value of the Bond or the amount to be paid at the maturity of the Bond using the Market Discount Rate.
#'@details
#'The par value of a bond is the amount that the issuer agrees to repay the holder of the debt instrument by the maturity date. This amount is also referred to as the principal amount, face value, or maturity value. Bonds can have any par value. Because debt instruments can have a different par value, the practice is to quote the price of a debt instrument as a percentage of its par value. A value of 100 means 100% of par value. So, for example, if a debt instrument has a par value of $1,000 and is selling for $900, it would be said to be selling at 90. If a debt instrument with a par value of 5,000 is selling for $5,500, it is said to be selling for 110 (Fabozzi, 2008).
#'Based on this, the method  \code{disMaturityValBond()} is developed to compute the present value or discounted value of the Par Value of the Bond. So,  \code{disMaturityValBond()} gives the discounted value of Par Value, Face Value, Maturity Value, or Principal amount for values passed to its three arguments. Here \code{bondMaturityVal} is amount of Par Value or Maturity Value of the Bond, \code{n} represents number of years till maturity and \code{r} is Market Discount Rate or Required Rate of return. The output is rounded off to three decimal places.
#'@param bondMaturityVal A number.
#'@param n A number.
#'@param r A number.
#'@return Input values to three arguments  \code{bondMaturityVal} , \code{n} and \code{r}.
#'@importFrom Rdpack reprompt
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'\insertRef{fifa}{bondAnalyst}
#'@examples
#'disMaturityValBond (bondMaturityVal=100,n=5,r=0.06)
#'disMaturityValBond (100,5,0.06)
#'@export
disMaturityValBond <-function (bondMaturityVal,n, r)
{
  pv_bondMaturityVal <- bondMaturityVal/(1 + r)^n
  (Pv_bondMaturityVal = round(pv_bondMaturityVal, digits=3))
}


