#'Valuing a share of stock using Gordon Growth Model with Negative Growth.
#'@description
#'The company named Afton Mines is a profitable venture that is expected to pay a $4.25 dividend next year. Because it is depleting its mining properties, the best estimate is that dividends will decline forever at a rate of 4 percent. The required rate of return on Afton stock is 9 percent. Compute the value of Afton share (Jerald E. Pinto, 2020).
#'@details
#'According to information provided by Jerald E. Pinto (2020), the method \code{shareValueGGMNegativeGrowth} is developed for Valuing a share of stock using Gordon Growth Model with Negative Growth for the values passed to its three arguments. Here, \code{dividend} is dollar value of the dividend, \code{r} is required rate of return and, \code{negG} represents the rate of decline in dividend.
#'@param dividend A number.
#'@param r A number.
#'@param negG A number.
#'@return Input values to three arguments  \code{dividend} , \code{r} and \code{negG}.
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'Pinto, J. E. (2020). Equity Asset Valuation (4th ed.). Wiley Professional Development (P&T). https://bookshelf.vitalsource.com/books/9781119628194
#'@examples
#'shareValueGGMNegativeGrowth(dividend=4.25,r=0.12,negG=-0.10)
#'shareValueGGMNegativeGrowth(dividend=4.25,r=0.12,negG=0.10)
#'@export
shareValueGGMNegativeGrowth <-function(dividend,r,negG){
  if (negG <= 0 )
   share_value <- (dividend/(r-negG))
  else
    share_value <- (dividend/(r+negG))
  (share_value = round(share_value, digits=2))
}

