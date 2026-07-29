#'Calculates DDM value of share under the assumption that Dividends are to grow at constant rate.
#'@description
#'The simplest pattern that can be assumed in forecasting future dividends is that dividends will grow at a constant rate. So, \code{DividendN1}  is equal to \code{dividendNt} multiplied with (1 + g). Here, \code{DividendN1} expected dividend to be paid after one year and \code{dividendNt} is current dividend (Jerald E. Pinto, 2020).
#'@details
#'According to information provided by Jerald E. Pinto (2020), the method \code{shareValueDDMconstantGrowth} is developed to compute DDM value of share under the assumption that Dividends are to grow at constant rate for the values passed to its four arguments.Here, \code{dividend} is current dividend, \code{g} is rate of constant growth, \code{r} is the required rate of return on the stock ,and \code{divN} lets you make choice between D0 or D1 (that is either using current dividend (D0) or Dividend in one year (D1) as \code{dividend} in the first argument of \code{shareValueDDMconstantGrowth}).
#'@param dividend A number.
#'@param r A number.
#'@param g A number.
#'@param divN A number.
#'@return Input values to four arguments  \code{dividend}, \code{r} and  \code{g} and \code{divN}.
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'Pinto, J. E. (2020). Equity Asset Valuation (4th ed.). Wiley Professional Development (P&T). https://bookshelf.vitalsource.com/books/9781119628194
#'@examples
#'shareValueGGMconstantGrowth(dividend=1.1024,r=0.101,g=0.06,divN=1)
#'shareValueGGMconstantGrowth(dividend=1.04,r=0.101,g=0.06,divN=0)
#'@export
shareValueGGMconstantGrowth <-function(dividend,r,g,divN){
  if (divN != 1 )
    share_value <- (dividend*(1+g)/(r-g))
  else
    share_value <- (dividend/(r-g))
  (share_value = round(share_value, digits=2))
}

