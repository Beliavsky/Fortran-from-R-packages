#'Calculates Earning to Price Ratio, also known as Earning Yield.
#'@description
#'If an analyst is interested in a ranking, however, one solution (applicable to any ratio involving a quantity that can be negative or zero) is the use of an inverse price ratio which is the reciprocal of the original ratio (which places price in the denominator). The use of inverse price multiples addresses the issue of consistent ranking because price is never negative. In the case of the PE, the inverse price ratio is earnings to price (EP), known as the earnings yield. Ranked by earnings yield from highest to lowest, the securities are correctly ranked from cheapest to most costly in terms of the amount of earnings one unit of currency buys (Jerald E. Pinto, 2020).
#'@details
#'According to information provided by Jerald E. Pinto (2020), the method \code{earningYieldEP} is developed for computing Earning to Price Ratio, also known as Earning Yield, for the values passed to its two arguments. Here, \code{currentShPr} is  current Share Price and \code{TTMdilutedEPS} is trailing 12 month (TTM) diluted EPS. Output of 0.0638 represents an Earning Yield of 6.38 percent.
#'@param  currentShPr number.
#'@param  TTMdilutedEPS vector.
#'@return Input values to two arguments  \code{currentShPr} and  \code{TTMdilutedEPS}.
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'Pinto, J. E. (2020). Equity Asset Valuation (4th ed.). Wiley Professional Development (P&T). https://bookshelf.vitalsource.com/books/9781119628194
#'@examples
#'earningYieldEP(currentShPr=49.19,TTMdilutedEPS=3.14)
#'@export
earningYieldEP<-function (currentShPr,TTMdilutedEPS){
  earning_yield<- TTMdilutedEPS/currentShPr
  (earning_yield = round(earning_yield, digits=4))
}
