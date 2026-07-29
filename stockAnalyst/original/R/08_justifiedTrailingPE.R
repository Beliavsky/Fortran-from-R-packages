#'Calculates Justified Trailing P/E Based on the Gordon Growth Model.
#'@description
#'The price-to-earnings ratio (P/E) is one of the most widely recognized valuation indicator and is familiar to readers of newspaper financial tables and institutional research reports. Using the Gordon growth model, an expression for P/E in terms of the fundamentals can be developed. Because P/E is so widely recognized, this method may be an effective way to communicate the analysis. Leading and trailing justified P/E expressions can be developed from the Gordon growth model. Assuming that the model can be applied to valuation of a particular stock, the dividend payout ratio is considered fixed. In trailing P/E, current price  is divided by trailing (current year) earnings (Jerald E. Pinto, 2020).
#'@details
#'According to information provided by Jerald E. Pinto (2020), the method \code{justifiedTrailingPE} is developed for computing Justified Trailing P/E Based on the Gordon Growth Model for the values passed to its three arguments. Here, \code{rCAPM} is required rate of return based on CAPM (Capital Asset Pricing Model), \code{payoutRatio} is payout ration and \code{g} is dividend growth rate.
#'@param rCAPM A number.
#'@param payoutRatio A number.
#'@param g A number.
#'@return Input values to three arguments  \code{rCAPM} , \code{payoutRatio} and \code{g}.
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'Pinto, J. E. (2020). Equity Asset Valuation (4th ed.). Wiley Professional Development (P&T). https://bookshelf.vitalsource.com/books/9781119628194
#'@examples
#'justifiedTrailingPE(rCAPM=0.09,payoutRatio=0.32,g=0.07)
#'@export
justifiedTrailingPE <-function(rCAPM,payoutRatio,g){
  trailing_PE<- payoutRatio*(1+g)/(rCAPM-g)
  (trailing_PE = round(trailing_PE, digits=1))
}
