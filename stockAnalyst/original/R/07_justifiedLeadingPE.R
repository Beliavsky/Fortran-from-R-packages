#'Calculates Justified Leading P/E based on the Gordon Growth Model.
#'@description
#'The price-to-earnings ratio (P/E) is perhaps the most widely recognized valuation indicator, familiar to readers of newspaper financial tables and institutional research reports. Using the Gordon growth model, an expression for P/E in terms of the fundamentals can be developed. When used with forecasts of the inputs to the model, the analyst obtains a justified (fundamental) P/E ; the P/E that is fair, warranted, or justified on the basis of fundamentals (given that the valuation model is appropriate). The analyst can then state his or her view of value in terms not of the Gordon growth model value but of the justified P/E. Because P/E is so widely recognized, this method may be an effective way to communicate the analysis. Leading and trailing justified P/E expressions can be developed from the Gordon growth model. Assuming that the model can be applied to valuation of a particular stock, the dividend payout ratio is considered fixed. In leading P/E, current price is divided by earnings of next year (Jerald E. Pinto, 2020).
#'@details
#'According to information provided by Jerald E. Pinto (2020), the method \code{justifiedLeadingPE} is developed for computing Justified Leading P/E Based on the Gordon Growth Model for the values passed to its three arguments. Here, \code{rCAPM} is required rate of return based on CAPM (Capital Asset Pricing Model), \code{payoutRatio} is payout ration, and \code{g} is dividend growth rate.
#'@param rCAPM A number.
#'@param payoutRatio A number.
#'@param g A number.
#'@return Input values to three arguments  \code{rCAPM} , \code{payoutRatio} and \code{g}.
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'Pinto, J. E. (2020). Equity Asset Valuation (4th ed.). Wiley Professional Development (P&T). https://bookshelf.vitalsource.com/books/9781119628194
#'@examples
#'justifiedLeadingPE(rCAPM=0.09,payoutRatio=0.32,g=0.07)
#'justifiedLeadingPE(rCAPM=0.125,payoutRatio=0.90,g=0.03)
#'@export
justifiedLeadingPE <-function(rCAPM,payoutRatio,g){
  leading_PE<- payoutRatio/(rCAPM-g)
  (leading_PE = round(leading_PE, digits=1))
}
