#'Calculates predicted value of Price to Earning Multiple based on yields on bonds.
#'@description
#'The US FED model based on a paper written by three analysts, Lander, Orphanides, and Douvogiannis in 1997, at the US Federal Reserve, predicts the return on the S&P 500 on the basis of the relationship between forecasted earnings yields and yields on bonds (as cited in Jerald E. Pinto, 2020).
#'@details
#'According to information provided by Jerald E. Pinto (2020), the method \code{predictedPEbyFEDmodel} is developed for computing predicted value of Price to Earning Multiple based on yields on bonds.
#'@param  tenYrBondYield number.
#'@return Input values to \code{tenYrBondYield}.
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'Pinto, J. E. (2020). Equity Asset Valuation (4th ed.). Wiley Professional Development (P&T). https://bookshelf.vitalsource.com/books/9781119628194
#'@examples
#'predictedPEbyFEDmodel(tenYrBondYield=0.0293)
#'@export
predictedPEbyFEDmodel<-function (tenYrBondYield){
  disp_val<- 1/tenYrBondYield
  (disp_val = round(disp_val, digits=1))
}
