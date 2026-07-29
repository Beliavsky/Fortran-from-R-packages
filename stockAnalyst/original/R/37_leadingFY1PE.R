#'Calculates Leading Price to Earning Multiple based on the mean of the current fiscal year (FY1 = Fiscal Year 1) forecasts.
#'@description
#'Applying the fiscal-year concept, Leading PE can be computed in two ways: first, based on the mean of the current fiscal year (FY1 = Fiscal Year 1) forecasts, for which analysts may have actual EPS in hand for some quarters; second, based on the following fiscal year (FY2 = Fiscal Year 2) forecasts, which must be based entirely on forecasts by analysts (Jerald E. Pinto, 2020).
#'@details
#'According to information provided by Jerald E. Pinto (2020), the method \code{leadingFY1PE} is developed for computing Leading PE Multiple based on the mean of the current fiscal year (FY1)  for the values passed to its two arguments. Here, \code{currentShPr} is the current Share Price and \code{FY1EPS} is the mean of the current fiscal year (FY1 = Fiscal Year 1) forecasts, for which analysts may have actual EPS in hand for some quarters.
#'@param  currentShPr number.
#'@param  FY1EPS number.
#'@return Input values to two arguments  \code{currentShPr} and \code{FY1EPS}.
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'Pinto, J. E. (2020). Equity Asset Valuation (4th ed.). Wiley Professional Development (P&T). https://bookshelf.vitalsource.com/books/9781119628194
#'@examples
#'leadingFY1PE(currentShPr=184.15,FY1EPS=16.19)
#'@export
leadingFY1PE<-function (currentShPr,FY1EPS){
  disp_val<- currentShPr/FY1EPS
  (disp_val = round(disp_val, digits=1))
}
