#'Calculates Leading Price to Earning Multiple based on the mean of the following fiscal year (FY2 = Fiscal Year 2) forecasts.
#'@description
#'Applying the fiscal-year concept, Leading PE can be computed in two ways: first, based on the mean of the current fiscal year (FY1 = Fiscal Year 1) forecasts, for which analysts may have actual EPS in hand for some quarters; second, based on the following fiscal year (FY2 = Fiscal Year 2) forecasts, which must be based entirely on forecasts by analysts (Jerald E. Pinto, 2020).
#'@details
#'According to information provided by Jerald E. Pinto (2020), the method \code{leadingFY2PE} is developed for computing Leading PE Multiple based on the mean of the following fiscal year (FY2 = Fiscal Year 2) forecasts for the values passed to its two arguments. Here, \code{currentShPr} is the current Share Price and \code{FY2EPS} is the mean of  following fiscal year (FY2 = Fiscal Year 2) forecasts by the analysts.
#'@param  currentShPr number.
#'@param  FY2EPS number.
#'@return Input values to two arguments  \code{currentShPr} and \code{FY2EPS}.
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'Pinto, J. E. (2020). Equity Asset Valuation (4th ed.). Wiley Professional Development (P&T). https://bookshelf.vitalsource.com/books/9781119628194
#'@examples
#'leadingFY2PE(currentShPr=184.15,FY2EPS=18.35)
#'@export
leadingFY2PE<-function (currentShPr,FY2EPS){
  disp_val<- currentShPr/FY2EPS
  (disp_val = round(disp_val, digits=1))
}
