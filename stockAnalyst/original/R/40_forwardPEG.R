#'Calculates PE-to-growth (PEG) ratio.
#'@description
#'A metric that appears to address the impact of earnings growth on PE is the PE-to-growth (PEG) ratio. PEG is calculated as the  PE of the stock divided by the expected earnings growth rate (in percentage terms). The ratio, in effect, is a calculation of  PE per percentage point of expected growth. Stocks with lower PEGs are more attractive than stocks with higher PEGs, all else being equal. Some consider that a PEG ratio less than 1 is an indicator of an attractive value level. PEG is useful but must be used with care ad PEG assumes a linear relationship between PE and growth. The model for PE in terms of the DDM shows that, in theory, the relationship is not linear (Jerald E. Pinto, 2020).
#'@details
#'According to information provided by Jerald E. Pinto (2020), the method \code{forwardPEG} is developed for computing PE-to-growth (PEG) ratio for the values passed to its two arguments. Here, \code{leadingPE} is leading PE Multiple and \code{percentEPSgrowth} is five-year EPS growth forecast (in percentage terms).
#'@param  leadingPE number.
#'@param  percentEPSgrowth number.
#'@return Input values to two arguments \code{leadingPE} and \code{percentEPSgrowth}.
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'Pinto, J. E. (2020). Equity Asset Valuation (4th ed.). Wiley Professional Development (P&T). https://bookshelf.vitalsource.com/books/9781119628194
#'@examples
#'forwardPEG(leadingPE=43.97,percentEPSgrowth=25.30)
#'@export
forwardPEG<-function (leadingPE,percentEPSgrowth){
  disp_val<- leadingPE/percentEPSgrowth
  (disp_val = round(disp_val, digits=2))
}
