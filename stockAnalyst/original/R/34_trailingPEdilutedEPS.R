#'Calculates trailing Price to Earnings Multiple based on diluted Earnings Per Share (EPS).
#'@description
#'Companies are themselves required to present both basic EPS and diluted EPS. Diluted earnings per share reflects division by the number of shares that would be outstanding if holders of securities such as executive stock options, equity warrants, and convertible bonds exercised their options to obtain common stock. The diluted EPS measure also reflects the effect of such conversion on the numerator, earnings. Because companies present both EPS numbers, the analyst does not need to make the computation. Companies also typically report details of the EPS computation in a footnote to the financial statements (Jerald E. Pinto, 2020).
#'@details
#'According to information provided by Jerald E. Pinto (2020), the method \code{trailingPEdilutedEPS} is developed for computing trailing Price to Earnings Multiple based on diluted EPS for the values passed to its two arguments. Here, \code{currentShPr} is current Share Price and \code{dilutedEPS} is diluted EPS as defined in the description above.
#'@param  currentShPr number.
#'@param  dilutedEPS vector.
#'@return Input values to two arguments  \code{currentShPr} and \code{dilutedEPS}.
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'Pinto, J. E. (2020). Equity Asset Valuation (4th ed.). Wiley Professional Development (P&T). https://bookshelf.vitalsource.com/books/9781119628194
#'@examples
#'trailingPEdilutedEPS(currentShPr=596.5,dilutedEPS=15.7)
#'@export
trailingPEdilutedEPS<-function (currentShPr,dilutedEPS){
  trailing_PE<- currentShPr/dilutedEPS
  (trailing_PE = round(trailing_PE, digits=1))
}
