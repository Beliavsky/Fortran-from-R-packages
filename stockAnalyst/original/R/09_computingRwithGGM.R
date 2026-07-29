#'Calculates Required Rate of Return using the Gordon Growth Model.
#'@description
#'Under the assumption of efficient prices, the Gordon growth model has been used to estimate a stock’s required rate of return, or equivalently, the market-price-implied expected return (Jerald E. Pinto, 2020).
#'@details
#'According to information provided by Jerald E. Pinto (2020), the method \code{computingRwithGGM} is developed for computing Required Rate of Return using the Gordon Growth Model for the values passed to its three arguments. Here, \code{divN1} is dollar value of the dividend in one year, \code{g} is dividend growth rate, and \code{spNot} is current share price.
#'@param divN1 A number.
#'@param g A number.
#'@param spNot A number.
#'@return Input values to three arguments  \code{divN1} , \code{g} and \code{spNot}.
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'Pinto, J. E. (2020). Equity Asset Valuation (4th ed.). Wiley Professional Development (P&T). https://bookshelf.vitalsource.com/books/9781119628194
#'@examples
#'computingRwithGGM(divN1=2.363,g=0.055,spNot=56.60)
#'@export
computingRwithGGM <-function(divN1,g,spNot){
  r <- (divN1/spNot + g)
  (r = round(r, digits=4))
}

