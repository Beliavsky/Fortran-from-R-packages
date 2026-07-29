#'Calculates Annualized Holding Period Return of a Stock.
#'@description
#'The holding period rate of return (for short, the holding period return) is the return earned from investing in an asset over a specified time period. The specified time period is the holding period under examination, whether it is one day, two weeks, four years, or any other length of time (Jerald E. Pinto, 2020).
#'@details
#'In the example given by Jerald E. Pinto (2020), it is assumed that a share is purchased at the price of 10 dollars each, and held for three years. The company paid a per share dividend of 0.10 dollars each of the three years. So, the total per share dividend for the Holding period of 3 years comes out to be 0.30 dollars. At the end of the three years unit share price was 12 dollars. So, the total dollar value of return per share over the holding period is 2.30 dollars (0.30 as total dividend yield plus 2.00 dollars as price appreciation return). In percentage terms, HPR for 3 years is 23 percent (2.30 dollars of total return divided by 10 dollars which was unit share price in the beginning of the investment). This return of 23 percent when annualized  works out to be 7.14 percent. Based on this understanding, the method \code{annulizedHPR} is developed for computing annualized Holding Period Return of the Stock for the values passed to its four arguments. Here, \code{totalPershareDividendHP} is the total dollar value of per share dividend for the Holding period, \code{spH} is unit share price at the end of holding period, and \code{spNot} represents unit share price in the beginning of the investment and \code{n} is number of years of the holding period.
#'@param totalPershareDividendHP A number.
#'@param spH A number.
#'@param spNot A number.
#'@param n A number.
#'@return Input values to four arguments  \code{totalPershareDividendHP} , \code{spH},\code{spNot} and \code{n}.
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'Pinto, J. E. (2020). Equity Asset Valuation (4th ed.). Wiley Professional Development (P&T). https://bookshelf.vitalsource.com/books/9781119628194
#'@examples
#'annulizedHPR(totalPershareDividendHP=0.30,spH=12,spNot=10,n=3)
#'@export
annulizedHPR <-function(totalPershareDividendHP,spH,spNot,n){
  (hpr <- (totalPershareDividendHP/spNot + (spH-spNot)/spNot))
  annualizedhpr <- ((1+hpr)^(1/n)-1)
  (annualizedhpr = round(annualizedhpr, digits=4))
}
