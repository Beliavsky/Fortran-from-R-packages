#'Calculates Future Value of Money Market Instrument using Add-on Rate (AOR)
#'@details
#'As explained by Adams and Smith (2019), suppose that a Canadian pension fund buys a 180-day bankers`'` acceptance (BA) with a quoted add-on rate of 4.38 percent for a 365-day year. If the initial principal amount is 10 million dollars, the redemption amount due at maturity is found is  10,216,000 dollars which is calculated as the principal of 10 million dollars plus interest of 216,000 dollars. The interest is calculated as the principal times the fraction of the year times the annual add-on rate. It is added to the principal to determine the redemption amount.
#'In light of the available information, the method \code{fvMoneyMarketInstrUsingAOR()} is developed to compute the Future Value of Money Market Instrument using Add-on Rate (AOR) for the values passed to its four arguments. Here, \code{pvMmi} is present value of the Money Market Instrument, \code{daysToMaturity} is number of days till the maturity, \code{daysInYear} is taken to be 365, and \code{AOR} is Add-on Rate.
#'@param pvMmi A number.
#'@param daysToMaturity A number.
#'@param daysInYear A number.
#'@param AOR A number.
#'@return Input values to four arguments  \code{pvMmi} ,\code{daysToMaturity}, \code{daysInYear} and  \code{AOR}.
#'@importFrom Rdpack reprompt
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'Adams,J.F. & Smith,D.J.(2019). Introduction to fixed-income valuation. In CFA Program Curriculum 2020 Level I Volumes 1-6. (Vol. 5, pp. 107-151). Wiley Professional Development (P&T). ISBN 9781119593577, <https://bookshelf.vitalsource.com/books/9781119593577>
#'@examples
#'fvMoneyMarketInstrUsingAOR(pvMmi=10000000,daysToMaturity=180,daysInYear=365,AOR=0.0438)
#'@export
fvMoneyMarketInstrUsingAOR <-function (pvMmi,daysToMaturity, daysInYear,AOR){
  pvMmi <-abs(pvMmi)
  mmi_fv <- (pvMmi+(pvMmi*(daysToMaturity/daysInYear)*AOR))
  mmi_fv= round(mmi_fv)
  (mmi_fv=format(mmi_fv, scientific = FALSE))
}


