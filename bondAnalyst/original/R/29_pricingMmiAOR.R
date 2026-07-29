#'Calculates Price of Money Market Instruments using Add-on Rate (AOR)
#'@details
#'Money market instruments are short-term debt securities. They range in time-to-maturity from overnight sale and repurchase agreements (repos) to one-year bank certificates of deposit. Money market instruments also include commercial paper, bankers`'` acceptances, and time deposits based on such indexes as Libor or Euribor. Money market mutual funds are a major investor in such securities. These mutual funds can invest only in certain eligible money market securities. Bank certificates of deposit, repos, and such indexes as Libor and Euribor are quoted on an add-on rate basis (Adams & Smith, 2019).
#'In light of the information given, the method \code{pricingMoneyMarketInstrUsingAOR()} is developed to compute the Price of Money Market Instruments using Add-on Rate (AOR) for the values passed to its four arguments. Here \code{maturityVal} is maturity value of the Money Market Instrument, \code{daysToMaturity} is number of days till the maturity, and \code{daysInYear} is taken to be 365, and \code{AOR} is Add-on Rate.
#'@param maturityVal A number.
#'@param daysToMaturity A number.
#'@param daysInYear A number.
#'@param AOR A number.
#'@return Input values to four arguments  \code{maturityVal} ,\code{daysToMaturity}, \code{daysInYear} and  \code{AOR}.
#'@importFrom Rdpack reprompt
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'Adams,J.F. & Smith,D.J.(2019). Introduction to fixed-income valuation. In CFA Program Curriculum 2020 Level I Volumes 1-6. (Vol. 5, pp. 107-151). Wiley Professional Development (P&T). ISBN 9781119593577, <https://bookshelf.vitalsource.com/books/9781119593577>
#'@examples
#'pricingMoneyMarketInstrUsingAOR(maturityVal=10216000,daysToMaturity=180,daysInYear=365,AOR=0.0438)
#'pricingMoneyMarketInstrUsingAOR(maturityVal=10216000,daysToMaturity=135,daysInYear=365,AOR=0.0417)
#'@export
pricingMoneyMarketInstrUsingAOR <-function (maturityVal,daysToMaturity, daysInYear,AOR){
  maturityVal <-abs(maturityVal)
  mmi_price <- (maturityVal/(1+(daysToMaturity/daysInYear)*AOR))
  mmi_price= round(mmi_price)
  (mmi_price=format(mmi_price, scientific = FALSE))
}


