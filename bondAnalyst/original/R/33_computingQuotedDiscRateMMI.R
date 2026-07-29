#'Calculates Discount Rate of Money Market Instrument.
#'@details
#'It is important to understand that discount rate has a unique meaning in the money market. In general, discount rate means interest rate used to calculate a present value such as market discount rate. In the money market however, discount rate is a specific type of quoted rate (Adams & Smith, 2019).
#'Based on the given information, the method \code{computingQuotedDiscRateMMIr()} is developed to compute Discount Rate of Money Market Instrument for the values passed to its four arguments. Here, \code{pvMmi} is present value of the Money Market Instrument, \code{fvMmi} is future value of the Money Market Instrument (that is sale price and not the redemption amount), \code{daysToMaturity} is number of days till the maturity, and \code{daysInYear} is taken to be 360. For example, an output value of 0.0225 represents a quoted discount rate of 2.25 percent for an assumed 360-day year.
#'@param pvMmi A number.
#'@param fvMmi A number.
#'@param daysToMaturity A number.
#'@param daysInYear A number.
#'@return Input values to four arguments  \code{pvMmi} ,\code{fvMmi}, \code{daysToMaturity} and  \code{daysInYear}.
#'@importFrom Rdpack reprompt
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'Adams,J.F. & Smith,D.J.(2019). Introduction to fixed-income valuation. In CFA Program Curriculum 2020 Level I Volumes 1-6. (Vol. 5, pp. 107-151). Wiley Professional Development (P&T). ISBN 9781119593577, <https://bookshelf.vitalsource.com/books/9781119593577>
#'@examples
#'computingQuotedDiscRateMMI(pvMmi=9943125,fvMmi=10000000,daysToMaturity=91,daysInYear=360)
#'@export
computingQuotedDiscRateMMI <- function(pvMmi,fvMmi,daysToMaturity,daysInYear){
quoted_dr=(daysInYear/daysToMaturity)*((fvMmi-pvMmi)/fvMmi)
(quoted_dr = round(quoted_dr, digits=5))
}



