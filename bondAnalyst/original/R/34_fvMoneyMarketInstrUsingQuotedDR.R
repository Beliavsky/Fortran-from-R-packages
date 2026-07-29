#'Calculates Future Value of Money Market Instruments using the given Discount Rate.
#'@details
#'According to information provided by Adams and Smith (2019), the method \code{fvMmiUsingQuotedDiscRate()} is developed to calculate Future Value of Money Market Instruments using the given Discount Rate for the values passed to its four arguments. Here, \code{pvMmi} is present value of the Money Market Instrument, \code{daysToMaturity} is number of days till the maturity, \code{daysInYear} is taken to be 360, and \code{mmQuotedDiscRate} is Money Market Quoted Discount Rate. For example, an output value of 10,000,000 dollars denotes the value of redemption amount to be paid at maturity for a Money Market Instrument that had the present value of 9,943,125 dollars.
#'@param pvMmi A number.
#'@param daysToMaturity A number.
#'@param daysInYear A number.
#'@param mmQuotedDiscRate A number.
#'@return Input values to four arguments  \code{pvMmi} , \code{daysToMaturity},\code{daysInYear}  and  \code{mmQuotedDiscRate}.
#'@importFrom Rdpack reprompt
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'Adams,J.F. & Smith,D.J.(2019). Introduction to fixed-income valuation. In CFA Program Curriculum 2020 Level I Volumes 1-6. (Vol. 5, pp. 107-151). Wiley Professional Development (P&T). ISBN 9781119593577, <https://bookshelf.vitalsource.com/books/9781119593577>
#'@examples
#'fvMmiUsingQuotedDiscRate(pvMmi=9943125,daysToMaturity=91,daysInYear=360,mmQuotedDiscRate=0.0225)
#'@export
fvMmiUsingQuotedDiscRate<-function (pvMmi,daysToMaturity, daysInYear,mmQuotedDiscRate){
  pvMmi <-abs(pvMmi)
  mmi_price <- (pvMmi/(1-(daysToMaturity/daysInYear) * mmQuotedDiscRate))
  mmi_price= round(mmi_price)
  (mmi_price=format(mmi_price, scientific = FALSE))
}

