#'Calculates Price of a Treasury bill (T-bill).
#'@details
#'As explained by Adams and Smith (2019), suppose that a 91-day Treasury bill (T-bill) with a face value of USD 10 million is quoted at a discount rate of 2.25 percent for an assumed 360-day year. Here, the maturity Value is 10,000,000 (that stand for 10 million US dollars), days to Maturity are 91, days in a year are taken as 360, and money market quoted discount rate is 0.0225. When these values are passed to the method , \code{pricingTbill}, the price of the T-bill works out to be  9,943,125 US dollars. In light of the information given, the method \code{pricingTbill} is developed to compute the Price of a Treasury bill (T-bill) for the values passed to its four arguments. Here, \code{maturityVal} is face value of the T-Bill, \code{daysToMaturity} is number of days till the maturity, \code{daysInYear} are taken to be 360, and \code{mmQuotedDiscRate} is money market quoted Discount Rate.
#'@param maturityVal A number.
#'@param daysToMaturity A number.
#'@param daysInYear A number.
#'@param mmQuotedDiscRate A number.
#'@return Input values to four arguments  \code{maturityVal} ,\code{daysToMaturity}, \code{daysInYear}, and  \code{mmQuotedDiscRate}.
#'@importFrom Rdpack reprompt
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'Adams,J.F. & Smith,D.J.(2019). Introduction to fixed-income valuation. In CFA Program Curriculum 2020 Level I Volumes 1-6. (Vol. 5, pp. 107-151). Wiley Professional Development (P&T). ISBN 9781119593577, <https://bookshelf.vitalsource.com/books/9781119593577>
#'@examples
#'pricingTbill(maturityVal=10000000,daysToMaturity=91,daysInYear=360,mmQuotedDiscRate=0.0225)
#'@export
pricingTbill <-function (maturityVal,daysToMaturity, daysInYear,mmQuotedDiscRate)
  {
    maturityVal <-abs(maturityVal)
    tbill_price <- (maturityVal*(1-(daysToMaturity/daysInYear) * mmQuotedDiscRate))
    tbill_price= round(tbill_price)
    (tbill_price=format(tbill_price, scientific = FALSE))
}




