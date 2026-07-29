#'Calculates Price of Commercial Paper.
#'@details
#'According to information provided by Adams and Smith (2019), the method \code{pricingCommercialPaper()} is developed to calculate the price of Commercial Paper for the values passed to its four arguments. Here, \code{maturityVal} is maturity value of the Commercial Paper, \code{daysToMaturity} is number of days till the maturity, \code{daysInYear} is taken to be 360, and \code{mmQuotedDiscRate} is money market quoted Discount Rate. For example, an output value of 98.56 means that the price of the commercial paper is 98.56 dollars per 100 dollars of face value.
#'@param maturityVal A number.
#'@param daysToMaturity A number.
#'@param daysInYear A number.
#'@param mmQuotedDiscRate A number.
#'@return Input values to four arguments  \code{maturityVal}, \code{daysToMaturity} ,  \code{daysInYear} and  \code{mmQuotedDiscRate}.
#'@importFrom Rdpack reprompt
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'Adams,J.F. & Smith,D.J.(2019). Introduction to fixed-income valuation. In CFA Program Curriculum 2020 Level I Volumes 1-6. (Vol. 5, pp. 107-151). Wiley Professional Development (P&T). ISBN 9781119593577, <https://bookshelf.vitalsource.com/books/9781119593577>
#'@examples
#'pricingCommercialPaper(maturityVal=100,daysToMaturity=90,daysInYear=360,mmQuotedDiscRate=0.0576)
#'@export
pricingCommercialPaper<-function (maturityVal, daysToMaturity, daysInYear,mmQuotedDiscRate){
  maturityVal <-abs(maturityVal)
  cpaper_price <- (maturityVal*(1-(daysToMaturity/daysInYear) * mmQuotedDiscRate))
  (cpaper_price = round(cpaper_price, digits=2))
}


