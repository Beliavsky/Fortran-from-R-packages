#'Calculates Effective Annual Rate (EAR) of a Zero-Coupon Bond for various values of Periodicity.
#'@details
#'This method, \code{earZcbVariousPeriodicity()} is developed to compute an Effective Annual Rate (EAR) of a Zero-Coupon Bond for various values of Periodicity, for the values passed to its four arguments. Here, \code{maturityVal} is Maturity Value of the Bond, \code{yearsToMaturity} represents years to maturity, \code{ZCBprice} represents price of Zero-Coupon Bond, and \code{desiredPeriodicity} desired periodicity for which the Effective Annual Rate is to be computed. The output is rounded off to six decimal places.
#'@param maturityVal A number.
#'@param yearsToMaturity A number.
#'@param ZCBprice A number.
#'@param desiredPeriodicity A number.
#'@return Input values to four arguments  \code{maturityVal} , \code{yearsToMaturity},\code{ZCBprice} and  \code{desiredPeriodicity}.
#'@importFrom Rdpack reprompt
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'Adams,J.F. & Smith,D.J.(2019). Introduction to fixed-income valuation. In CFA Program Curriculum 2020 Level I Volumes 1-6. (Vol. 5, pp. 107-151). Wiley Professional Development (P&T). ISBN 9781119593577, <https://bookshelf.vitalsource.com/books/9781119593577>
#'@examples
#'earZcbVariousPeriodicity(maturityVal=100, yearsToMaturity=5, ZCBprice=80, desiredPeriodicity=1)
#'earZcbVariousPeriodicity(maturityVal=100, yearsToMaturity=5, ZCBprice=80, desiredPeriodicity=12)
#'earZcbVariousPeriodicity(maturityVal=100, yearsToMaturity=5, ZCBprice=80, desiredPeriodicity=4)
#'earZcbVariousPeriodicity(maturityVal=100, yearsToMaturity=5, ZCBprice=80, desiredPeriodicity=2)
#'@export
earZcbVariousPeriodicity <- function(maturityVal,yearsToMaturity, ZCBprice,desiredPeriodicity){
  ytm_zcb <- ((maturityVal/ZCBprice) ^ (1/(yearsToMaturity*desiredPeriodicity))-1)
  (ytm_zcb=(round(ytm_zcb, digits=6)))
}


