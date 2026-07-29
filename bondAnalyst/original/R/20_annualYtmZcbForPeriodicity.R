#'Calculates annual Yield-To-Maturity (YTM) of Zero-Coupon Bond with given Price and given Maturity Value for various values of Periodicity.
#'@details
#'The periodicity of the annual market discount rate for a zero-coupon bond is arbitrary because there are no coupon payments. As in the given examples, for semi-annual compounding, the annual yield-to-maturity on the five-year, zero-coupon bond priced at 80 percent 100 of par value is stated to be 4.5130 percent. This annual rate has a periodicity of two. For quarterly compounding, the annual yield-to-maturity is stated to be 4.4880 percent. This annual rate has a periodicity of four. For monthly compounding, the annual yield-to-maturity is stated to be 4.4712 percent. This annual rate has a periodicity of 12. For annual compounding, the yield-to-maturity is stated to be 4.5640 percent. This annual rate has a periodicity of one (Adams & Smith, 2019).
#'With regard to the given information ,\code{annualYtmZcbForPeriodicity()} is developed to compute annual Yield-To-Maturity (YTM) of Zero-Coupon Bond with given Price and given Maturity Value for different values of Periodicity for values passed to its four arguments. Here, \code{maturityVal} is Maturity Value of the Bond, \code{yearsToMaturity} represents years to maturity, \code{ZCBprice} represents price of Zero-Coupon Bond, and \code{desiredPeriodicity} desired periodicity for which the YTM is to be computed. The output is rounded off to six decimal places.
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
#'annualYtmZcbForPeriodicity(maturityVal=100,yearsToMaturity=5,ZCBprice=80,desiredPeriodicity=2)
#'annualYtmZcbForPeriodicity(maturityVal=100,yearsToMaturity=5,ZCBprice=80,desiredPeriodicity=4)
#'annualYtmZcbForPeriodicity(maturityVal=100,yearsToMaturity=5,ZCBprice=80,desiredPeriodicity=12)
#'annualYtmZcbForPeriodicity(maturityVal=100,yearsToMaturity=5,ZCBprice=80,desiredPeriodicity=1)
#'@export
annualYtmZcbForPeriodicity <- function(maturityVal,yearsToMaturity, ZCBprice,desiredPeriodicity)
{
   ytm_zcb <- ((maturityVal/ZCBprice) ^ (1/(yearsToMaturity*desiredPeriodicity))-1)*desiredPeriodicity
    (ytm_zcb=(round(ytm_zcb, digits=6)))
}






