#'Calculates estimated Return on Floating-Rate Note (FRN) for a given Index, Quoted Margin, Maturity Value, and Periodicity.
#'@details
#'Floating-rate notes are very different from a fixed-rate bond. The interest payments on a floating-rate note, which often is called a floater or an FRN, are not fixed. Instead, they vary from period to period depending on the current level of a reference interest rate. The interest payments can't go up or down; that is why they `“`float.`"` The intent of an Floating-Rate Note (FRN) is to offer the investor a security that has less market price risk than a fixed-rate bond when market interest rates fluctuate. On these floaters, a specified yield spread is added to, or subtracted from, the reference rate. For example, the floater might reset its interest rate quarterly at three-month Libor plus 0.50 percent. This specified yield spread over the reference rate is called the quoted margin on the FRN. The role of the quoted margin is to compensate the investor for the difference in the credit risk of the issuer and that implied by the reference rate. For example, a company with a stronger credit rating than that of the banks included in Libor may be able to obtain a “sub-Libor” cost of borrowed funds, which results in a negative quoted margin. An AAA-rated company might be able to issue an FRN that pays three-month Libor minus 0.25 percent (Adams & Smith, 2019).
#'Based on the information provided, the method \code{returnIncomeFRN()} is developed to compute estimated Return on Floating-Rate Note (FRN) for the values passed to its four arguments. Here, \code{index} is reference rate, stated as an annual percentage rate, \code{qtdMargin} represents quoted Margin, \code{maturityVal} represents Maturity Value, and \code{periodicity} represents periodicity of the floating-rate note that is the number of payment periods per year.
#'@param index A number.
#'@param qtdMargin A number.
#'@param maturityVal A number.
#'@param periodicity A number.
#'@return Input values to four arguments  \code{index} , \code{qtdMargin},\code{maturityVal} and  \code{periodicity}.
#'@importFrom Rdpack reprompt
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'Adams,J.F. & Smith,D.J.(2019). Introduction to fixed-income valuation. In CFA Program Curriculum 2020 Level I Volumes 1-6. (Vol. 5, pp. 107-151). Wiley Professional Development (P&T). ISBN 9781119593577, <https://bookshelf.vitalsource.com/books/9781119593577>
#'@examples
#'returnIncomeFRN(index=0.0125,qtdMargin=0.0050,maturityVal=100,periodicity=2)
#'@export
returnIncomeFRN <-function (index, qtdMargin,maturityVal,periodicity){
  numeratorApproxIncome <- ((index+qtdMargin) *maturityVal)/periodicity
  (numeratorApproxIncome= round(numeratorApproxIncome, digits=6))
  }


