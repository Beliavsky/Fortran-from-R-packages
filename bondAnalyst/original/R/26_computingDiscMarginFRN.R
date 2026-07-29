#'Calculates Discount Margin of a Floating-Rate Note (FRN).
#'@details
#'According to information provided by Adams and Smith (2019), the method \code{discMarginFRN()} is developed to compute Discount Margin of a Floating-Rate Note (FRN) for the values passed to its six arguments. Here, \code{index} is reference rate, stated as an annual percentage rate,\code{estRtrn} is a vector of estimated returns on FRN (this does not include repayment of the principal), \code{mvFRN} represents Maturity Value, \code{priceFRN} price of FRN, \code{maturityYears} is number of years to Maturity, and \code{periodicity} is periodicity of interest payments. The given examples show various ways in which the arguments can be passed to \code{discMarginFRN()}.
#'@param index A number.
#'@param estRtrn A number.
#'@param mvFRN A number.
#'@param priceFRN A number.
#'@param maturityYears A number.
#'@param periodicity A number.
#'@return Input values to six arguments  \code{index} ,\code{estRtrn}, \code{mvFRN}, \code{priceFRN}, \code{maturityYears}, and  \code{periodicity}.
#'@importFrom Rdpack reprompt
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'Adams,J.F. & Smith,D.J.(2019). Introduction to fixed-income valuation. In CFA Program Curriculum 2020 Level I Volumes 1-6. (Vol. 5, pp. 107-151). Wiley Professional Development (P&T). ISBN 9781119593577, <https://bookshelf.vitalsource.com/books/9781119593577>
#'@examples
#'discMarginFRN(index=0.0200,estRtrn=0.8125,mvFRN=100,priceFRN=98,maturityYears=4,periodicity=4)
#'discMarginFRN(0.0200,0.8125, 100, 98, 4, 4)
#'@export
discMarginFRN <- function(index, estRtrn, mvFRN, priceFRN, maturityYears, periodicity){
  maturityPeriods=maturityYears*periodicity
  cfv <- mvFRN + estRtrn
  vec <- c(-priceFRN, rep(estRtrn,maturityPeriods-1), cfv)
  x <- polyroot(vec)
  r <- 1/x - 1
  roots <- Re(r[abs(Im(r))<1e-8])
  roots <- roots[roots>0]
  roots=(round(roots, digits=6))
  roots=as.numeric(roots)
  (discMargin =  round(((roots * periodicity)-index),digits=6))
}


