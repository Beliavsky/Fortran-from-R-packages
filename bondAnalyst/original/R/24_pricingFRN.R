#'Calculates Price of a Floating-Rate Note (FRN).
#'@details
#'Floating-rate notes are very different from a fixed-rate bond. The interest payments on a floating rate note, which often is called a floater or an FRN, are not fixed. Instead, they vary from period to period depending on the current level of a reference interest rate. The interest payments could go up or down; that is why they `“`float.`"` In principle, a floater has a stable price even in a period of volatile interest rates. With a traditional fixed-income security, interest rate volatility affects the price because the future cash flows are constant. With a floating rate note, interest rate volatility affects future interest payments. The valuation of a floating rate note needs a pricing model. Suppose that the yield spread required by investors is 40 bps over the reference rate, DM = 0.0040. The assumed discount rate per period is 0.825 percent. For N = 4, the FRN is priced at 100.196 per 100 of par value. This floater is priced at a premium above par value because the quoted margin is greater than the discount margin (Adams & Smith, 2019).
#'Based on the information provided, the method \code{pricingFRN()} is developed to compute Price of a Floating-Rate Note (FRN) for the values passed to its five arguments. Here, \code{estRtrn} is a vector of estimated returns on FRN (this does not include repayment of the principal), \code{t} is a vector of number of years ranging from 1 to any specified number of periods, \code{mv} represents Maturity Value, \code{maturityPeriod} is number of evenly spaced periods to maturity, and \code{estDisc} is assumed discount rate per period. The given examples show various ways in which the arguments can be passed to \code{pricingFRN()}.
#'@param estRtrn A vector.
#'@param t A vector.
#'@param mv A number.
#'@param maturityPeriod A number.
#'@param estDisc A number.
#'@return Input values to five arguments  \code{estRtrn} ,\code{t}, \code{mv}, \code{maturityPeriod}, and  \code{estDisc}.
#'@importFrom Rdpack reprompt
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'Adams,J.F. & Smith,D.J.(2019). Introduction to fixed-income valuation. In CFA Program Curriculum 2020 Level I Volumes 1-6. (Vol. 5, pp. 107-151). Wiley Professional Development (P&T). ISBN 9781119593577, <https://bookshelf.vitalsource.com/books/9781119593577>
#'@examples
#'pricingFRN(estRtrn=c(0.875,0.875,0.875,0.875),t=c(1,2,3,4),mv=100,maturityPeriod=4,estDisc=0.00825)
#'pricingFRN(estRtrn=c(0.875,0.875,0.875,0.875),t=c(1:4),mv=100,maturityPeriod=4,estDisc=0.00825)
#'pricingFRN(estRtrn=c(rep(0.875,4)), t=c(1:4),mv=100,maturityPeriod=4,estDisc=0.00825)
#'pricingFRN(c(rep(0.875,4)), c(1:4),100,4,0.00825)
#'@export
pricingFRN<-function (estRtrn, t,mv,maturityPeriod, estDisc)
{
  all = list(estRtrn, t, estDisc)
  if (!is.vector(estRtrn) | !is.numeric(estRtrn))
    stop("Note:estmReturn payments must be a numeric values.")
  if (any(t <= 0))
    stop("Note: occurance of estmReturn cannot have negative values in times.")
  if (estDisc < 0)
    stop("Note: estimated Discount  (estmDiscount) cannot be negative.")
  if (length(estRtrn) == 0 | length(t) == 0)
    stop("Note: Please enter Complete Information.")
  if (length(estRtrn) != length(t))
    stop("Note: Number of estmReturn Payments not equal to times.")
  pv_estRtrn <- sum(estRtrn/(1 + estDisc)^t)
  pv_mv <- mv/(1 + estDisc)^maturityPeriod
  bond_price <- pv_estRtrn + pv_mv
  (bond_price= round(bond_price, digits=3))
}


