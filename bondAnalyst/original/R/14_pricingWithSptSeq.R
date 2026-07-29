#'Calculate Present Value or the Price of the Bond using two different Sequences of Spot Rates.
#'@details
#'The method  \code{pricingWithSptSeq()} is developed to compute the Price of Bond using Sequences of the Spot Rates. So,  \code{pricingWithSptSeq()} gives the Price of Bond using two different series or sequences of Spot Rates for values passed to its five arguments. Here \code{cpns} represents the dollar values of all the coupon payments,\code{sp} is a vector of spot rates, \code{t} is a vector of number of years ranging from 1 to any specified number of years till maturity, \code{mv} is Maturity Value, and \code{n} is number of years till maturity. The output is rounded off to three decimal places. The given examples show various ways in which the arguments can be passed to \code{pricingWithSptSeq()}.
#'@param cpns A vector.
#'@param sp A vector.
#'@param t A vector.
#'@param mv A number.
#'@param n A number.
#'@return Input values to five arguments  \code{cpns} , \code{sp},\code{t}, \code{mv}, and \code{n}.
#'@importFrom Rdpack reprompt
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'Adams,J.F. & Smith,D.J.(2019). Introduction to fixed-income valuation. In CFA Program Curriculum 2020 Level I Volumes 1-6. (Vol. 5, pp. 107-151). Wiley Professional Development (P&T). ISBN 9781119593577, <https://bookshelf.vitalsource.com/books/9781119593577>
#'@examples
#'pricingWithSptSeq (cpns=c(3,3,3,3), sp=c(0.0408,0.0401,0.0370,0.0350), t=c(1,2,3,4),mv=100,n=4)
#'pricingWithSptSeq (cpns=c(3,3,3,3), sp=c(0.0039,0.0140,0.0250,0.0360), t=c(1,2,3,4),mv=100,n=4)
#'pricingWithSptSeq (c(3,3,3,3), c(0.0408,0.0401,0.0370,0.0350), c(1,2,3,4),100,4)
#'pricingWithSptSeq (c(3,3,3,3), c(0.0039,0.0140,0.0250,0.0360), c(1,2,3,4),100,4)
#'pricingWithSptSeq (cpns=c(3,3,3,3), sp=c(0.0408,0.0401,0.0370,0.0350), t=c(1:4),mv=100,n=4)
#'pricingWithSptSeq (cpns=c(3,3,3,3), sp=c(0.0039,0.0140,0.0250,0.0360), t=c(1:4),mv=100,n=4)
#'pricingWithSptSeq (c(3,3,3,3), c(0.0408,0.0401,0.0370,0.0350), c(1:4),100,4)
#'pricingWithSptSeq (c(3,3,3,3), c(0.0039,0.0140,0.0250,0.0360), c(1:4),100,4)
#'pricingWithSptSeq (cpns=c(rep(3,4)), sp=c(0.0408,0.0401,0.0370,0.0350), t=c(1:4),mv=100,n=4)
#'pricingWithSptSeq (cpns=c(rep(3,4)), sp=c(0.0039,0.0140,0.0250,0.0360), t=c(1:4),mv=100,n=4)
#'pricingWithSptSeq (c(rep(3,4)), c(0.0408,0.0401,0.0370,0.0350), c(1:4),100,4)
#'pricingWithSptSeq (c(rep(3,4)), c(0.0039,0.0140,0.0250,0.0360), c(1:4),100,4)
#'@export
pricingWithSptSeq <- function (cpns, sp, t,mv,n ){
  {
  all = list(cpns, t)
  if (!is.vector(cpns) | !is.numeric(cpns))
    stop("Note: coupon payments must be a numeric Values.")
  if (!is.vector(t) | !is.numeric(t))
    stop("Note:t represents number of coupon payments must be 1 to number of years as a numeric vector.")
    if (any(t <= 0))
    stop("Note:occurance of coupons cannot have negative values in times.")
  if (length(cpns) == 0 | length(t) == 0)
    stop("Note:Please enter Complete Information.")
  if (length(cpns) != length(t))
    stop("Note: Number of Coupon Payments not equal to times.")
  }
  pv_coupon <- sum(cpns/(1 + sp)^t)
  pv_mv <- mv/(1 + (sp[n]))^n
  bond_price = pv_coupon + pv_mv
  (bond_price = round(bond_price, digits=3))
  }


