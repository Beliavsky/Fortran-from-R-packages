#'Calculates Present Value or the Price of the Bond paying semi-annual Coupons.
#'@details
#'The method  \code{pricingSaCpnBond()} is developed to compute the Price of Bond making semi-annual Coupon Payments. So,  \code{pricingSaCpnBond()} gives the Price of Bond making semi-annual Coupon Payments for values passed to its five arguments. Here, \code{saCoupons} represents the dollar values of all the semi-annual coupon payments till maturity, \code{times} is a vector of number of years ranging from 1 to any specified number of semi-annual periods till maturity, that is from 1 to times (n * 2), \code{maturityVal} is Maturity Value, \code{n} is number of years till maturity, and \code{r} is annual Market Discount Rate or Required Rate of return. The output is rounded off to two decimal places. The given examples show various ways in which the arguments can be passed to \code{pricingSaCpnBond()}.
#'@param saCoupons A vector.
#'@param times A vector.
#'@param maturityVal A number.
#'@param n A number.
#'@param r A number.
#'@return Input values to five arguments  \code{saCoupons} , \code{times}, \code{maturityVal},\code{n} and \code{r}.
#'@importFrom Rdpack reprompt
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'Adams,J.F. & Smith,D.J.(2019). Introduction to fixed-income valuation. In CFA Program Curriculum 2020 Level I Volumes 1-6. (Vol. 5, pp. 107-151). Wiley Professional Development (P&T). ISBN 9781119593577, <https://bookshelf.vitalsource.com/books/9781119593577>
#'@examples
#'pricingSaCpnBond(saCoupons=c(4,4,4,4,4,4), times=c(1,2,3,4,5,6) ,maturityVal=100,n=3,r=0.07)
#'pricingSaCpnBond(saCoupons=c(4,4,4,4,4,4), times=c(1:6) ,maturityVal=100,n=6,r=0.06)
#'pricingSaCpnBond(saCoupons=c(rep(4,6)), times=c(1:6) ,maturityVal=100,n=6,r=0.06)
#'pricingSaCpnBond(c(rep(4,6)), c(1:6) ,100,6,0.06)
#'@export
pricingSaCpnBond <-function (saCoupons, times,maturityVal,n, r)
{
  all = list(saCoupons, times, r)
  if (!is.vector(saCoupons) | !is.numeric(saCoupons))
    stop("Note: Semi-Annual coupon payments Amount must be a numeric vector.")
  if (!is.vector(times) | !is.numeric(times))
    stop("Note: times(1 : n*2) represents number of semi-Annual coupon payments as a numeric vector.")
  if (any(times <= 0))
    stop("Note: Occurance of semi-Annual oupon cannot have negative values in times.")
  if (r < 0)
    stop("Note: r is Annual market Discount rate and cannot be negative.")
  if (length(saCoupons) == 0 | length(times) == 0)
    stop("Note: Insufficient information.")
  if (length(saCoupons) != length(times))
    stop("Note: Number of Coupon Payments not equal to times(n*2).")
  pv_coupon <- sum(saCoupons/(1 + r/2)^times)
  pv_maturityVal <- maturityVal/(1 + r/2)^(n*2)
  bond_price <- pv_coupon + pv_maturityVal
  (bond_price = round(bond_price, digits=2))
}

