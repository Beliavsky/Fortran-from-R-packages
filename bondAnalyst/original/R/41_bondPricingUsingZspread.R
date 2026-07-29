#'Calculates Bond Price using the given value of a Z-Spread and spot rates taken from the spots curve.
#'@details
#'According to information provided by Adams and Smith (2019), the method \code{pricingWithZspread()} is developed to compute Bond Price using the given value of a Z-Spread and spot rates taken from the spots curve. Here, \code{cpns} is vector of Coupon Payments,  \code{spots} is a vector of spot rates taken from the spots curve, \code{t} is a vector of number of years ranging from 1 to any specified number of years under consideration, \code{mv} is maturity value of the bond, \code{n} is number of years for which spots are available, and \code{zSprd} is given value of a Z-spread.
#'@param cpns A vector.
#'@param spots A vector.
#'@param t A vector.
#'@param mv A number.
#'@param n A number.
#'@param zSprd A number.
#'@return Input values to six arguments \code{cpns},  \code{spots} ,\code{t},\code{mv}, \code{n} and \code{zSprd}.
#'@importFrom Rdpack reprompt
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'Adams,J.F. & Smith,D.J.(2019). Introduction to fixed-income valuation. In CFA Program Curriculum 2020 Level I Volumes 1-6. (Vol. 5, pp. 107-151). Wiley Professional Development (P&T). ISBN 9781119593577, <https://bookshelf.vitalsource.com/books/9781119593577>
#'@examples
#'pricingWithZspread(cpns=c(6,6),spots=c(0.0210,0.03635),t=c(1,2),mv=100,n=2,zSprd=0.023422)
#'pricingWithZspread(cpns=c(5,5,5),spots=c(0.0486,0.0495,0.0565),t=c(1,2,3),mv=100,n=3,zSprd=0.0234)
#'pricingWithZspread(cpns=c(rep(5,3)),spots=c(0.0486,0.0495,0.0565),t=c(1:3),mv=100,n=3,zSprd=0.0234)
#'pricingWithZspread(c(rep(5,3)),c(0.0486,0.0495,0.0565), c(1:3),100,3, zSprd=0.0234)
#'@export
pricingWithZspread<-function (cpns, spots, t,mv,n, zSprd )
{
  all = list(cpns,spots, t,mv,n )
  if (!is.vector(cpns) | !is.numeric(cpns))
    stop("coupon payments must be a numeric vector.")
  if (!is.vector(t) | !is.numeric(t))
    stop("t represents number of coupon payments must be 1 to number of years as a numeric vector.")
  if (any(t <= 0))
    stop("occurance of coupons cannot have negative values in times.")
  if (length(cpns) == 0 | length(t) == 0)
    stop("Insufficient information.")
  if (length(cpns) != length(t))
    stop("Number of Coupon Payments not equal to times.")
  pv_coupon <- sum(cpns/(1 + zSprd+  spots)^t)
  pv_mv <- mv/(1 + zSprd+ (spots[n]))^n
  bond_price <- pv_coupon + pv_mv
  ( bond_price= round(bond_price, digits=3))
}


