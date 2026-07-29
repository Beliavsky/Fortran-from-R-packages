#'Calculates Bond Price using given values of G-Spread and yield-to-maturity for the government benchmark bond.
#'@details
#'According to information provided by Adams and Smith (2019),  the method \code{pricingWithGspread} is developed to compute Bond Price using the given values of G-Spread and the yield-to-maturity for the government benchmark bond. Here, \code{coupons} is vector of Coupon Payments,\code{t} is a vector of number of years ranging from 1 to any specified number of years under consideration, \code{mv} is maturity value of the bond, \code{n} is number of years, \code{ytmBenchGovtBond} is the yield-to-maturity for the government benchmark bond, and \code{Gspread} is given value of G-spread.
#'@param coupons A vector.
#'@param t A vector.
#'@param mv A number.
#'@param n A number.
#'@param ytmBenchGovtBond A number.
#'@param Gspread A number.
#'@return Input values to six arguments  \code{coupons} ,\code{t}, \code{mv},\code{n}, \code{ytmBenchGovtBond}, and \code{Gspread}.
#'@importFrom Rdpack reprompt
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'Adams,J.F. & Smith,D.J.(2019). Introduction to fixed-income valuation. In CFA Program Curriculum 2020 Level I Volumes 1-6. (Vol. 5, pp. 107-151). Wiley Professional Development (P&T). ISBN 9781119593577, <https://bookshelf.vitalsource.com/books/9781119593577>
#'@examples
#'pricingWithGspread(coupons=c(6,6), t=c(1,2),mv=100,n=2,ytmBenchGovtBond=0.03605, Gspread=0.02327)
#'pricingWithGspread(coupons=c(5,5,5),t=c(1,2,3),mv=100,n=3,ytmBenchGovtBond=0.01913,Gspread=0.0285)
#'pricingWithGspread(coupons=c(60,60),t=c(1,2),mv=1000,n=2,ytmBenchGovtBond=0.03605,Gspread=0.02327)
#'@export
pricingWithGspread <-function (coupons, t,mv,n, ytmBenchGovtBond, Gspread){
  {
    all = list(coupons, ytmBenchGovtBond)
  if (!is.vector(coupons) | !is.numeric(coupons))
    stop("Note:coupon payments must be a numeric values.")
  if (any(t <= 0))
    stop("Note: occurance of coupons cannot have negative values in time.")
  if (ytmBenchGovtBond < 0)
    stop("Note: Market Discount rate (r) cannot be negative.")
  if (length(coupons) == 0 | length(t) == 0)
    stop("Note: Please enter Complete Information.")
  if (length(coupons) != length(t))
    stop("Note: Number of Coupon Payments not equal to time.")
    }
   pv_coupon <- sum((coupons/(1 + ytmBenchGovtBond+ Gspread)^t))
  pv_mv <- mv/(1 +ytmBenchGovtBond+Gspread )^n
  bond_price <- pv_coupon + pv_mv
  (bond_price= round(bond_price, digits=2))
}

