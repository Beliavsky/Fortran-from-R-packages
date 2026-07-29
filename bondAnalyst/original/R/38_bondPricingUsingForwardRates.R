#'Calculates Bond Price using the Forward Rate Input.
#'@details
#'According to information provided by Adams and Smith (2019), the method \code{frPricing()} is developed to compute the Bond Price using the Forward Rate Input. If 0y1y is 1.88 percent, 1y1y is 2.77 percent and 2y1y is 3.54 percent then the bond price works out to be 102.965 dollars. Here, \code{cpns} is vector of coupon Payments for three years,\code{fri} is a vector of that contains three values that serve as forward rate input factor,\code{mv}is maturity Value, and \code{n} is number of years that is 3 years.
#'@param cpns A vector.
#'@param fri A vector.
#'@param mv A number.
#'@param n A number.
#'@return Input values to four arguments  \code{cpns} ,\code{fri}, \code{mv} and \code{n}.
#'@importFrom Rdpack reprompt
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'Adams,J.F. & Smith,D.J.(2019). Introduction to fixed-income valuation. In CFA Program Curriculum 2020 Level I Volumes 1-6. (Vol. 5, pp. 107-151). Wiley Professional Development (P&T). ISBN 9781119593577, <https://bookshelf.vitalsource.com/books/9781119593577>
#'@examples
#'frPricing(cpns=c(3.75,3.75,3.75),fri=c((1.0188),(1.0188*1.0277),(1.0188*1.0277*1.0354)),mv=100,n=3)
#'@export
frPricing<-function (cpns, fri,mv,n)
{
  all = list(cpns,fri,mv )
  if (!is.vector(cpns) | !is.numeric(cpns))
    stop("coupon payments must be a numeric vector.")

  if (length(cpns) == 0 | length(fri) == 0)
    stop("Insufficient information.")
  if (length(cpns) != length(fri))
    stop("Number of Coupon Payments not equal to Implied Forward Rates.")
  pv_coupon <- sum(cpns/fri)
  pv_mv <- mv/(fri[n])
  bond_price <- pv_coupon + pv_mv
  ( bond_price= round(bond_price, digits=3))
}

