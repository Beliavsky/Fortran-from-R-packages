#'Calculates value of a share using the given Residual Income.
#'@description
#'In the long term, companies that earn more than the cost of capital should sell for more than book value, and companies that earn less than the cost of capital should sell for less than book value. The residual income model of valuation analyzes the intrinsic value of equity as the sum of two components; first the current or the beginning book value of equity, and second, the present value of expected future residual income.
#'@details
#'According to information provided by Jerald E. Pinto (2020), the method \code{shareValueRI} is developed to compute value of share using Residual Income Model with given values of Residual Income for the values passed to its four arguments. Here, \code{bgnBVPS} is the beginning or current book value per share, \code{RI} is  a vector of given values of Residual Income for a specified number of years, \code{times} is a vector of number of years ranging from 1 to any specified number of years Residual Income Values are given, and \code{r} is the required rate of return on the stock.
#'@param bgnBVPS A number.
#'@param RI A vector.
#'@param r A number.
#'@param times A vector.
#'@return Input values to four arguments  \code{bgnBVPS} \code{RI}, \code{r}, , \code{times}.
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'Pinto, J. E. (2020). Equity Asset Valuation (4th ed.). Wiley Professional Development (P&T). https://bookshelf.vitalsource.com/books/9781119628194
#'@examples
#'shareValueRI(bgnBVPS=6,RI=c(1.40,1.80,3.175),r=0.10, times=c(1,2,3) )
#'@export
shareValueRI<-function (bgnBVPS,RI,r,times)
{
  pv_RI  <-sum(RI/(1 + r)^times)
  share_Value <- bgnBVPS + pv_RI
  ( share_Value = round( share_Value, digits=2))
}


