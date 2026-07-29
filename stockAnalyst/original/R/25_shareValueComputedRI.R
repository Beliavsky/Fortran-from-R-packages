#'Calculates value of a share using given values of Earnings Per Share (EPS) and beginning Book Values Per Share (bgnBVPS) for a specified number of years.
#'@description
#'This valuation is sum of two components; first, the current or the beginning book value of equity, and second, the present value of expected future residual income that is computed as Earnings Per Share (EPS) minus the required rate of return on equity that is multiplied with beginning book value of per share. So, in this method RI is computed as discussed above and then this dollar value of Residual Income is discounted at the required rate of return on the equity and then beginning Book Value per Share (bgnBVPS) is added to arrive at the share value.
#'@details
#'According to information provided by Jerald E. Pinto (2020), the method \code{shareValueComputedRI} is developed to compute value of share using Residual Income Model with given values of Earnings Per Share (EPS) and beginning Book Values Per Share (bgnBVPS) for a specified number of years for the values passed to its four arguments. Here, \code{bgnBVPS} is a vector the beginning or current book value per share for a specified number of years, \code{EPS} is  a vector of given values of Earnings Per Share for a specified number of years, \code{times} is a vector of number of years ranging from 1 to any specified number of years Residual Income Values are to be computed, and \code{r} is the required rate of return on the stock. As an internal step \code{shareValueComputedRI} computes Residual Incomes as EPS minus per share equity charge for specified number of years and then computes sum of discounted values of Residual Income that is added to current Book value per share to arrive at the share value.
#'@param bgnBVPS A number.
#'@param EPS A vector.
#'@param r A number.
#'@param times A vector.
#'@return Input values to four arguments \code{bgnBVPS} \code{EPS}, \code{r}, and \code{times}.
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'Pinto, J. E. (2020). Equity Asset Valuation (4th ed.). Wiley Professional Development (P&T). https://bookshelf.vitalsource.com/books/9781119628194
#'@examples
#'shareValueComputedRI(bgnBVPS=c(6,7,8.25),EPS=c(2,2.5,4),r=0.10, times=c(1,2,3) )
#'@export
shareValueComputedRI<-function (bgnBVPS,EPS,r,times)
{
  RI <- (EPS-r*bgnBVPS)
  pv_RI  <- sum(RI/(1 + r)^times)
  share_Value <- bgnBVPS[1] + pv_RI
  ( share_Value= round( share_Value, digits=2))
}


