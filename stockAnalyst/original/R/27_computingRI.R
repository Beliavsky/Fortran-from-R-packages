#'Calculates per share Residual Income using given values of Earnings Per Share (EPS) and beginning Book Values Per Share (bgnBVPS) for a specified number of years.
#'@details
#'According to information provided by Jerald E. Pinto (2020), the method \code{computingRI} is developed to compute value of share using Residual Income Model with given values of Earnings Per Share (EPS) and beginning Book Values Per Share (bgnBVPS) for a specified number of years for the values passed to its four arguments. Here, \code{bgnBVPS} is a vector of the beginning or current book value per share for a specified number of years,  \code{EPS} is a vector of the given values of Earnings Per Share for a specified number of years, and  \code{r} is the required rate of return on the stock. The \code{computingRI} computes Residual Incomes as EPS minus per share equity charge for specified number of years and then computes sum of discounted values of Residual Income that is added to current Book value per share to arrive at the share value.
#'@param bgnBVPS A number vector.
#'@param EPS A number vector.
#'@param r A number.
#'@return Input values to three arguments  \code{bgnBVPS} \code{EPS},and \code{r}.
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'Pinto, J. E. (2020). Equity Asset Valuation (4th ed.). Wiley Professional Development (P&T). https://bookshelf.vitalsource.com/books/9781119628194
#'@examples
#'computingRI(bgnBVPS=c(6,7,8.25),EPS=c(2,2.5,4),r=0.10 )
#'@export
computingRI<-function (bgnBVPS,EPS,r)
{
  RI <- (EPS-r*bgnBVPS)
    ( RI= round( RI, digits=4))
}


