#'Calculates value of a share using Feltham and Ohlson Model.
#'@description
#'The residual income model used here has its origins largely in the academic work of Feltham and Ohlson (as given by Feltham and Ohlson 1995, as cited in Jerald E. Pinto, 2020).
#'@details
#'According to information provided by Jerald E. Pinto (2020), the method \code{shareValueROE} is developed to compute value of share using Residual Income Model with given values of ROE and beginning Book Values Per Share(bgnBVPS) for a specified number of years for the values passed to its four arguments. Here, \code{bgnBVPS} is a vector of the beginning or current book value per share for a specified number of years, \code{ROE} is  a vector of given values of Return on Equity for a specified number of years, \code{r} is the required rate of return on the stock, and \code{times} is a vector of number of years ranging from 1 to any specified number of years Residual Income Values are to be computed. The \code{shareValueROE} computes Residual Incomes as EPS minus per share equity charge for specified number of years and then computes sum of discounted values of Residual Income that is added to current Book value per share to arrive at the share value.
#'@param bgnBVPS A number.
#'@param ROE A number vector.
#'@param r A number.
#'@param times A number vector.
#'@return Input values to four arguments  \code{bgnBVPS} \code{ROE}, \code{r},and \code{times}.
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'Pinto, J. E. (2020). Equity Asset Valuation (4th ed.). Wiley Professional Development (P&T). https://bookshelf.vitalsource.com/books/9781119628194
#'@examples
#'shareValueROE(ROE=c(0.3333,0.3571,0.4848), bgnBVPS=c(6,7,8.25),r=0.10,times=c(1,2,3))
#'@export
shareValueROE<-function (ROE,bgnBVPS,r,times){
  RI <- ((ROE-r)*bgnBVPS)
  pv_RI  <- sum(RI/(1 + r)^times)
  share_Value <- bgnBVPS[1] + pv_RI
  ( share_Value= round( share_Value, digits=2))
}
