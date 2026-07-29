#'Calculates value of a share based on single-stage (constant-growth) Residual Income model.
#'@description
#'The single-stage (constant-growth) residual income model assumes that a company has a constant return on equity and constant earnings growth rate through time.
#'@details
#'According to information provided by Jerald E. Pinto (2020), the method \code{singleStageR} is developed to compute value of a share based on single-stage (constant-growth) residual income model for the values passed to its four arguments. Here, \code{ROErate} is rate of Return on Equity, \code{g} is constant rate of growth under single stage constant growth model, \code{bgnBVPS} is beginning Book Value per Share, \code{r} is required rate of return on equity.
#'@param ROErate A number.
#'@param bgnBVPS A number.
#'@param r A number.
#'@param g A number.
#'@return Input values to four arguments  \code{bgnBVPS} \code{RI}, \code{r},and \code{g}.
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'Pinto, J. E. (2020). Equity Asset Valuation (4th ed.). Wiley Professional Development (P&T). https://bookshelf.vitalsource.com/books/9781119628194
#'@examples
#'singleStageR(ROErate=0.16, bgnBVPS=18.81,r=0.11,g=0.08)
#'@export
singleStageR<-function (ROErate,bgnBVPS,r,g){
  share_Value <- bgnBVPS + ((ROErate-r)/(r-g))* bgnBVPS
  ( share_Value= round( share_Value, digits=2))
}
