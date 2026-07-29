#'Calculates the share value from total Equity Value (based on single stage constant growth) that is divided by number of outstanding shares.
#'@details
#'According to information provided by Jerald E. Pinto (2020), the method \code{shareValConstantG} is developed to compute estimated value of the equity when FCFE is growing at a constant rate for the values passed to its three arguments. Here, \code{FCFE0} is given amount of future Free Cash Flow to the Equity in millions of dollars, \code{g} is constant rate of growth under single stage constant growth model, \code{WACC} is Weighted Average Cost of Capital, and \code{shares} is number of shares in millions.
#'@param FCFE0 A number.
#'@param g A number.
#'@param WACC A number.
#'@param shares A number.
#'@return Input values to four arguments  \code{FCFE0} \code{g}, \code{WACC}, and \code{shares}.
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'Pinto, J. E. (2020). Equity Asset Valuation (4th ed.). Wiley Professional Development (P&T). https://bookshelf.vitalsource.com/books/9781119628194
#'@examples
#'shareValConstantG(FCFE0=1.8,g=0.08,WACC=0.12,shares=1.5 )
#'shareValConstantG(FCFE0=700,g=0.05,WACC=0.102,shares=200)
#'@export
shareValConstantG<-function (FCFE0,g,WACC,shares)
{
  equityVal <- FCFE0*(1+g)/(WACC-g)
   shareValue <- equityVal/shares
  ( shareValue = round( shareValue, digits=2))
}
