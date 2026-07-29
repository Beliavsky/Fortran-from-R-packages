#'Calculates Residual Income using given values of Earnings Per Share (EPS) and beginning Book Values Per Share(bgnBVPS) for a specified number of years.
#'@details
#' Residual Income is computed in three steps. Here, Step 1 is to compute \code{preTaxIncome} as (EBIT minus \code{rd} multiplied with debt). Step 2 is to get \code{netIncome} as (\code{preTaxIncome} minus (\code{t} multiplied with \code{preTaxIncome}) ), and finally step 3 is to obtain the Residual Income (RI) as \code{netIncome} minus (\code{r} times \code{equity}).
#'According to information provided by Jerald E. Pinto (2020), the method \code{computingAbsoluteRI} is developed to compute absolute value of Residual Income. Here, \code{EBIT} is a number vector that hold values of EBIT in millions of dollars, \code{debt} is a number vector that has dollar value of debt (expressed in millions of dollars), \code{equity} is a number vector that holds dollar value of equity (expressed in millions of dollars), \code{r} required rate of return on equity (expressed in decimal terms), \code{rd} is cost of debt (expressed in decimal terms), and \code{t} is rate of taxes. Output gives dollar value of Residual Income (expressed in millions of dollars).
#'@param EBIT A number vector.
#'@param debt A number vector.
#'@param equity A number vector.
#'@param r A number.
#'@param rd A number.
#'@param t A number.
#'@return Input values to six arguments  \code{bgnBVPS} \code{RI}, \code{r}, , \code{times}.
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'Pinto, J. E. (2020). Equity Asset Valuation (4th ed.). Wiley Professional Development (P&T). https://bookshelf.vitalsource.com/books/9781119628194
#'@examples
#'computingAbsRI(EBIT=c(0.5,1.5,2.25),debt=c(1,2.2,2.5),equity=c(1,2.2,2.5), r=0.12,rd=0.07,t=0.30)
#'@export
computingAbsRI<-function (EBIT,debt,equity,r,rd,t)
{
   preTaxIncome <- (EBIT-rd*debt)
   netIncome <- (preTaxIncome- (t*preTaxIncome) )
   RI <- netIncome - r*equity
  ( RI= round( RI, digits=4))
}

