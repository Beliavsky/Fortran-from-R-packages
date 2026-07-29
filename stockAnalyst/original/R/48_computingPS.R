#'Calculates Price to Sales (PS) Multiple as trailing PS or GGM based PS.
#'@description
#'Certain types of privately held companies, including investment management companies and many types of companies in partnership form, have long been valued by a multiple of annual revenues. In recent decades, the ratio of price to sales has become well known as a valuation indicator for the equity of publicly traded companies as well. Analyst use PS Multiple as Sales are generally less subject to distortion or manipulation than are other fundamentals, such as EPS or book value. For example, through discretionary accounting decisions about expenses, company managers can distort EPS as a reflection of economic performance. In contrast, total sales, as the top line in the income statement, is prior to any expenses.Although the determination of sales is more straightforward than the determination of earnings, the analyst should evaluate a company’s revenue recognition practices, in particular, those tending to speed up the recognition of revenues—before relying on the P/S multiple. Trailing PS is calculated as price per share divided by annual net sales per share (net sales is total sales minus returns and customer discounts).Like other multiples, PS can be based on forecasted fundamentals like growth based on Gordon growth model (Jerald E. Pinto, 2020).
#'@details
#'According to information provided in Jerald E. Pinto (2020), the method \code{computingPS} is developed for computing Price to Sales (PS) Multiple as trailing PS or GGM based PS for the values passed to its seven arguments. Here, \code{PS} is character string, either trialing or GGM , \code{currentShPrice} is current Share Price , \code{payout} is payout ratio,\code{EPS0} is current earnings per share, \code{S0} is sales per share, \code{g} is earnings growth rate, and  \code{r} is required rate of return on equity.
#'@param  PS  character vector.
#'@param  currentShPrice number.
#'@param  payout number.
#'@param  EPS0 number.
#'@param  S0 number.
#'@param  g number.
#'@param  r number.
#'@return Input values to seven arguments  \code{PS}, \code{currentShPrice}, \code{payout}, \code{EPS0},\code{S0}, \code{g}, and \code{r}.
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'Pinto, J. E. (2020). Equity Asset Valuation (4th ed.). Wiley Professional Development (P&T). https://bookshelf.vitalsource.com/books/9781119628194
#'@examples
#'computingPS("trailing", currentShPrice=20,payout=0.35,EPS0=0.9,S0=10,g=0.07,r=0.09)
#'computingPS("GGM",  currentShPrice=20,payout=0.35,EPS0=0.9,S0=10,g=0.07,r=0.09)
#'@export
computingPS<-function (PS=c("trialing","GGM"),currentShPrice,payout,EPS0,S0,g,r){
  if (PS != "GGM" )
  { PS_ratio<- currentShPrice/S0
  (PS_ratio = round(PS_ratio, digits=2))}
  else
  { PS_ratio<- ((EPS0/S0)*payout*(1+g))/(r-g)
  (PS_ratio = round(PS_ratio, digits=1))}
}
