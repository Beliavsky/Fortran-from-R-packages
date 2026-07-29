#'Calculates value of a share that is held for multiple holding periods (for n years) using the Dividend Discount Model (DDM).
#'@description
#'If an investor plans to hold a share of stock for two years, the value of the share is the present value of the expected dividend in Year 1, plus the present value of the expected dividend in Year 2, plus the present value of the expected selling price at the end of Year 2. For an n-period model, the value of a stock is the present value of the expected dividends for the n periods plus the present value of the expected price at the end of nth period (Jerald E. Pinto, 2020).
#'@details
#'According to information provided by Jerald E. Pinto (2020), the method \code{shareValueUsingDDMnYrs} is developed to compute DDM value of share with multiple holding periods (that is for n years) for the values passed to its five arguments. Here, \code{dividend} is the expected dividend per share for n years, assumed to be paid at the end each year,  \code{expSharePriceNyr} is  the expected price per share at the end of nth year, \code{times} is a vector of number of years ranging from 1 to any specified number of years for which share is being held ,\code{n}, for example, n with value of 2 represents that share is held for two years, and  \code{r} is the required rate of return on the stock.
#'@param dividend A vector.
#'@param expSharePriceNyr A number.
#'@param times A vector.
#'@param n A number.
#'@param r A number.
#'@return Input values to five arguments  \code{dividend} , \code{expSharePriceNyr}, \code{times}, \code{n} and  \code{r}.
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'Pinto, J. E. (2020). Equity Asset Valuation (4th ed.). Wiley Professional Development (P&T). https://bookshelf.vitalsource.com/books/9781119628194
#'@examples
#'shareValueUsingDDMnYrs(dividend=c(3,3.15),expSharePriceNyr=40,times=c(1,2),n=2,r=0.08)
#'shareValueUsingDDMnYrs(dividend=c(2,3),expSharePriceNyr=48,times=c(1,2),n=2,r=0.10)
#'shareValueUsingDDMnYrs(dividend=c(2,2.10,2.20),expSharePriceNyr=20,times=c(1,2,3),n=3,r=0.10)
#'@export
shareValueUsingDDMnYrs <-function(dividend,expSharePriceNyr,times,n, r){
  pv_dividend  <-sum(dividend/(1 + r)^times)
  pv_expSharePrice <- expSharePriceNyr/(1 + r)^n
  shareValue <-pv_dividend + pv_expSharePrice
  (shareValue = round(shareValue, digits=2))
}
