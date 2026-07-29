#'Calculates value of a share that is held for a single period (that is one year) using the Dividend Discount Model(DDM).
#'@description
#'From the perspective of a shareholder who buys and holds a share of stock, the cash flows he or she will obtain are the dividends paid on it and the market price of the share when he or she sells it. The future selling price should in turn reflect expectations about dividends subsequent to the sale. Assuming an investor wishes to buy a share of stock and hold it for one year, the value of that share of stock today is the present value of the expected dividend to be received on the stock plus the present value of the expected selling price at the end of one year (Jerald E. Pinto, 2020).
#'@details
#'According to information provided by Jerald E. Pinto (2020), the method \code{shareValueUsingDDM1yr} is developed to compute DDM value of share with a single holding period (that is one year) for the values passed to its four arguments. Here, \code{dividend1yr} is the expected dividend per share for Year 1, assumed to be paid at the end of the one year,  \code{expSharePriceIn1yr} is  the expected price per share at the end of one year, \code{n} is 1 representing that share is held for one year, and  \code{r} is the discount rate.
#'@param dividend1yr A number.
#'@param expSharePriceIn1yr A number.
#'@param n A number.
#'@param r A number.
#'@return Input values to four arguments  \code{dividend1yr}, \code{expSharePriceIn1yr}, \code{n} and  \code{r}.
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'Pinto, J. E. (2020). Equity Asset Valuation (4th ed.). Wiley Professional Development (P&T). https://bookshelf.vitalsource.com/books/9781119628194
#'@examples
#'shareValueUsingDDM1yr(dividend1yr=0.20,expSharePriceIn1yr=50,n=1, r=0.08)
#'shareValueUsingDDM1yr(dividend1yr=1.10,expSharePriceIn1yr=53.55,n=1, r=0.09)
#'@export
shareValueUsingDDM1yr <-function (dividend1yr,expSharePriceIn1yr,n, r)
{
  pv_dividend1yr <- dividend1yr/(1 + r)^n
  pv_expSharePriceIn1yr <- expSharePriceIn1yr/(1 + r)^n
  shareValue <-pv_dividend1yr + pv_expSharePriceIn1yr
  (shareValue = round(shareValue, digits=2))
}
