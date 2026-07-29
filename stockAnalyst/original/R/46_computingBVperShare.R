#'Calculates the Book Value (BV) per share.
#'@description
#'To compute book value per share, we need to refer to the balance sheet, which has a shareholders (or stockholders) equity section. The computation of book value is done through the following formula: Shareholders equity minus Total value of equity claims that are senior to common stock is equal to Common shareholders equity. After this, Common shareholders equity is divided by the number of common shares outstanding to get the Book value per share. Possible claims senior to the claims of common stock, which would be subtracted from shareholders’ equity, include the value of preferred stock and the dividends in arrears on preferred stock (Jerald E. Pinto, 2020).
#'@details
#'According to information provided in Jerald E. Pinto (2020), the method \code{computingBVperShare} is developed for computing the Book Value (BV) per share for the values passed to its three arguments. Here, \code{totalEquity} is  total market value of Common Equity, \code{prefStockMV} is market value of Preference Stock, and \code{outstdCommShares} is number of common stock shares that are outstanding.
#'@param  totalEquity number.
#'@param  prefStockMV number.
#'@param  outstdCommShares number.
#'@return Input values to three arguments  \code{totalEquity}, \code{prefStockMV}, and \code{outstdCommShares}.
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@examples
#'computingBVperShare(totalEquity=49000,prefStockMV=3396,outstdCommShares=918.2)
#'@export
computingBVperShare<-function (totalEquity, prefStockMV,outstdCommShares){
  Disp_Val<- (totalEquity-prefStockMV)/outstdCommShares
  (Disp_Val = round(Disp_Val, digits=2))
}
