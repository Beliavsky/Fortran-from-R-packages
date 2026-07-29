#'Calculates Price to Book Value (PB) Multiple as trailing PB or GGM based PB.
#'@description
#'The ratio of market price per share to book value per share (PB), like PE, has along history of use in valuation practice as discussed by Graham and Dodd in 1934 (as cited in Jerald E. Pinto, 2020).In the measure of value in the PB denominator (book value per share) is a stock or level variable coming from the balance sheet. (Book refers to the fact that the measurement of value comes from accounting records or books, in contrast to market value.) Analysts use  PB because book value is a cumulative balance sheet amount, book value is generally positive even when EPS is zero or negative. An analyst can generally use PB when EPS is zero or negative, whereas P/E based on a zero or negative EPS is not meaningful.
#'@details
#'According to information provided by Jerald E. Pinto (2020), the method \code{computingPB} is developed for computing Price to Book Value (PB) Multiple as trailing PB or GGM based PB for the values passed to its six arguments. Here, \code{PB} is character string, either trailing or GGM , \code{currentShPrice} is current Share Price , \code{BV0} is initial Book Value,\code{ROE} is return on equity, \code{g} is sustainable growth rate under the Gordon growth model, and  \code{r} is required rate of return on equity.
#'@param  PB  character vector.
#'@param  BV0 number.
#'@param  currentShPrice number.
#'@param  ROE number.
#'@param  g number.
#'@param  r number.
#'@return Input values to six arguments  \code{PB}, \code{currentShPrice}, \code{ROE},\code{BV0}, \code{g}, and \code{r}.
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@examples
#'computingPB("trailing", currentShPrice=81.23,BV0=49.67,ROE=0.12,g=0.07,r=0.10)
#'computingPB("GGM",  currentShPrice=81.23,BV0=49.67,ROE=0.12,g=0.07,r=0.10)
#'@export
computingPB<-function (PB=c("trailing","GGM"),BV0,currentShPrice,ROE,g,r){
  if (PB != "GGM" )
  { PB_ratio<- currentShPrice/BV0
  (PB_ratio = round(PB_ratio, digits=2))}
  else
  { PB_ratio<- (ROE-g)/(r-g)
  (PB_ratio = round(PB_ratio, digits=2))}
}
