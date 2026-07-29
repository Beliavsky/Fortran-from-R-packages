#'Calculate value of a share using the two-stage Dividend Discount Model (DDM).
#'@description
#'Two-stage DDM provides for a high growth rate for the initial period, followed by a sustainable and usually lower growth rate thereafter. The two-stage DDM is based on the multiple-period model. The two-stage model assumes that the first n dividends grow at an extraordinary short-term rate(gS) and after time n, the annual dividend growth rate changes to a normal long-term rate (gL). The two-stage DDM is useful because many scenarios exist in which a company can achieve a super-normal growth rate for a few years, after which time the growth rate falls to a more sustainable level. For example, a company may achieve super-normal growth through possession of a patent, first-mover advantage, or another factor that provides a temporary lead in a specific marketplace. Subsequently, earnings will most likely descend to a level that is more consistent with competition and growth in the overall economy. Accordingly, that is why in the two-stage model, extraordinary growth is often forecast for a few years and normal growth is forecast thereafter. A possible limitation of the two-stage model is that the transition between the initial abnormal growth period and the final steady-state growth period is abrupt (Jerald E. Pinto, 2020).
#'@details
#'According to information provided by Jerald E. Pinto (2020), the method \code{shareValUsingtwoStageHmodel} is developed to compute value of share using two stage H-Model for the values passed to its five arguments. Here, \code{divNot} is dollar value of the current dividend, \code{r} is required rate of return on equity, \code{n} is number of years of super-normal growth period, \code{gS} is initial short-term dividend growth rate, and \code{gL} is normal long-term dividend growth rate.
#'@param divNot A number.
#'@param r A number.
#'@param n A number.
#'@param gS A number.
#'@param gL A number.
#'@return Input values to five arguments  \code{divNot}  , \code{r}, \code{n}, \code{gS} and \code{gL}.
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'Pinto, J. E. (2020). Equity Asset Valuation (4th ed.). Wiley Professional Development (P&T). https://bookshelf.vitalsource.com/books/9781119628194
#'@examples
#'shareValUsingTwoStageDDM(divNot=0.14, r=0.097,n=10,gS=0.15,gL=0.08)
#'shareValUsingTwoStageDDM(divNot=0.40, r=0.071,n=10,gS=0.09,gL=0.05)
#'@export
shareValUsingTwoStageDDM <-function(divNot,r,n,gS,gL){
  shareVal<- (divNot*((1+gS)^n)*(1+gL))/(r-gL)
  (shareVal = round(shareVal, digits=4))
}
