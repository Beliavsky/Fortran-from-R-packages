#'Calculates Price-to-Earnings Multiple by Yardeni Model that incorporates the impact of long-term expected growth rate of earnings on PE.
#'@description
#'The Long-term Earning Growth Model given by Yardeni in 2000 (as cited in Jerald E. Pinto, 2020) incorporates the expected growth rate in earnings, a variable that is missing in the Fed Model. This model is known as Yardeni Model and it incorporates the impact of long-term expected growth rate of earnings on PE and thereby overcomes the issue that was limitation of the US FED Model (Jerald E. Pinto, 2020).
#'@details
#'According to information provided by Jerald E. Pinto (2020), the method \code{impliedPEbyYardeniModel} is developed for computing Price to Earnings Multiple by Yardeni Model that that incorporates the expected growth rate of earnings for values passed to its four arguments.  Here, \code{CBY} is corporate bond yield, \code{b} is given coefficient of LTEG.The coefficient \code{b} measures the weight the market gives to five- year earnings projections, \code{LTEG} is Long Term Earning Growth.LTEG is taken as the consensus five- year earnings growth rate forecast for the market index and \code{residualVal} is residual value of the estimator that tends to zero.
#'@param  CBY number.
#'@param  b number.
#'@param  LTEG number.
#'@param  residualVal number.
#'@return Input values to four arguments  \code{CBY} \code{b}, \code{LTEG},and \code{residualVal}.
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@examples
#'impliedPEbyYardeniModel(CBY=0.06,b=0.2,LTEG=0.025,residualVal=0)
#'@export
impliedPEbyYardeniModel<-function (CBY,b,LTEG,residualVal){
  disp_val<- 1/(CBY-(b*LTEG)-residualVal)
  (disp_val = round(disp_val, digits=1))
}
