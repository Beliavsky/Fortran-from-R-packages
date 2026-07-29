#'Calculates PE Multiple of the companies with different abilities to  pass through the inflation to its customers through higher prices.
#'@description
#'While studying PE through cross-country comparisons the main differences in inflation rates and in the ability of companies to pass through inflation in their costs in the form of higher prices to their customers plays a vital role. For two companies with the same inflation pass-through ability, the company operating in the environment with higher inflation will have a lower justified PE; if the inflation rates are equal but pass-through rates differ, the justified PE should be lower for the company with the lower pass-through rate (Jerald E. Pinto, 2020).
#'@details
#'According to information obtained from Jerald E. Pinto (2020), the method \code{PEforPassThroughInflation} is developed for computing PE Multiple of the companies with different abilities to  pass through the inflation to customers for values passed to its three arguments.  Here, \code{realROR} is real Rate of Return, \code{I} is rate of Inflation, and \code{passThruRate} is percentage of inflation in costs that the company can pass through to its customers through higher prices.
#'@param  realROR number.
#'@param  I number.
#'@param  passThruRate number.
#'@return Input values to three arguments  \code{realROR}, \code{I}, and \code{passThruRate}.
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@examples
#'PEforPassThroughInflation(realROR=0.03,I=0.06,passThruRate=0.70)
#'PEforPassThroughInflation(realROR=0.03,I=0.06,passThruRate=0.90)
#'@export
PEforPassThroughInflation<-function (realROR,I,passThruRate){
 Disp_Val<- 1/(realROR+I*(1-passThruRate))
  (Disp_Val = round(Disp_Val, digits=1))
}
