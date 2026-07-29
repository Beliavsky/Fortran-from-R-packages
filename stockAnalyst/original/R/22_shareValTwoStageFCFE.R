#'Calculates share value using two-stage Free Cash Flow Model.
#'@details
#'The version of a two-stage model used here assumes constant high growth in first stage and low rate of constant growth in the second stage. According to information provided by Jerald E. Pinto (2020), the method, \code{shareValTwoStage} is developed to compute the share value using two stage Free Cash Flow Model for the values passed to its five arguments. Here, \code{FCFE} is a vector of given amounts of future Free Cash Flow to the Equity (FCFE) in millions of dollars. The example given here uses values for four years only. However, it is more practical to you have values for say 7 years where first stage of high constant growth continues for let us say four years, followed by second stage of low constant growth of three years. In this case, \code{t} is vector of number of years ranging from 1 to any specified number of years used for computing the cumulative value of given Free Cash Flows and \code{G} is a vector of Growth rates in two stages. Here, high growth of 20 percent is in stage one that continues for three years and the second stage of low growth at 6 percent and after that  \code{r} is required rate of return on equity (WACC can be used as \code{r} here), and \code{s} is number of shares in millions so a value of 0.5 means 500,000 outstanding shares.
#'@param FCFE A vector.
#'@param t A vector.
#'@param G A vector.
#'@param r A number.
#'@param s A number.
#'@return Input values to five arguments  \code{FCFE}, \code{t}, \code{G}  \code{r} and  \code{s}.
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'Pinto, J. E. (2020). Equity Asset Valuation (4th ed.). Wiley Professional Development (P&T). https://bookshelf.vitalsource.com/books/9781119628194
#'@examples
#'shareValTwoStage(FCFE=c(1.8,1.8,1.8,1.8),t=c(1,2,3,4),G=c(0.20,0.20,0.20,0.06),r=0.124,s=0.5)
#'@export
shareValTwoStage<-function (FCFE,t,G,r,s)
{
  FCFE_val <- (FCFE*G)
  equityVal <- sum(FCFE_val/(1 + r)^t)
  shareValu <-equityVal/s
  (shareValu = round(shareValu, digits=2))
}
