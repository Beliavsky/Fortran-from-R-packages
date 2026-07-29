#'Calculates share value using three-stage Free Cash Flow Model.
#'@description
#'Three-stage models are a straightforward extension of the two-stage models (Jerald E. Pinto, 2020). The example used in this method uses cash flow values for four years only. In the given example, stage one lasts for first two years and has cash flows of 2.8 million dollars each year with growth rate of 8.8 percent. The second stage, in the given example, lasts just for one year and has cash flows of 2.8 million dollars with growth rate of 7.4 percent, and the third stage, in the given example, has cash flows of 2.8 million dollars with low growth rate of 6.6 percent. However, it is more practical to you have values for, let us say 8 years, where first stage of high constant growth continues let us say for four years, followed by second stage of declining growth for three years, and then third stage of low constant growth thereafter.
#'@details
#'The version of a three-stage model used here assumes constant high growth in first stage, declining growth in second transitory stage and low constant growth in third stage.The method \code{shareValThreeStg} is developed to compute share value using three stage Free Cash Flow Model for the values passed to its five arguments. Here, \code{FCFE} is a vector of given amounts of future Free Cash Flow to the Equity (FCFE) in millions of dollars ,  \code{t} is vector of number of years ranging from 1 to any specified number of years used for computing the cumulative value of given Free Cash Flows , \code{G} is a vector of Growth rates in all the three stages,  \code{r} is required rate of return on equity (WACC can be used as r here), and  \code{s} is number of shares in millions, so a value of 0.5 means 500,000 outstanding shares.
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
#'shareValThreeStg(FCFE=c(2.8,2.8,2.8,2.8),t=c(1,2,3,4),G=c(0.088,0.088,0.074,0.066),r=0.1,s=0.5)
#'@export
shareValThreeStg<-function (FCFE,t,G,r,s)
{
  FCFE_val <- (FCFE*G)
  equityVal <- sum(FCFE_val/(1 + r)^t)
  shareValu <-equityVal/s
  (shareValu = round(shareValu, digits=2))
}
