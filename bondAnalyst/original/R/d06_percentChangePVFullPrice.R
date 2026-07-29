#'Calculates the percentage change in Full Price of the Bond for a given a change in its Yield-To-Maturity and Modified Duration statistic.
#'@details
#'According to information provided by Adams and Smith (2019), the method \code{estimatedPercentChangePVFullPrice()} is developed to calculate the percentage change in Full Price of the Bond for a given a change in its Yield-To-Maturity and Modified Duration statistic.
#'@param annualModifDuration A number.
#'@param changeInAnnualYtm A number.
#'@return Input values to two arguments  \code{annualModifDuration}  and \code{changeInAnnualYtm}.
#'@importFrom Rdpack reprompt
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'Adams,J.F. & Smith,D.J.(2019). Understanding Fixed‑Income Risk and Return. In CFA Program Curriculum 2020 Level I Volumes 1-6. (Vol. 5, pp. 237-299). Wiley Professional Development (P&T). ISBN 9781119593577, <https://bookshelf.vitalsource.com/books/9781119593577>
#'@examples
#'estimatedPercentChangePVFullPrice(annualModifDuration=6.126829,changeInAnnualYtm=0.01)
#'estimatedPercentChangePVFullPrice(annualModifDuration=6.3432,changeInAnnualYtm=-0.01)
#'@export
estimatedPercentChangePVFullPrice <- function(annualModifDuration,changeInAnnualYtm){
  estmPerChange= annualModifDuration*changeInAnnualYtm
  (estmPerChange = round(estmPerChange, digits=6))
}


