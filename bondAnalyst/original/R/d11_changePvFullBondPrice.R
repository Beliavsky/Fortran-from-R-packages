#'Calculates estimated change in the Full Price of the Bond (in currency units) for a given Money Duration and a given change in the Yield-To-Maturity.
#'@details
#'According to information provided by Adams and Smith (2019), the method \code{changePvFullBondPrice()} is developed to compute estimated change in the Full Price of the Bond (in currency units) for a given Money Duration and a given Change in the Yield-To-Maturity.
#'@param moneyDuration A number.
#'@param changeYtm A number.
#'@return Input values to two arguments  \code{moneyDuration} and \code{changeYtm}.
#'@importFrom Rdpack reprompt
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'Adams,J.F. & Smith,D.J.(2019). Understanding Fixed‑Income Risk and Return. In CFA Program Curriculum 2020 Level I Volumes 1-6. (Vol. 5, pp. 237-299). Wiley Professional Development (P&T). ISBN 9781119593577, <https://bookshelf.vitalsource.com/books/9781119593577>
#'@examples
#'changePvFullBondPrice(moneyDuration=542.3638,changeYtm=0.01)
#'changePvFullBondPrice(moneyDuration=542.3638,changeYtm=-0.01)
#'@export
changePvFullBondPrice<- function(moneyDuration,changeYtm){
  change_PvFullPrice <- (-moneyDuration*changeYtm)
  (change_PvFullPrice = round(change_PvFullPrice, digits=4))
}

