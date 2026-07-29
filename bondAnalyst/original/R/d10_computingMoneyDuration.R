#'Calculates Money Duration of a Bond.
#'@details
#'According to information provided by Adams and Smith (2019), the method \code{moneyDuration()} is developed to calculate Money Duration of a Bond.
#'@param macDuration A number.
#'@param ytm A number.
#'@param pvFullBondPrice A number.
#'@return Input values to three arguments  \code{macDuration},\code{ytm} and \code{pvFullBondPrice}.
#'@importFrom Rdpack reprompt
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'Adams,J.F. & Smith,D.J.(2019). Understanding Fixed‑Income Risk and Return. In CFA Program Curriculum 2020 Level I Volumes 1-6. (Vol. 5, pp. 237-299). Wiley Professional Development (P&T). ISBN 9781119593577, <https://bookshelf.vitalsource.com/books/9781119593577>
#'@examples
#'moneyDuration(macDuration=7.0029,ytm=0.104,pvFullBondPrice=85.5031)
#'@export
moneyDuration<- function(macDuration,ytm,pvFullBondPrice){
   moneyDuration = macDuration/(1+ytm) *pvFullBondPrice
  (moneyDuration = round(moneyDuration, digits=4))
}

