#'Calculates Modified Duration using the Macaulay Duration and Yield-To-Maturity.
#'@details
#'According to information provided by Adams and Smith (2019), the method \code{modifDurationUsingMacDuration()} is developed to calculate Modified Duration using the Macaulay Duration and Yield-To-Maturity.
#'@param macDuration A number.
#'@param ytm A number.
#'@return Input values to two arguments  \code{macDuration}  and \code{ytm}.
#'@importFrom Rdpack reprompt
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'Adams,J.F. & Smith,D.J.(2019). Understanding Fixed‑Income Risk and Return. In CFA Program Curriculum 2020 Level I Volumes 1-6. (Vol. 5, pp. 237-299). Wiley Professional Development (P&T). ISBN 9781119593577, <https://bookshelf.vitalsource.com/books/9781119593577>
#'@examples
#'modifDurationUsingMacDuration(macDuration=12.6213,ytm=0.06/2)
#'modifDurationUsingMacDuration(macDuration=7.0029,ytm=0.104)
#'@export
modifDurationUsingMacDuration <- function(macDuration,ytm){
  modif_duration= macDuration/(1+ytm)
  (modif_duration = round(modif_duration, digits=4))
}

