#'Calculates the Approximate Modified Duration.
#'@details
#'According to information provided by Adams and Smith (2019), the method \code{approxModifDuration()} is developed to calculate the Approximate Modified Duration.
#'@param pvBase A number.
#'@param pvPlus A number.
#'@param pvMinus A number.
#'@param percentChangeYtm A number.
#'@return Input values to four arguments  \code{pvBase},\code{pvPlus} ,\code{pvMinus}  and \code{percentChangeYtm}.
#'@importFrom Rdpack reprompt
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'Adams,J.F. & Smith,D.J.(2019). Understanding Fixed‑Income Risk and Return. In CFA Program Curriculum 2020 Level I Volumes 1-6. (Vol. 5, pp. 237-299). Wiley Professional Development (P&T). ISBN 9781119593577, <https://bookshelf.vitalsource.com/books/9781119593577>
#'@examples
#'approxModifDuration(pvBase=99.956780,pvPlus=100.631781,pvMinus=101.250227,percentChangeYtm=0.0005)
#'approxModifDuration (pvBase=42.223649,pvPlus=42.100694,pvMinus=42.346969,percentChangeYtm=0.0001)
#'@export
approxModifDuration <- function(pvBase, pvPlus,pvMinus,percentChangeYtm){
  (approxModifDuration = (pvMinus-pvPlus)/(2*percentChangeYtm*pvBase))
}


