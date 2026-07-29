#'Calculates the Approximated Macaulay duration using the Approximate Modified Duration and Yield-To-Maturity.
#'@details
#'According to information provided by Adams and Smith (2019), the method \code{approxMacDurationUsingApprModifDuration()} is developed to calculate the Approximated Macaulay duration using the Approximate Modified Duration and Yield-To-Maturity.
#'@param approxModifDuration A number.
#'@param periodicYtm A number.
#'@return Input values to two arguments  \code{approxModifDuration}  and \code{periodicYtm}.
#'@importFrom Rdpack reprompt
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'Adams,J.F. & Smith,D.J.(2019). Understanding Fixed‑Income Risk and Return. In CFA Program Curriculum 2020 Level I Volumes 1-6. (Vol. 5, pp. 237-299). Wiley Professional Development (P&T). ISBN 9781119593577, <https://bookshelf.vitalsource.com/books/9781119593577>
#'@examples
#'approxMacDurationUsingApprModifDuration(approxModifDuration=13.466,periodicYtm=0.0257)
#'@export
approxMacDurationUsingApprModifDuration<- function(approxModifDuration,periodicYtm){
  approxMacDuration <- approxModifDuration * (1+ periodicYtm)
  (approxMacDuration = round(approxMacDuration, digits=3))
}


