#'Calculates Price Value of a Basis Point (PVBP) for the Bond.
#'@details
#'According to information provided by Adams and Smith (2019), the method \code{computingBondPVBP()} is developed to compute Price Value of a Basis Point (PVBP) for the Bond.
#'@param pvPlus A number.
#'@param pvMinus A number.
#'@return Input values to two arguments  \code{pvPlus} and \code{pvMinus}.
#'@importFrom Rdpack reprompt
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'Adams,J.F. & Smith,D.J.(2019). Understanding Fixed‑Income Risk and Return. In CFA Program Curriculum 2020 Level I Volumes 1-6. (Vol. 5, pp. 237-299). Wiley Professional Development (P&T). ISBN 9781119593577, <https://bookshelf.vitalsource.com/books/9781119593577>
#'@examples
#'computingBondPVBP(pvPlus=100.594327,pvMinus=100.765123)
#'@export
computingBondPVBP<- function(pvPlus,pvMinus){
  pvbp = (pvMinus-pvPlus)/2
  (pvbp = round(pvbp, digits=4))
}


