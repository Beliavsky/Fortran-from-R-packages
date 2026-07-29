#'Calculates the Effective Duration statistic of a Callable Bond.
#'@details
#'According to information provided by Adams and Smith (2019), the method \code{effDurtnCallableBond()} is developed to calculate the Effective Duration statistic of a Callable Bond.
#'@param pvBase A number.
#'@param pvPlus A number.
#'@param pvMinus A number.
#'@param perChangeBenchYtm A number.
#'@return Input values to four arguments  \code{pvBase},\code{pvPlus} ,\code{pvMinus}  and \code{perChangeBenchYtm}.
#'@importFrom Rdpack reprompt
#'@author MaheshP Kumar, \email{maheshparamjitkumar@@gmail.com}
#'@references
#'Adams,J.F. & Smith,D.J.(2019). Understanding Fixed‑Income Risk and Return. In CFA Program Curriculum 2020 Level I Volumes 1-6. (Vol. 5, pp. 237-299). Wiley Professional Development (P&T). ISBN 9781119593577, <https://bookshelf.vitalsource.com/books/9781119593577>
#'@examples
#'effDurtnCallableBond(pvBase=101.060489,pvPlus=99.050120,pvMinus=102.890738,perChangeBenchYtm=0.0025)
#'@export
effDurtnCallableBond<- function(pvBase, pvPlus,pvMinus,perChangeBenchYtm){
  effDuration = (pvMinus-pvPlus)/(2*perChangeBenchYtm*pvBase)
  (effDuration = round(effDuration, digits=4))
}

