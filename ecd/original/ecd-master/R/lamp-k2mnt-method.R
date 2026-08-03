#' Conversion between cumulants and moments
#' 
#' Implements conversion between the first four cumulants and moments
#'
#' @param k numeric, first four cumulants.
#' @param m numeric, first four moments.
#' 
#' @return numeric
#'
#' @keywords moments
#'
#' @author Stephen H-T. Lihn
#'
#' @export k2mnt
#' @export mnt2k
#' 
### <======================================================================>
k2mnt <- function(k) {
    m1 <- k[1]
    m2 <- k[2]+k[1]^2
    m3 <- k[3]+3*k[2]*k[1]+k[1]^3
    m4 <- k[4]+4*k[3]*k[1]+3*k[2]^2+6*k[2]*k[1]^2+k[1]^4
    c(m1,m2,m3,m4)
}
### <---------------------------------------------------------------------->
#' @rdname k2mnt
mnt2k <- function(m) {
    k1 <- m[1]
    k2 <- m[2]-m[1]^2
    k3 <- m[3]-3*m[2]*m[1]+2*m[1]^3
    k4 <- m[4]-4*m[3]*m[1]-3*m[2]^2+12*m[2]*m[1]^2-6*m[1]^4
    c(k1,k2,k3,k4)
}
### <---------------------------------------------------------------------->
