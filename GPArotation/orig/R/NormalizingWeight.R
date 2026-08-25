# Kaiser normalization

#NormalizingWeight <- function(A, normalize=FALSE){
# if ("function" == mode(normalize)) normalize <- normalize(A)
# if (is.logical(normalize)){
#    if (normalize) normalize <- sqrt(rowSums(A^2))
#    else return(array(1, dim(A)))
#    }
# if (is.vector(normalize)) 
#    {if(nrow(A) != length(normalize))
#        stop("normalize length wrong in NormalizingWeight")
#     return(array(normalize, dim(A)))
#    }
# stop("normalize argument not recognized in NormalizingWeight")
#}
# 

#
# Version below submitted by Kim-Laura Speck, Uni Kassel, 25 October 2023
#
# avoid NaNs in matrix A by adding machine precision values to zeros


NormalizingWeight <- function(A, normalize = FALSE) {
    if (is.function(normalize)) 
        normalize <- normalize(A)
    
    if (is.character(normalize)) {
        normalize <- switch(normalize, 
                            Kaiser = pmax(sqrt(rowSums(A^2)), .Machine$double.eps),
                            
                            CM = {
                                # 1. Row communalities (h_i)
                                h <- pmax(sqrt(rowSums(A^2)), .Machine$double.eps)
                                
                                # 2. First unrotated factor projection (bounded to [-1, 1] for acos safety)
                                fpls <- pmin(pmax(A[, 1] / h, -1), 1)
                                
                                # 3. Angle calculations
                                m <- ncol(A)
                                acosi <- acos(m^(-1/2))
                                alpha_i <- acos(abs(fpls))
                                
                                # 4. Piecewise denominator
                                in_region1 <- abs(fpls) < (m^(-1/2))
                                dem <- acosi - ifelse(in_region1, pi / 2, 0)
                                num <- acosi - alpha_i
                                
                                # 5. Cureton-Mulaik weights (w_i)
                                wghts <- cos((num / dem) * (pi / 2))^2 + 0.001
                                
                                # 6. Return h / wghts as a 1D vector so A / (h / wghts) = (wghts / h) * A
                                h / wghts
                            },
                            stop("normalize string not recognized in NormalizingWeight", call. = FALSE)
        )
    }
    
    if (is.logical(normalize)) {
        if (normalize) 
            normalize <- pmax(sqrt(rowSums(A^2)), .Machine$double.eps)
        else 
            return(array(1, dim(A)))
    }
    
    if (is.vector(normalize)) {
        if (nrow(A) != length(normalize)) 
            stop("normalize length wrong in NormalizingWeight", call. = FALSE)
        return(matrix(normalize, nrow(A), ncol(A)))
    }
    
    stop("normalize argument not recognized in NormalizingWeight", call. = FALSE)
}

#function(A, normalize = FALSE) {

#  if (is.function(normalize)) normalize <- normalize(A)

  # String shortcuts --- convert to weight vector
#  if (is.character(normalize)) {
 #   normalize <- switch(normalize,
#      Kaiser = pmax(sqrt(rowSums(A^2)), .Machine$double.eps),
#      CM     = {
#        Dk    <- diag(1 / pmax(sqrt(diag(A %*% t(A))),
#                               .Machine$double.eps)) %*% A
#        wghts <- rep(0, nrow(A))
#        fpls  <- Dk[, 1]
#        acosi <- acos(ncol(A)^(-1/2))
#        for (i in 1:nrow(A)) {
#          num      <- acosi - acos(abs(fpls[i]))
 #         dem      <- acosi - (function(a, m)
#                        ifelse(abs(a) < (m^(-1/2)), pi/2, 0))(fpls[i], ncol(A))
#          wghts[i] <- cos(num / dem * pi/2)^2 + 0.001
#        }
#        wghts / pmax(sqrt(diag(A %*% t(A))), .Machine$double.eps)
#      },
#      stop("normalize string not recognized in NormalizingWeight", call. = FALSE)
#    )
#  }
#
  # Logical shortcut --- TRUE maps to Kaiser
 # if (is.logical(normalize)) {
 #   if (normalize)
 #     normalize <- pmax(sqrt(rowSums(A^2)), .Machine$double.eps)
 #   else
 #     return(array(1, dim(A)))
 # }

  # Weight vector --- apply to all columns
  #if (is.vector(normalize)) {
  #  if (nrow(A) != length(normalize))
  #    stop("normalize length wrong in NormalizingWeight", call. = FALSE)
  #  return(matrix(normalize, nrow(A), ncol(A)))
  #}
#
#  stop("normalize argument not recognized in NormalizingWeight", call. = FALSE)
#}
  