## Mostly for ease of comparison of  our R code  with R's  C code :
##  (*not* exported)

DBL_MAX     <- .Machine$double.xmax
DBL_MIN     <- .Machine$double.xmin # 2.225e-308
DBL_EPSILON <- .Machine$double.eps
DBL_MANT_DIG<- .Machine$double.digits     # 53

ML_POSINF <-  Inf
ML_NEGINF <- -Inf
ML_NAN <- NaN

## from R's    <Rsrc>/src/nmath/dpq.h :
##
M_LN2 <- log(2) # 0.693147180559945309417232121458 # ln(2)
M_SQRT2 <- sqrt(2)
M_PI   <- pi   # 3.141592653589793238462643383280 # pi
M_1_PI <- 1/pi # 0.318309886183790671537767526745 # 1/pi
M_PI_2 <- pi/2 # 1.570796326794896619231321691640 # pi/2
M_2PI  <- 2*pi # 6.283185307179586476925286766559 # 2*pi
## (cannot use ldexp() , as dyn.load() has not happened at pkg build time)

