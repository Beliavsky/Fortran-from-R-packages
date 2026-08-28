# tests using normalization
# All tests below use Kaiser normalization
# A few other tests also use normalization when comparing varimax and Varimax

# Following examples are from SPSS
# See https://psych.unl.edu/psycrs/statpage/pc_rot.pdf

require("stats")  
require("GPArotation")  

fuzz <- 1e-3  #less strict; differences in 4rd decimal compared to SPSS
all.ok <- TRUE  

# unrotated matrix
L <- matrix(c(.758, .413, 1.164E-03, .693, .489, -.199, .362, .656, -.204, 
.826, 6.589E-02, .235, .540, -.510, .441, .654, -.335, .507,
-.349, .539, .669, -.580, .450, .551), byrow=T, ncol=3)

# quartimax, Kaiser normalization
# uses the print command to get the right order of factors
v <- print(quartimax(L, normalize = TRUE, eps = 1e-6))$loadings 

tst <- matrix(c(.814, .285, -4.99E-02, .856, 8.321E-02, -.135,
.746, -.203, 7.244E-02, .576, .634, -8.73E-02,
-6.10E-02, .850, -.142, .129, .882, -3.86E-02,
2.063E-02, -4.15E-02, .927, -.181, -.220, .873), byrow=T, ncol=3)

  if( fuzz < max(abs(v - tst))) {
    cat("Calculated value is not the same as test value in test rotations 1. Value:\n")
    print(v, digits=18)
    cat("difference:\n")
    print(v - tst, digits=18)
    all.ok <- FALSE  
    } 

# oblimin, Kaiser normalization
# Pattern Matrix
vw <- print(oblimin(L, normalize = TRUE, eps = 1e-7))
v <- vw$loadings

tst <- matrix(c(.241, .787, -1.36E-02, 1.783E-02, .848, -.119,
-.240, .779, 6.824E-02, .608, .507, -2.52E-02,
.858, -.163, -7.26E-02, .896, 3.050E-02, 3.954E-02,
9.405E-02, 7.397E-02, .949, -8.61E-02, -.113, .875), byrow=T, ncol=3)
tst <- tst %*% matrix(c(0,1,0,1,0,0,0,0,1), 3) # Needed to line up factors correctly

fuzz <- 3e-3  #less strict; differences in 4th decimal compared to SPSS; 0.003 or smaller diff

  if( fuzz < max(abs(v - tst))) {
    cat("Calculated value is not the same as test value in test rotations 1. Value:\n")
    print(v, digits=18)
    cat("difference:\n")
    print(v - tst, digits=18)
    all.ok <- FALSE  
    } 

# oblimin, Kaiser normalization
# Structure Matrix
v <- vw$loadings %*% vw$Phi

tst <- matrix(c(.379, .829, -.146, .191, .862, -.203,
-.123, .731, .051, .701, .613, -.218, .847, -.010, -.261,
.891, .180, -.176, -.118, .000, .919, -.313, -.211, .906), byrow=T, ncol=3)
tst <- tst %*% matrix(c(0,1,0,1,0,0,0,0,1), 3) # Needed to line up factors correctly

fuzz <- 4e-3  #less strict; differences in 4th decimal compared to SPSS; 0.004 or smaller diff

  if( fuzz < max(abs(v - tst))) {
    cat("Calculated value is not the same as test value in test rotations 1. Value:\n")
    print(v, digits=18)
    cat("difference:\n")
    print(v - tst, digits=18)
    all.ok <- FALSE  
    } 

#################################################################
#
# Confirmation that a row of zeroes will not break the normalization function
# Normalizing with a column of zeroes was not affected
# based on example from Kim-Laura Speck (25 October 2023)
# Only affects Normalize=TRUE settings 

fuzz <- 1e-6

D <- matrix(c(0,0,0, 1,2,3, 2,3,4, 5,2,5, 1,2,1, 3,4,5),ncol=3,byrow=T)
set.seed(1000)  #set seed becasuse some variance is observed in converged values
v <- geominQ(D, normalize = TRUE, maxit = 10000)$loadings

tst <- matrix(c(
  0.00000000,  0.00000000, 0.00000000,
 -0.36979732, -0.13603325, 3.99622380,
  0.03102554,  0.76678245, 4.68896063,
  3.28926158,  0.01317821, 5.35447764,
 -0.02755956,  2.40311582, 0.06247234,
  0.43184841,  1.66959816, 5.38169746), ncol = 3, byrow = TRUE)

  if( fuzz < max(abs(v - tst))) {
    cat("Calculated value is not the same as test value in test rotations 1. Value:\n")
    print(v, digits=18)
    cat("difference:\n")
    print(v - tst, digits=18)
    all.ok <- FALSE  
    } 
    
