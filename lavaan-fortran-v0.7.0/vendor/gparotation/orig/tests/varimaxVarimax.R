# Also see the first test in rotations.R
# compares varimax to Varimax to 0.001 discrepancy

 require("GPArotation")

### note that this is a slightly lower bar than other tests
### to correct for the built-in varimax function working differently
### than GPA, and to ensure Varimax convergence
### these are differences in the 4th decimal or better
 fuzz <- 1e-4
###
 all.ok <- TRUE  

sortFac <- function(x){  # Based on Fungible faSort 
  vx <- order(colSums(x$loadings^2), decreasing = TRUE)
  Dsgn <- diag(sign(colSums(x$loadings^3))) [ , vx]
  x$Th <- x$Th %*% Dsgn 
  x$loadings <- x$loadings %*% Dsgn
  if ("Phi" %in% names(x)) { 
    x$Phi <- diag(1/diag(Dsgn)) %*% x$Phi %*% Dsgn
    }
  x 
} 
data(Thurstone, package="GPArotation")
yv1 <- varimax(box26, normalize = FALSE, eps = 1e-7) #built-in R
names(yv1) <- c("loadings","Th")
yv1 <- sortFac(yv1)
yv2 <- sortFac(Varimax(box26, normalize = FALSE, maxit = 10000, eps = 1e-7)) #GPArotation version
# yv.diff <- unclass(yv1$loadings) - unclass(yv2$loadings)
# max(abs(yv.diff))


 if( fuzz < max(abs(yv1$loadings - yv2$loadings))) {
    cat("Calculated varimax is not the same as Varimax:\n")
    # print(yv2$loadings, digits=18)
    cat("difference:\n")
    print(unclass(yv1$loadings - yv2$loadings), digits=18)
    all.ok <- FALSE  
    } 

