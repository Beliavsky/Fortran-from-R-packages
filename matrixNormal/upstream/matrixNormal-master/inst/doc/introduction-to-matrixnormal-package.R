## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,  # If TRUE, all output would be in the code chunk.
  results = "markup",
  comment = NA,
  prompt = TRUE,
  strip.white = TRUE,
  tidy = TRUE,
  tidy.opts = list(width.cutoff = 90), # options for tidy to remove blank lines [blank = FALSE] and set the approximate line width to be 80.
  fig.show = "asis",
  fig.height = 4.5,  # inches
  fig.width = 4.5   # inches
)
library("formatR")

## -----------------------------------------------------------------------------
library(matrixNormal)

## -----------------------------------------------------------------------------

library(datasets)
data(USArrests)
X <- cbind(USArrests$Assault, USArrests$Murder)
Y <- cbind(USArrests$UrbanPop, USArrests$Rape)
cor(Y)

## -----------------------------------------------------------------------------
# Y is n = 50 x p = 2 that follows a matrix normal with mean matrix M, which is product of
M <- (100 * toeplitz(50:1))[, 1:2]
dim(M)
head(M)
U <- I(50) # Covariance across states: Assumed to be independent
U[1:5, 1:5]
V <- cov(X) # Covariance across predictors
V

# Find the density if Y has the density with these arguments.
matrixNormal::dmatnorm(Y, M, U, V)

## -----------------------------------------------------------------------------
# Generate a random matrix from this prior.
# The prior mean of regression matrix
J(2, 3)
# The prior variance between rape and population
t(X) %*% X
# The prior variance between regression parameters
I(3)

# Random draw for prior would have these values
A <- matrixNormal::rmatnorm(M = J(2, 3), U = t(X) %*% X, V = I(3))
A

# Predicted Counts for y can be given as:
ceiling(rowSums(X %*% A))

## -----------------------------------------------------------------------------
# Make a 3 x 3 Identity matrix
I(3)
# Make a 3 x 4 J matrix
J(3, 4)
# Make a 3 x 3 J matrix
J(3, 3)

# Calculate the trace of a J matrix
tr(J(3, 3))  # Should be 3

# Stack a matrix (used in distribution functions)
A <- matrix(c(1:4), nrow = 2, dimnames = list(NULL, c("A", "B")))
A
vec(A)

# Test if matrix is symmetric (used in distribution function)
is.symmetric.matrix(A)

## ----echo=FALSE---------------------------------------------------------------
# sessioninfo::session_info() # makes a mess! Instead

cat(" -- Session info ---------------------------------------------------")
sessioninfo::platform_info()
cat("--  Packages -------------------------------------------------------")
tmp.df <- sessioninfo::package_info(
  pkgs =c("LaplacesDemon", "MBSP", "matrixsampling", "matrixcalc"), dependencies = FALSE
  )
print(tmp.df)

