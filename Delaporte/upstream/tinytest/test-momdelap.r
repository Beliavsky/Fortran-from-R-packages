# Copyright (c) 2013, Avraham Adler All rights reserved
# SPDX-License-Identifier: BSD-2-Clause

tol <- 1e-12
testData <- c(5, 7, 9, 9, 10, 11, 11, 13, 17, 24)
inapp <- paste("Method of moments not appropriate for this data;",
               "results include non-positive parameters.")
MTD <- mean(testData)
VTD <- var(testData)
vmMTD <- VTD - MTD
SK1 <- 1.1944183061921252 # Three values from skewness in e1071 package
SK2 <- 1.4164058724628172
SK3 <- 1.0198122281732285

# Function accuracy - type 1 explicit
mom <- MoMdelap(testData, type = 1L)
P2 <- 0.5 * (SK1 * (VTD ^ 1.5) - MTD - 3 * vmMTD) / vmMTD
P1 <- vmMTD / (P2 ^ 2)
P3 <- MTD - P1 * P2
expect_equal(mom[[1]], P1, tolerance = tol)
expect_equal(mom[[2]], P2, tolerance = tol)
expect_equal(mom[[3]], P3, tolerance = tol)

# Function accuracy - type 2 implicit
mom <- MoMdelap(testData)
P2 <- 0.5 * (SK2 * (VTD ^ 1.5) - MTD - 3 * vmMTD) / vmMTD
P1 <- vmMTD / (P2 ^ 2)
P3 <- MTD - P1 * P2
expect_equal(mom[[1]], P1, tolerance = tol)
expect_equal(mom[[2]], P2, tolerance = tol)
expect_equal(mom[[3]], P3, tolerance = tol)

# Function accuracy - type 2 explicit
mom <- MoMdelap(testData, type = 2)
expect_equal(mom[[1]], 0.88342721893491116, tolerance = tol)
expect_equal(mom[[2]], 4.51388888888888928, tolerance = tol)
expect_equal(mom[[3]], 7.61230769230769155, tolerance = tol)

# Function accuracy - type 3 explicit
mom <- MoMdelap(testData, type = 3)
P2 <- 0.5 * (SK3 * (VTD ^ 1.5) - MTD - 3 * vmMTD) / vmMTD
P1 <- vmMTD / (P2 ^ 2)
P3 <- MTD - P1 * P2
expect_equal(mom[[1]], P1, tolerance = tol)
expect_equal(mom[[2]], P2, tolerance = tol)
expect_equal(mom[[3]], P3, tolerance = tol)

# MoMdelap traps bad types
expect_error(MoMdelap(testData, type = 4),
             "Skew type must be one of 1, 2, or 3.")

# MoMdelap traps bad parameters
testData <- c(3,  2, 12, 11,  1,  7,  1,  4,  0, 4)
expect_error(MoMdelap(testData), inapp)

# Non-double vector converted
expect_equal(MoMdelap(c(30L, 32L, 39L, 50L), type = 2L),
             MoMdelap(c(30, 32, 39, 50), type = 2L), tolerance = tol)
