library(testthat)
library(MultiATSM)


# 1) Test basic label generation
test_that("LabFac generates correct labels for typical input", {
  N <- 2
  DomVar <- c("inflation", "Output gap")
  GlobalVar <- "Commodity Prices"
  Economies <- c("U.S.", "Canada")
  ModelType <- "JPS original"
  res <- LabFac(N, DomVar, GlobalVar, Economies, ModelType)
  expect_type(res, "list")
  expect_equal(res$Spanned, c("Level", "Slope"))
  expect_equal(res$Domestic, c("inflation", "Output gap", "Level", "Slope"))
  expect_equal(res$Global, "Commodity Prices")
  expect_equal(res$Tables$AllCountries, c(
    "inflation U.S.", "Output gap U.S.", "Level U.S.", "Slope U.S.",
    "inflation Canada", "Output gap Canada", "Level Canada", "Slope Canada"
  ))
})


# 2) Test error for invalid N
test_that("LabFac issues an error for invalid (large) N", {
  N <- 10
  DomVar <- c("inflation")
  GlobalVar <- "Commodity Prices"
  Economies <- c("U.S.")
  ModelType <- "JPS original"
  expect_error(LabFac(N, DomVar, GlobalVar, Economies, ModelType))
})


# 3) Test JLL model type labels
test_that("LabFac generates JLL labels when ModelType is JLL", {
  N <- 1
  DomVar <- c("inflation")
  GlobalVar <- "Commodity Prices"
  Economies <- c("U.S.")
  ModelType <- "JLL original"
  res <- LabFac(N, DomVar, GlobalVar, Economies, ModelType)
  expect_type(res, "list")
  expect_equal(res$Spanned, c("Level"))
  expect_equal(res$Domestic, c("inflation", "Level"))
  expect_equal(res$Tables$AllCountriesJLL, c("inflation U.S. JLL", "Level U.S. JLL"))
})
