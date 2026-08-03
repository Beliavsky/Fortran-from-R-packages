library(testthat)
library(MultiATSM)

#  1) Load all available data from excel files
# a) Macro and bond yield data

MacroDataAll <- Load_Excel_Data(system.file("extdata", "MacroData.xlsx", package = "MultiATSM"))
YieldsDataAll <- Load_Excel_Data(system.file("extdata", "YieldsData.xlsx", package = "MultiATSM"))
Trade_Alldata <- Load_Excel_Data(system.file("extdata", "TradeData.xlsx", package = "MultiATSM"))

# Adjust trade data
DataConnect <- lapply(Trade_Alldata, function(df) {
  countries <- df[[1]]
  df <- as.data.frame(df[-1])
  rownames(df) <- countries
  df
})


# 2) Select desired inputs
Economies <- c("China", "Brazil", "Russia")
DomVar <- c("Eco_Act")
GlobalVar <- c("GBC", "CPI_OECD")
N <- 1
t0 <- "2007-09-01"
tF <- "2019-01-01"


############ 1) Test output structure ############
test_that("DatabasePrep returns correct structure", {
  ModelType <- "GVAR multi"
  FactorLabels <- LabFac(N, DomVar, GlobalVar, Economies, ModelType)
  Wgvar <- Transition_Matrix(t_First = "2006", t_Last = "2019", Economies, type = "Sample Mean", DataConnect)
  res <- DatabasePrep(t0, tF, Economies, N, FactorLabels, ModelType, MacroDataAll, YieldsDataAll, Wgvar)
  expect_type(res, "list")
  expect_true(all(Economies %in% names(res)))
  expect_true("Global" %in% names(res))
})


############ 2) Test error for missing transition matrix ############
test_that("DatabasePrep issue an error for missing Wgvar a GVAR-based model", {
  ModelType <- "GVAR multi"
  FactorLabels <- LabFac(N, DomVar, GlobalVar, Economies, ModelType)
  expect_error(DatabasePrep(t0, tF, Economies, N, FactorLabels, ModelType, MacroDataAll, YieldsDataAll, Wgvar = NULL))
})


############# 3) Test for all model types, but GVAR-based one ############
test_that("DatabasePrep checks overall structure of non GVAR-based models", {
  ModelTypeSet <- c("JPS original", "JPS global", "JPS multi", "JLL original", "JLL No DomUnit", "JLL joint Sigma")

  for (l in 1:length(ModelTypeSet)) {
    FactorLabels <- LabFac(N, DomVar, GlobalVar, Economies, ModelTypeSet[l])
    res <- DatabasePrep(t0, tF, Economies, N, FactorLabels, ModelTypeSet[l], MacroDataAll, YieldsDataAll)
    expect_type(res, "list")
    expect_true(all(Economies %in% names(res)))
    expect_true("Global" %in% names(res))

    # Check that "star" lists are all NULL
    for (i in 1:length(Economies)) {
      FacList <- res[[Economies[i]]]$Factors
      star_names <- grep("\\.Star$", names(FacList), value = TRUE)
      expect_true(any(sapply(FacList[star_names], is.null)))
    }
  }
})
