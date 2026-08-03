library(testthat)
library(LifeInsureR)
library(MortalityTables)
mortalityTables.load("Austria_Census")
mortalityTables.load("Austria_Annuities_AVOe2005R")

test_check("LifeInsureR")
