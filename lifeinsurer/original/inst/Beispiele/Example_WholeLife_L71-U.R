library(magrittr)
library(MortalityTables)
library(LifeInsureR)
mortalityTables.load("Austria_Census")

Tarif.L71U = InsuranceTarif$new(
    name = "L71-U",
    type = "wholelife",
    tarif = "DeathPlus - Short Term Life Insurance",
    desc = "Term Life insurance (5 years) with constant sum insured and regular premiums",
    policyPeriod = 5, premiumPeriod = 5,  # premiumPeriod not needed, defaults to maturity

    mortalityTable = mortalityTable.mixed(
        table1 = mort.AT.census.2011.male, weight1 = 0.65,
        table2 = mort.AT.census.2011.female, weight2 = 0.35
    ),
    i = 0.005,
    tax = 0.04,
    costs = initializeCosts(alpha = 0.05, gamma = 0.01, gamma.paidUp = 0.01, unitcosts = 10),
    surrenderValueCalculation = function(surrenderReserve, params, values) {
        surrenderReserve * 0.9
    }
);

contract.L71U  = InsuranceContract$new(
    Tarif.L71U,
    age = 35,
    contractClosing = as.Date("2020-08-18"),
    sumInsured = 100000);

contract.L71U$Values$premiums
contract.L71U$Values$reserves

contract.L71U = contract.L71U$premiumWaiver(t = 3)
contract.L71U$Values$reserves
contract.L71U$Values$premiumCoefficients



exportInsuranceContract.xlsx(contract.L71U, "L71U.xlsx")
openXL("L71U.xlsx")
