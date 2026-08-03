test_that("Calculation of sumInsured from premium", {
    library(MortalityTables)
    mortalityTables.load("Austria_Census")

    Tarif.EndowmentSI = InsuranceTarif$new(
        name = "Example Tariff - Endowment",
        type = "endowment",
        tarif = "E1-RP",
        desc = "An endowment with regular premiums (standard tariff)",

        mortalityTable = mort.AT.census.2011.unisex,
        cost = initializeCosts(alpha = 0.04, gamma.contract = 0.0005, unitcosts = 10),
        i = 0.03
    )
    Contract.sumInsured = InsuranceContract$new(
        tarif = Tarif.EndowmentSI,
        age = 40, policyPeriod = 20,
        sumInsured = 10000,
        contractClosing = as.Date("2020-09-01")
    )
    Contract.premium = InsuranceContract$new(
        tarif = Tarif.EndowmentSI,
        age = 40, policyPeriod = 20,
        premium = Contract.sumInsured$Values$premiums[["written"]],
        contractClosing = as.Date("2020-09-01")
    )
    Contract.premium_beforetax = InsuranceContract$new(
        tarif = Tarif.EndowmentSI,
        age = 40, policyPeriod = 20,
        premium = c(written_beforetax = Contract.sumInsured$Values$premiums[["written_beforetax"]]),
        contractClosing = as.Date("2020-09-01")
    )
    Contract.premium_gross = InsuranceContract$new(
        tarif = Tarif.EndowmentSI,
        age = 40, policyPeriod = 20,
        premium = c(gross = Contract.sumInsured$Values$premiums[["gross"]]),
        contractClosing = as.Date("2020-09-01")
    )

    # All four contracts above should result in the same sumInsured:
    expect_equal(Contract.premium$Parameters$ContractData$sumInsured, Contract.sumInsured$Parameters$ContractData$sumInsured)
    expect_equal(Contract.premium_beforetax$Parameters$ContractData$sumInsured, Contract.sumInsured$Parameters$ContractData$sumInsured)
    expect_equal(Contract.premium_gross$Parameters$ContractData$sumInsured, Contract.sumInsured$Parameters$ContractData$sumInsured)
    expect_equal(Contract.sumInsured$Parameters$ContractData$sumInsured, Contract.premium$Parameters$ContractData$sumInsured)

})
