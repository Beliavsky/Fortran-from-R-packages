test_that("Endowment Cash Flows", {
    library(MortalityTables)
    mortalityTables.load("Austria_Census")

    Tarif.Endowment = InsuranceTarif$new(
        name = "Example Tariff - Endowment",
        type = "endowment",
        tarif = "E1-RP",
        desc = "An endowment with regular premiums (standard tariff)",

        mortalityTable = mort.AT.census.2011.unisex,
        i = 0.005
    )
    Contract.Endowment = InsuranceContract$new(
        tarif = Tarif.Endowment,
        age = 40, policyPeriod = 20,
        sumInsured = 10000,
        contractClosing = as.Date("2020-09-01"),
        calculate = "cashflows"
    )
    expect_equal(Contract.Endowment$Parameters$ContractData$policyPeriod, 20)
    expect_equal(Contract.Endowment$Parameters$ContractData$deferralPeriod, 0)
    expect_equal(Contract.Endowment$Parameters$ContractData$premiumPeriod, 20)


    expect_true(all(Contract.Endowment$Values$cashFlows %>% select(-premiums_advance, -survival_advance, -death_SumInsured, -death_PremiumFree) == 0))

    expect_equal(Contract.Endowment$Values$cashFlows$premiums_advance, c(rep(1, 20), 0))
    expect_equal(Contract.Endowment$Values$cashFlows$survival_advance, c(rep(0, 20), 1))
    expect_equal(Contract.Endowment$Values$cashFlows$death_SumInsured, c(rep(1, 20),0))
    expect_equal(Contract.Endowment$Values$cashFlows$death_PremiumFree, c(rep(1, 20), 0))
})
