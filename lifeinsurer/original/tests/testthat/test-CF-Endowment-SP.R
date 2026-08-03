test_that("Endowment Cash Flows", {
    library(MortalityTables)
    mortalityTables.load("Austria_Census")

    Tarif.Endowment.SP = InsuranceTarif$new(
        name = "Example Tariff - Endowment Single Premium",
        type = "endowment",
        tarif = "E1-SP",
        desc = "An endowment with single premiums (standard tariff)",
        premiumPeriod = 1,
        mortalityTable = mort.AT.census.2011.unisex,
        i = 0.005
    )
    Contract.Endowment.SP = InsuranceContract$new(
        tarif = Tarif.Endowment.SP,
        age = 40, policyPeriod = 20,
        sumInsured = 10000,
        contractClosing = as.Date("2020-09-01"),
        calculate = "cashflows"
    )
    expect_equal(Contract.Endowment.SP$Parameters$ContractData$policyPeriod, 20)
    expect_equal(Contract.Endowment.SP$Parameters$ContractData$deferralPeriod, 0)
    expect_equal(Contract.Endowment.SP$Parameters$ContractData$premiumPeriod, 1)


    expect_true(all(Contract.Endowment.SP$Values$cashFlows %>% select(-premiums_advance, -survival_advance, -death_SumInsured, -death_PremiumFree) == 0))

    expect_equal(Contract.Endowment.SP$Values$cashFlows$premiums_advance, c(1, rep(0, 20)))
    expect_equal(Contract.Endowment.SP$Values$cashFlows$survival_advance, c(rep(0, 20), 1))
    expect_equal(Contract.Endowment.SP$Values$cashFlows$death_SumInsured, c(rep(1, 20),0))
    expect_equal(Contract.Endowment.SP$Values$cashFlows$death_PremiumFree, c(rep(1, 20), 0))
})
