test_that("Dread-Disease Cash Flows", {
    library(MortalityTables)
    mortalityTables.load("Austria_Census")

    # An example dread-disease tariff, morbidity is assumed linearly increasing with age
    ddTable = mortalityTable.period(name = "Linear dread-disease table",
                                    ages = 0:100, deathProbs = 0:100/500)
    Tarif.DreadDisease = InsuranceTarif$new(
        name = "Example Tariff - Dread-Disease",
        type = "dread-disease",
        tarif = "DD1",
        desc = "A dread disease insurance with a lump-sum payment upon diagnosis",

        sumInsured = 50000,
        mortalityTable = mort.AT.census.2011.unisex,
        invalidityTable = ddTable,
        i = 0.005
    )

    Contract.DreadDisease = InsuranceContract$new(
        tarif = Tarif.DreadDisease,
        age = 40, policyPeriod = 20,
        sumInsured = 10000,
        contractClosing = as.Date("2020-09-01"),
        calculate = "cashflows"
    )
    expect_equal(Contract.DreadDisease$Parameters$ContractData$policyPeriod, 20)
    expect_equal(Contract.DreadDisease$Parameters$ContractData$deferralPeriod, 0)
    expect_equal(Contract.DreadDisease$Parameters$ContractData$premiumPeriod, 20)


    expect_true(all(Contract.DreadDisease$Values$cashFlows %>% select(-premiums_advance, -disease_SumInsured) == 0))

    expect_equal(Contract.DreadDisease$Values$cashFlows$premiums_advance, c(rep(1, 20), 0))
    expect_equal(Contract.DreadDisease$Values$cashFlows$disease_SumInsured, c(rep(1, 20), 0))
})

""
test_that("Endowment + Dread-Disease Cash Flows", {
    library(MortalityTables)
    mortalityTables.load("Austria_Census")

    # An example dread-disease tariff, morbidity is assumed linearly increasing with age
    ddTable = mortalityTable.period(name = "Linear dread-disease table",
                                    ages = 0:100, deathProbs = 0:100/500)
    Tarif.EndowDreadDisease = InsuranceTarif$new(
        name = "Example Tariff - Endowment + Dread-Disease",
        type = "endowment + dread-disease",
        tarif = "EnDD1",
        desc = "An endowment and dread disease insurance with a lump-sum payment upon diagnosis",

        sumInsured = 50000,
        mortalityTable = mort.AT.census.2011.unisex,
        invalidityTable = ddTable,
        invalidityEndsContract = TRUE,
        i = 0.005
    )

    Contract.EndowDreadDisease = InsuranceContract$new(
        tarif = Tarif.EndowDreadDisease,
        age = 40, policyPeriod = 20,
        sumInsured = 10000,
        contractClosing = as.Date("2020-09-01"),
        calculate = "cashflows"
    )
    expect_equal(Contract.EndowDreadDisease$Parameters$ContractData$policyPeriod, 20)
    expect_equal(Contract.EndowDreadDisease$Parameters$ContractData$deferralPeriod, 0)
    expect_equal(Contract.EndowDreadDisease$Parameters$ContractData$premiumPeriod, 20)


    expect_true(all(Contract.EndowDreadDisease$Values$cashFlows %>% select(-premiums_advance, -survival_advance, -death_SumInsured, -death_PremiumFree, -disease_SumInsured) == 0))

    expect_equal(Contract.EndowDreadDisease$Values$cashFlows$premiums_advance, c(rep(1, 20), 0))
    expect_equal(Contract.EndowDreadDisease$Values$cashFlows$survival_advance, c(rep(0, 20), 1))
    expect_equal(Contract.EndowDreadDisease$Values$cashFlows$death_SumInsured, c(rep(1, 20), 0))
    expect_equal(Contract.EndowDreadDisease$Values$cashFlows$disease_SumInsured, c(rep(1, 20), 0))
    expect_equal(Contract.EndowDreadDisease$Values$cashFlows$death_PremiumFree, c(rep(1, 20), 0))

})

