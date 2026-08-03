test_that("Annuity Cash Flows", {
    library(MortalityTables)
    mortalityTables.load("Austria_Annuities_AVOe2005R")

    Tarif.Annuity = InsuranceTarif$new(
        name = "Example Tariff - Immediate Annuity",
        type = "annuity",
        tarif = "Ann1",
        desc = "An annuity with single-premium",
        premiumPeriod = 1,

        mortalityTable = AVOe2005R.unisex,
        i = 0.005
    )
    Contract.Annuity = InsuranceContract$new(
        tarif = Tarif.Annuity,
        age = 65, YOB = 1955,
        sumInsured = 1200,
        policyPeriod = 55,
        contractClosing = as.Date("2020-09-01"),
        calculate = "cashflows"
    )
    expect_equal(Contract.Annuity$Parameters$ContractData$policyPeriod, 55)
    expect_equal(Contract.Annuity$Parameters$ContractData$deferralPeriod, 0)
    expect_equal(Contract.Annuity$Parameters$ContractData$premiumPeriod, 1)


    expect_true(all(Contract.Annuity$Values$cashFlows %>% select(-premiums_advance, -survival_advance) == 0))

    # 1 year premium cash flow
    expect_equal(Contract.Annuity$Values$cashFlows$premiums_advance, c(1, rep(0, 55)))
    # premium payment start immediately for the whole contract period
    expect_equal(Contract.Annuity$Values$cashFlows$survival_advance, c(rep(1, 55), 0))
})
