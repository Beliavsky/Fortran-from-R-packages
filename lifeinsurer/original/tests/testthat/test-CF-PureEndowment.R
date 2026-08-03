test_that("Pure Endowment Cash Flows", {
    library(MortalityTables)
    mortalityTables.load("Austria_Census")

    Tarif.PureEnd = InsuranceTarif$new(
        name = "Example Tariff - Pure Endowment",
        type = "pureendowment",
        tarif = "PE1-RP",
        desc = "A pure endowment with regular premiums (standard tariff)",

        mortalityTable = mort.AT.census.2011.unisex,
        i = 0.005,
        premiumRefund = 1
    )
    Contract.PureEnd = InsuranceContract$new(
        tarif = Tarif.PureEnd,
        age = 40, policyPeriod = 20,
        sumInsured = 10000,
        contractClosing = as.Date("2020-09-01"),
        calculate = "cashflows"
    )
    expect_equal(Contract.PureEnd$Parameters$ContractData$policyPeriod, 20)
    expect_equal(Contract.PureEnd$Parameters$ContractData$deferralPeriod, 0)
    expect_equal(Contract.PureEnd$Parameters$ContractData$premiumPeriod, 20)


    expect_true(all(Contract.PureEnd$Values$cashFlows %>% select(-premiums_advance, -survival_advance, -death_GrossPremium, -death_Refund_past) == 0))

    expect_equal(Contract.PureEnd$Values$cashFlows$premiums_advance, c(rep(1, 20), 0))
    expect_equal(Contract.PureEnd$Values$cashFlows$survival_advance, c(rep(0, 20), 1))
    expect_equal(Contract.PureEnd$Values$cashFlows$death_GrossPremium, c(1:20,0))
    expect_equal(Contract.PureEnd$Values$cashFlows$death_Refund_past, c(rep(1, 20), 0))
})
