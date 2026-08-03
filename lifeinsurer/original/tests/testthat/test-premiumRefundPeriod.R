test_that("premiumRefundPeriod", {
    library(MortalityTables)
    mortalityTables.load("Austria_Annuities_AVOe2005R")

    #--------------------------------------------------------------------------- -
    # For deferred contracts, premium refund applies during deferral only by default
    #--------------------------------------------------------------------------- -
    Tarif.DefAnnuity = InsuranceTarif$new(
        type = "annuity",

        policyPeriod = function(params, values) { 120 - params$ContractData$age},
        deferralPeriod = function(params, values) { 65 - params$ContractData$age},
        premiumPeriod = function(params, values) { 65 - params$ContractData$age},
        premiumRefund = 1,
        mortalityTable = AVOe2005R.unisex
    )
    Contract.DefAnnuity = InsuranceContract$new(
        tarif = Tarif.DefAnnuity,
        age = 40, YOB = 1980,
        sumInsured = 1200,
        contractClosing = as.Date("2020-09-01"),
        calculate = "cashflows"
    )
    expect_equal(Contract.DefAnnuity$Parameters$ContractData$policyPeriod, 80)
    expect_equal(Contract.DefAnnuity$Parameters$ContractData$deferralPeriod, 25)
    expect_equal(Contract.DefAnnuity$Parameters$ContractData$premiumPeriod, 25)
    expect_equal(Contract.DefAnnuity$Parameters$ContractData$premiumRefundPeriod, 25)

    # premium refund only during the frist 25 years (linearly increasing), then 0
    expect_equal(Contract.DefAnnuity$Values$cashFlows[,"death_GrossPremium"], c(1:25, rep(0, 81-25)))
    expect_equal(Contract.DefAnnuity$Values$cashFlows[,"death_Refund_past"], c(rep(1, 25), rep(0, 81-25)))



    #--------------------------------------------------------------------------- -
    # For all other contracts without deferral period, refund period is whole contract
    #--------------------------------------------------------------------------- -
    Tarif.PureEnd = InsuranceTarif$new(
        type = "pureendowment",

        policyPeriod = 25,
        premiumRefund = 1,
        mortalityTable = AVOe2005R.unisex
    )
    Contract.PureEnd = InsuranceContract$new(
        tarif = Tarif.PureEnd,
        age = 40, YOB = 1980,
        sumInsured = 10000,
        contractClosing = as.Date("2020-09-01"),
        calculate = "cashflows"
    )
    # premium refund during the whole contract (25 years), the last year is the survival date without any death benefits
    expect_equal(Contract.PureEnd$Values$cashFlows[,"death_GrossPremium"], c(1:25, 0))
    expect_equal(Contract.PureEnd$Values$cashFlows[,"death_Refund_past"], c(rep(1, 25), 0))
})
