test_that("Rounding Helper", {

    rounding = RoundingHelper$new(test = 2, gross.premium = 0, "Sum Insured" = -2)

    rounding$rounding

    expect_equal(
        rounding$round("test", c(2*10^(-8:5), 987654321.987654321)),
        c(0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.02, 0.20, 2.00, 20.00, 200.00, 2000.00, 20000.00, 200000.00, 987654321.99)
    )

    expect_equal(
        rounding$round("gross.premium", c(2*10^(-8:5), 987654321.987654321)),
        c(0, 0, 0, 0, 0, 0, 0, 0, 2, 20, 200, 2000, 20000, 200000, 987654322)
    )

    expect_equal(
        rounding$round("Sum Insured", c(2*10^(-8:5), 987654321.987654321)),
        c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 200, 2000, 20000, 200000, 987654300)
    )

    expect_equal(
        rounding$round("NotExisting", c(2*10^(-8:5), 987654321.987654321)),
        c(2*10^(-8:5), 987654321.987654321)
    )

})

test_that("Calculation of sumInsured from premium", {
    library(MortalityTables)
    mortalityTables.load("Austria_Census")

    Tarif.EndowmentSI = InsuranceTarif$new(
        name = "Example Tariff - Endowment",
        type = "endowment",
        tarif = "E1-RP",
        desc = "An endowment with regular premiums (standard tariff)",
        age = 40, policyPeriod = 20,

        mortalityTable = mort.AT.census.2011.unisex,
        cost = initializeCosts(alpha = 0.04, gamma.contract = 0.0005, unitcosts = 10),
        i = 0.03,
        sumInsured = 10000,
        contractClosing = as.Date("2020-09-01")
    )
    Tarif.EndowmentSI.rounded = Tarif.EndowmentSI$createModification(
        Rounding = list("Premium gross" = 0, "Premium net" = 2)
    )
    Contract.sumInsured = InsuranceContract$new(
        tarif = Tarif.EndowmentSI
    )
    Contract.sumInsured.rounded = InsuranceContract$new(
        tarif = Tarif.EndowmentSI.rounded
    )

    expect_equal(Contract.sumInsured.rounded$Values$premiums[["gross"]], round(Contract.sumInsured$Values$premiums[["gross"]], 0))
    expect_equal(Contract.sumInsured.rounded$Values$premiums[["net"]], round(Contract.sumInsured$Values$premiums[["net"]], 2))

})
