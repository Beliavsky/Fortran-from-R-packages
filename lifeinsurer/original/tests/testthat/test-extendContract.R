test_that("Extend contract by $addExtension", {
    library(MortalityTables)
    mortalityTables.load("Austria_Census")

    Tarif.EndowmentA = InsuranceTarif$new(
        name = "Example Tariff - Endowment 1981",
        type = "endowment",
        tarif = "E1-RP81",
        desc = "An endowment with regular premiums (standard tariff)",

        mortalityTable = mort.AT.census.1981.male,
        cost = initializeCosts(alpha = 0.04, gamma.contract = 0.0005, unitcosts = 10),
        i = 0.03
    )
    Tarif.EndowmentB = Tarif.EndowmentA$createModification(
        name = "Example Tariff - Endowment 2001",
        tarif = "E1-RP01",
        mortalityTable = mort.AT.census.2001.male,
        cost = initializeCosts(alpha = 0.024, gamma.contract = 0.00075, unitcosts = 20),
        i = 0.01)

    ContractA = InsuranceContract$new(
        tarif = Tarif.EndowmentA,
        age = 40, policyPeriod = 20,
        sumInsured = 10000,
        contractClosing = as.Date("2000-07-01")
    )


    # premium-free extension
    ContractB = ContractA$clone()$addExtension(id = "Verlaengerung1", contractPeriod = 5, premiumPeriod = 0)
    expect_equal(ContractB$blocks$Verlaengerung1$Parameters$ContractData$sumInsured, 15117.03896)


    # extension with given sumInsured resulting in 0 premiums
    ContractC = ContractA$clone()$addExtension(id = "Verlaengerung1", contractPeriod = 5, sumInsured = 15117.03896)
    expect_equal(ContractC$blocks$Verlaengerung1$Values$premiums[["gross"]], 0, tolerance = 1e-06)

    # extension with increased sumInsured: real premiums are charged, reserves start from the existing reserve:
    ContractD = ContractA$clone()$addExtension(id = "Verlaengerung1", contractPeriod = 5, sumInsured = 20000)
    expect_equal(ContractD$blocks$Verlaengerung1$Values$premiums[["written"]], 315.109)
    expect_equal(ContractD$blocks$Verlaengerung1$Values$reserves[["0", "contractual"]], 10000)

    # extension with increased sumInsured and different tariff: check whether interest rate has really changed
    ContractE = ContractA$clone()$addExtension(id = "Verlaengerung1", contractPeriod = 5, sumInsured = 20000, tarif = Tarif.EndowmentB)
    expect_equal(unname(ContractE$Values$basicData[c(1,20,30,40), "InterestRate"]), c(0.03, 0.03, 0.01, 0.01))
})
