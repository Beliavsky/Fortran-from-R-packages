library(MortalityTables)
library(LifeInsureR)
mortalityTables.load("Austria_Census")
mortalityTables.load("Austria_Annuities_AVOe2005R")
# Costs: 4% acquisition, where 2.5% are zillmered, 5\% of each premium as beta costs,
#        1%o acquisition costs of the sum insured over the whole contract period
example.Costs = initializeCosts(
    alpha = 0.04, Zillmer = 0.025,
    beta = 0.05,
    gamma.contract = 0.001, gamma.paidUp = 0.001
)
example.Surrender = function(surrenderReserve, params, values) {
    n = params$ContractData$policyPeriod
    # Surrender Penalty is 10% at the beginning and decreases linearly to 0%
    surrenderReserve * (0.9 + 0.1 * (0:n)/n)
}


Tarif.Endowment = InsuranceTarif$new(
    name = "Example Tariff - Endowment",
    type = "endowment",
    tarif = "EN1",
    desc = "An endowment with regular premiums",

    mortalityTable = mort.AT.census.2011.unisex,
    i = 0.005,
    costs = example.Costs,
    unitcosts = 10,
    tax = 0.04,         # 4% insurance tax
    surrenderValueCalculation = example.Surrender
)

Tarif.Life = InsuranceTarif$new(
    name = "Example Tariff - Whole/Term Life",
    type = "wholelife",
    tarif = "Life1",
    desc = "A whole or term life insurance with regular premiums",

    mortalityTable = mort.AT.census.2011.unisex,
    i = 0.005,
    costs = example.Costs,
    unitcosts = 10,
    tax = 0.04,         # 4% insurance tax
    surrenderValueCalculation = example.Surrender
)

Tarif.ImmAnnuity = InsuranceTarif$new(
    name = "Example Tariff - Immediate Annuity",
    type = "annuity",
    tarif = "Ann1",
    desc = "An annuity with single-premium",
    premiumPeriod = 1,

    mortalityTable = AVOe2005R.unisex,
    i = 0.005,
    costs = example.Costs,
    tax = 0.04         # 4% insurance tax
)


# Premium periods and deferral periods can also be given as a function of other
# contract parameters (like the age at contract inception, etc.)
Tarif.DefAnnuity = InsuranceTarif$new(
    name = "Example Tariff - Deferred Annuity",
    type = "annuity",
    tarif = "Life1",
    desc = "A deferred annuity (life-long payments start at age 65) with reg. premiums",

    contractPeriod = function(params, values) { 120 - params$ContractData$age},
    deferralPeriod = function(params, values) { 65 - params$ContractData$age},
    premiumPeriod = function(params, values) { 65 - params$ContractData$age},

    mortalityTable = AVOe2005R.unisex,
    i = 0.005,
    costs = example.Costs,
    tax = 0.04,         # 4% insurance tax
    surrenderValueCalculation = example.Surrender
)

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
    i = 0.005,
    costs = example.Costs,
    unitcosts = 10,
    tax = 0.04,         # 4% insurance tax
    surrenderValueCalculation = example.Surrender
)


Tarif.PureEnd = InsuranceTarif$new(
    name = "Example Tariff - Pure Endowment",
    type = "pureendowment",
    tarif = "PE1-RP",
    desc = "A pure endowment with regular premiums (standard tariff)",

    mortalityTable = mort.AT.census.2011.unisex,
    i = 0.005,
    # Costs: 4% acquisition, where 2.5% are zillmered, 5\% of each premium as beta costs,
    #        1%o administration costs of the sum insured over the whole contract period
    costs = initializeCosts(alpha = 0.04, Zillmer = 0.025, beta = 0.05, gamma.contract = 0.001, gamma.paidUp = 0.001),
    unitcosts = 10,

    # Yearly premiums get no surcharge, monthly premiums add +4%
    premiumFrequencyLoading = list("1" = 0, "2" = 0, "4" = 0, "12" = 0.04),
    premiumRefund = 1,  # Full gross premium refund upon death
    tax = 0.04,         # 4% insurance tas

    surrenderValueCalculation = function(surrenderReserve, params, values) {
        n = params$ContractData$policyPeriod
        # Surrender Penalty is 10% at the beginning and decreases linearly to 0%
        surrenderReserve * (0.9 + 0.1 * (0:n)/n)
    }
)


contract.PureEnd = InsuranceContract$new(
    Tarif.PureEnd,
    age = 50, policyPeriod = 20,
    premiumFrequency = 2,
    sumInsured = 100000,
    initialCapital = 50000,
    contractClosing = as.Date("2020-10-01"),
    balanceSheetDate = as.Date("2020-09-30")
)

contract.PureEnd$Values$cashFlows
contract.PureEnd$Values$premiums
contract.PureEnd$Values$reservesBalanceSheet

exportInsuranceContract.xlsx(contract.PureEnd, "test.xlsx")
openxlsx::openXL("text.xlsx")
