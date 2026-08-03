test_that("multiplication works", {
    c1 = initializeCosts(
        alpha = 0.04,
        Zillmer = 0.025,
        beta = 0.05,
        gamma = 0.0005,
        gamma.paidUp = 0.001,
        gamma.premiumfree = 0.00075,
        gamma.contract = 0.002,
        unitcosts = 10,
        unitcosts.PolicyPeriod = 12
    )

    # the above is the short form of:
    c2 = initializeCosts()
    c2[["alpha", "SumPremiums", "once"]] = 0.04
    c2[["Zillmer", "SumPremiums", "once"]] = 0.025
    c2[["beta", "GrossPremium", "PremiumPeriod"]] = 0.05
    c2[["gamma", "SumInsured", "PremiumPeriod"]] = 0.0005
    c2[["gamma_nopremiums", "SumInsured", "PolicyPeriod"]] = 0.001
    c2[["gamma", "SumInsured", "PremiumFree"]] = 0.00075
    c2[["gamma", "SumInsured", "PolicyPeriod"]] = 0.002
    c2[["unitcosts", "Constant", "PremiumPeriod"]] = 10
    c2[["unitcosts", "Constant", "PolicyPeriod"]] = 12

    expect_equal(c1, c2)
    expect_equal(c1, initializeCosts(c1))
})
