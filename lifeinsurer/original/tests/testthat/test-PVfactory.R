test_that("multiplication works", {
    library(MortalityTables)
    mortalityTables.load("Austria_Census")

    # Sample data: Austrian population mortality, interest 3%
    ags = 40:100
    v = 1/1.03
    qx = deathProbabilities(mort.AT.census.2011.unisex, ages = ags)
    qq = data.frame(x = ags, qx = qx, ix = 0, px = 1 - qx)
    pvf = PVfactory$new(qx = qq, m = 1, v = v)
    # For invalidity, simply use half the death prob (split qx into qx and ix, leave px unchanged!)
    qqI = qq
    qqI$ix = qqI$qx/2
    qqI$qx = qqI$ix
    pvfI = PVfactory$new(qx = qqI, m = 1, v = v)

    #################################### #
    # Cash Flow Definitions
    #################################### #
    # Single payment after 10 years
    cf1 = c(rep(0, 10), 1)
    # Annual payments for 10 years
    cfN = rep(1, 10)

    #################################### #
    # Guaranteed  PV
    #################################### #
    PV1adv = v^(10:0)
    PV1arr = v^(11:1)
    PVNadv = (1 - v^(10:1))/(1-v)
    PVNarr = v * PVNadv

    # Check basic present values for correctness (one-dimensional CF vector)
    expect_equal(as.vector(pvf$guaranteed(advance = cf1)), PV1adv)
    expect_equal(as.vector(pvf$guaranteed(arrears = cf1)), PV1arr)
    expect_equal(as.vector(pvf$guaranteed(advance = cfN)), PVNadv)
    expect_equal(as.vector(pvf$guaranteed(arrears = cfN)), PVNarr)

    # Same cash flows, either understood paid in advance at next timestep or in arrears of previous timestep => same present values
    expect_equal(c(pvf$guaranteed(advance = c(1,cfN))), 1 + c(pvf$guaranteed(arrears = cfN), 0))

    # Check cash flow arrays
    # PV of single payment is v^(n-t), PV of annuity is (1-v^(n-t))/(1-v)
    # Use CF array with those two cash flows
    cf2d = array(c(cf1, cfN, 0), dim = c(length(cf1),2))
    expect_equal(pvf$guaranteed(advance = cf2d), array(c(PV1adv, PVNadv, 0), dim = c(length(cf1), 2)))

    # two-dimensional cashflows at each time => 3-dimensional tensor
    cf3d = array(c(cf1, cfN, 0, cf1, cfN, 0, cfN, 0, cf1), dim = c(length(cf1), 2, 3))
    expect_equal(pvf$guaranteed(advance = cf3d), array(c(PV1adv, PVNadv, 0, PV1adv, PVNadv, 0, PVNadv, 0, PV1adv), dim = c(length(cf1), 2, 3)))



    #################################### #
    # Survival PV
    #################################### #
    PV1adv.sv = Reduce(`*`, 1-qx[1:10], init = 1, right = TRUE, accumulate = TRUE) * (v^(10:0))
    PV1arr.sv = head(Reduce(`*`, 1-qx[1:11], init = 1, right = TRUE, accumulate = TRUE), -1) * (v^(11:1))
    PVNadv.sv = head(Reduce(function(p, pv1) {pv1 * v * p + 1}, 1-qx[1:10], init = 0, right = TRUE, accumulate = TRUE), -1)
    PVNarr.sv = head(Reduce(function(p, pv1) { v * p * (1 + pv1)}, 1-qx[1:10], init = 0, right = TRUE, accumulate = TRUE), -1)

    # check basic PV
    expect_equal(as.vector(pvf$survival(advance = cf1)), PV1adv.sv)
    expect_equal(as.vector(pvf$survival(arrears = cf1)), PV1arr.sv)
    expect_equal(as.vector(pvf$survival(advance = cfN)), PVNadv.sv)
    expect_equal(as.vector(pvf$survival(arrears = cfN)), PVNarr.sv)

    # Check cash flow arrays
    expect_equal(pvf$survival(advance = cf2d), array(c(PV1adv.sv, PVNadv.sv, 0), dim = c(length(cf1), 2)))
    expect_equal(pvf$survival(arrears = cf2d), array(c(PV1arr.sv, PVNarr.sv, 0), dim = c(length(cf1), 2)))

    # two-dimensional cashflows at each time => 3-dimensional tensor
    expect_equal(pvf$survival(advance = cf3d), array(c(PV1adv.sv, PVNadv.sv, 0, PV1adv.sv, PVNadv.sv, 0, PVNadv.sv, 0, PV1adv.sv), dim = c(length(cf1), 2, 3)))
    expect_equal(pvf$survival(arrears = cf3d), array(c(PV1arr.sv, PVNarr.sv, 0, PV1arr.sv, PVNarr.sv, 0, PVNarr.sv, 0, PV1arr.sv), dim = c(length(cf1), 2, 3)))



    #################################### #
    # Death PV
    #################################### #
    PVN.death = head(Reduce(function(q, pv1) {q * v + (1 - q) * v * pv1}, qx[1:10], init = 0, right = TRUE, accumulate = TRUE), -1)
    expect_equal(as.vector(pvf$death(benefits = cfN)), PVN.death)

    #################################### #
    # Disease PV
    #################################### #
    # Death and disease probabilities are equal, so the PV should be equal. If death() is implemented correctly, this can detect errors in
    expect_equal(pvfI$disease(benefits = cfN), pvfI$death(benefits = cfN))

})
