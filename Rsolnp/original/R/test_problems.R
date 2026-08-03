hs01_problem <- function()
{
    fn <- function(x) {
        100 * (x[2] - x[1]^2)^2 + (1 - x[1])^2
    }
    gr <- function(x) {
        g2 <- 200 * (x[2] - x[1]^2)
        g1 <- -2 * (x[1] * (g2 - 1) + 1)
        c(g1, g2)
    }
    lower <- c(-1000, -1.5)
    upper <- c(1000, 1000)
    start <- c(-2, 1)
    eq_fn <- NULL
    eq_jac <- NULL
    eq_b <- NULL
    ineq_fn <- NULL
    ineq_jac <- NULL
    ineq_lower <- NULL
    ineq_upper <- NULL
    best_fn <- 0
    best_par <- c(1, 1)
    list(
        name = "hs01",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = best_fn,
        best_par = best_par
    )
}

hs02_problem <- function()
{
    fn <- function(x) {
        100 * (x[2] - x[1]^2)^2 + (1 - x[1])^2
    }
    gr <- function(x) {
        dfdx1 <- -400 * x[1] * (x[2] - x[1]^2) - 2 * (1 - x[1])
        dfdx2 <- 200 * (x[2] - x[1]^2)
        c(dfdx1, dfdx2)
    }
    lower <- c(-1000, 1.5)
    upper <- c( 1000,  1.5)
    start <- c(-2, 1)
    # Known solution from Fortran (XEX and FEX)
    w1 <- sqrt(598/1200)
    best_par <- c(2 * w1 * cos(acos(2.5e-3 / w1^3) / 3), 1.5)
    best_fn <- 100 * (best_par[2] - best_par[1]^2)^2 + (1 - best_par[1])^2
    list(
        name = "hs02",
        fn = fn,
        gr = gr,
        eq_fn = NULL,
        eq_b = NULL,
        eq_jac = NULL,
        ineq_fn = NULL,
        ineq_jac = NULL,
        ineq_lower = NULL,
        ineq_upper = NULL,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = best_fn,
        best_par = best_par
    )
}



hs03_problem <- function()
{
    fn <- function(x) {
        x[2] + 0.00001 * (x[2] - x[1])^2
    }
    gr <- function(x) {
        dfdx1 <- -2e-5 * (x[2] - x[1])
        dfdx2 <- 1 + 2e-5 * (x[2] - x[1])
        c(dfdx1, dfdx2)
    }
    lower <- c(-1000, 0)
    upper <- c( 1000, 1000)
    start <- c(10, 1)
    eq_fn <- NULL
    eq_jac <- NULL
    ineq_fn <- NULL
    ineq_jac <- NULL
    eq_b <- NULL
    ineq_lower <- NULL
    ineq_upper <- NULL
    list(
        name = "hs03",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = 0,
        best_par = c(0, 0)
    )
}

hs04_problem <- function()
{
    fn <- function(x) {
        ((x[1] + 1)^3) / 3 + x[2]
    }
    gr <- function(x) {
        dfdx1 <- (x[1] + 1)^2
        dfdx2 <- 1
        c(dfdx1, dfdx2)
    }
    lower <- c(1, 0)
    upper <- c(1000, 1000)
    start <- c(1.125, 0.125)
    eq_fn <- NULL
    eq_jac <- NULL
    ineq_fn <- NULL
    ineq_jac <- NULL
    eq_b <- NULL
    ineq_lower <- NULL
    ineq_upper <- NULL
    list(
        name = "hs04",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = 8/3,
        best_par = c(1, 0)
    )
}

hs05_problem <- function()
{
    fn <- function(x) {
        sin(x[1] + x[2]) + (x[1] - x[2])^2 - 1.5 * x[1] + 2.5 * x[2] + 1
    }
    gr <- function(x) {
        dfdx1 <- cos(x[1] + x[2]) + 2 * (x[1] - x[2]) - 1.5
        dfdx2 <- cos(x[1] + x[2]) - 2 * (x[1] - x[2]) + 2.5
        c(dfdx1, dfdx2)
    }
    lower <- c(-1.5, -3)
    upper <- c(4, 3)
    start <- c(0, 0)
    eq_fn <- NULL
    eq_jac <- NULL
    ineq_fn <- NULL
    ineq_jac <- NULL
    eq_b <- NULL
    ineq_lower <- NULL
    ineq_upper <- NULL
    list(
        name = "hs05",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = -(sqrt(3)/2 + pi/3),
        best_par = c(0.5 - pi/3, -0.5 - pi/3)
    )
}

hs06_problem <- function()
{
    fn <- function(x) {
        (1 - x[1])^2
    }
    gr <- function(x) {
        dfdx1 <- -2 * (1 - x[1])
        dfdx2 <- 0
        c(dfdx1, dfdx2)
    }
    eq_fn <- function(x) {
        10 * (x[2] - x[1]^2)
    }
    eq_jac <- function(x) {
        matrix(c(-20 * x[1], 10), nrow = 1)
    }
    lower <- c(-1000, -1000)
    upper <- c( 1000,  1000)
    start <- c(-1.2, 1)
    eq_b <- 0
    ineq_fn <- NULL
    ineq_jac <- NULL
    ineq_lower <- NULL
    ineq_upper <- NULL
    list(
        name = "hs06",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = 0,
        best_par = c(1, 1)
    )
}

hs07_problem <- function()
{
    fn <- function(x) {
        log(1 + x[1]^2) - x[2]
    }
    gr <- function(x) {
        dfdx1 <- 2 * x[1] / (1 + x[1]^2)
        dfdx2 <- -1
        c(dfdx1, dfdx2)
    }
    eq_fn <- function(x) {
        (1 + x[1]^2)^2 + x[2]^2 - 4
    }
    eq_jac <- function(x) {
        matrix(c(4 * x[1] * (1 + x[1]^2), 2 * x[2]), nrow = 1)
    }
    lower <- c(-1000, -1000)
    upper <- c(1000, 1000)
    start <- c(2, 2)
    eq_b <- 0
    ineq_fn <- NULL
    ineq_jac <- NULL
    ineq_lower <- NULL
    ineq_upper <- NULL
    list(
        name = "hs07",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = -sqrt(3),
        best_par = c(0, sqrt(3))
    )
}

s08_problem <- function()
{
    fn <- function(x) {
        -1
    }
    gr <- function(x) {
        c(0, 0)
    }
    eq_fn <- function(x) {
        c(x[1]^2 + x[2]^2 - 25, x[1] * x[2] - 9)
    }
    eq_jac <- function(x) {
        matrix(c(2 * x[1], 2 * x[2], x[2], x[1]), nrow = 2, byrow = TRUE)
    }
    lower <- c(-1000, -1000)
    upper <- c(1000, 1000)
    start <- c(2, 1)
    eq_b <- c(0, 0)
    ineq_fn <- NULL
    ineq_jac <- NULL
    ineq_lower <- NULL
    ineq_upper <- NULL
    # List all four solutions as best_par if you want:
    A <- sqrt((25 + sqrt(301))/2)
    B <- sqrt((25 - sqrt(301))/2)
    best_par <- list(
        c( A,  9/A),
        c(-A, -9/A),
        c( B,  9/B),
        c(-B, -9/B)
    )
    list(
        name = "hs08",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = -1,
        best_par = NULL
    )
}

hs09_problem <- function()
{
    fn <- function(x) {
        sin(pi * x[1] / 12) * cos(pi * x[2] / 16)
    }
    gr <- function(x) {
        v3 <- pi / 12
        v4 <- pi / 16
        v1 <- v3 * x[1]
        v2 <- v4 * x[2]
        dfdx1 <- v3 * cos(v1) * cos(v2)
        dfdx2 <- -v4 * sin(v1) * sin(v2)
        c(dfdx1, dfdx2)
    }
    eq_fn <- function(x) {
        4 * x[1] - 3 * x[2]
    }
    eq_jac <- function(x) {
        matrix(c(4, -3), nrow = 1)
    }
    lower <- c(-1000, -1000)
    upper <- c(1000, 1000)
    start <- c(0, 0)
    eq_b <- 0
    ineq_fn <- NULL
    ineq_jac <- NULL
    ineq_lower <- NULL
    ineq_upper <- NULL
    list(
        name = "hs09",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = -0.5,
        best_par = c(-3, -4)
    )
}

hs10_problem <- function()
{
    fn <- function(x) {
        x[1] - x[2]
    }
    gr <- function(x) {
        c(1, -1)
    }
    ineq_fn <- function(x) {
        -3 * x[1]^2 + 2 * x[1] * x[2] - x[2]^2 + 1
    }
    ineq_jac <- function(x) {
        matrix(c(-6 * x[1] + 2 * x[2], 2 * (x[1] - x[2])), nrow = 1)
    }
    lower <- c(-1000, -1000)
    upper <- c(1000, 1000)
    start <- c(-10, 10)
    eq_fn <- NULL
    eq_jac <- NULL
    eq_b <- NULL
    ineq_lower <- 0
    ineq_upper <- 1e8
    list(
        name = "hs10",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = -1,
        best_par = c(0, 1)
    )
}

hs11_problem <- function()
{
    fn <- function(x) {
        (x[1] - 5)^2 + x[2]^2 - 25
    }
    gr <- function(x) {
        c(2 * (x[1] - 5), 2 * x[2])
    }
    ineq_fn <- function(x) {
        -x[1]^2 + x[2]
    }
    ineq_jac <- function(x) {
        matrix(c(-2 * x[1], 1), nrow = 1)
    }
    lower <- c(-1000, -1000)
    upper <- c(1000, 1000)
    start <- c(4.9, 0.1)
    eq_fn <- NULL
    eq_jac <- NULL
    eq_b <- NULL
    ineq_lower <- 0
    ineq_upper <- 1e10
    # Solution per Fortran code
    AEX <- 7.5 * sqrt(6)
    AW <- (sqrt(AEX^2 + 1) + AEX)^(1/3)
    QAW <- AW^2
    x1sol <- (AW - 1/AW) / sqrt(6)
    x2sol <- (QAW - 2 + 1/QAW) / 6
    best_par <- c(x1sol, x2sol)
    best_fn <- (x1sol - 5)^2 + x2sol^2 - 25
    list(
        name = "hs11",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = best_fn,
        best_par = best_par
    )
}

hs12_problem <- function()
{
    fn <- function(x) {
        0.5 * x[1]^2 + x[2]^2 - x[1] * x[2] - 7 * x[1] - 7 * x[2]
    }
    gr <- function(x) {
        c(x[1] - x[2] - 7, 2 * x[2] - x[1] - 7)
    }
    ineq_fn <- function(x) {
        25 - 4 * x[1]^2 - x[2]^2
    }
    ineq_jac <- function(x) {
        matrix(c(-8 * x[1], -2 * x[2]), nrow = 1)
    }
    lower <- c(-1000, -1000)
    upper <- c(1000, 1000)
    start <- c(0, 0)
    eq_fn <- NULL
    eq_jac <- NULL
    eq_b <- NULL
    ineq_lower <- 0
    ineq_upper <- 1e8
    list(
        name = "hs12",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = -30,
        best_par = c(2, 3)
    )
}

hs13_problem <- function()
{
    fn <- function(x) {
        (x[1] - 2)^2 + x[2]^2
    }
    gr <- function(x) {
        c(2 * (x[1] - 2), 2 * x[2])
    }
    ineq_fn <- function(x) {
        (1 - x[1])^3 - x[2]
    }
    ineq_jac <- function(x) {
        matrix(c(-3 * (1 - x[1])^2, -1), nrow = 1)
    }
    lower <- c(0, 0)
    upper <- c(1000, 1000)
    start <- c(0, 0)
    eq_fn <- NULL
    eq_jac <- NULL
    eq_b <- NULL
    ineq_lower <- 0
    ineq_upper <- 1e8
    list(
        name = "hs13",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = 1,
        best_par = c(1, 0)
    )
}
hs14_problem <- function()
{
    fn <- function(x) {
        (x[1] - 2)^2 + (x[2] - 1)^2
    }
    gr <- function(x) {
        c(2 * (x[1] - 2), 2 * (x[2] - 1))
    }
    ineq_fn <- function(x) {
        1 - 0.25 * x[1]^2 - x[2]^2
    }
    ineq_jac <- function(x) {
        matrix(c(-0.5 * x[1], -2 * x[2]), nrow = 1)
    }
    eq_fn <- function(x) {
        x[1] - 2 * x[2] + 1
    }
    eq_jac <- function(x) {
        matrix(c(1, -2), nrow = 1)
    }
    lower <- c(-1000, -1000)
    upper <- c(1000, 1000)
    start <- c(2, 2)
    ineq_lower <- 0
    ineq_upper <- 1e8
    eq_b <- 0
    w7 <- sqrt(7)
    x1sol <- (w7 - 1) * 0.5
    x2sol <- (w7 + 1) * 0.25
    best_par <- c(x1sol, x2sol)
    best_fn <- 9 - 23 * w7 / 8
    list(
        name = "hs14",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = best_fn,
        best_par = best_par
    )
}

hs15_problem <- function()
{
    fn <- function(x) {
        (x[2] - x[1]^2)^2 + 0.01 * (1 - x[1])^2
    }
    gr <- function(x) {
        g2 <- 2 * (x[2] - x[1]^2)
        g1 <- -0.02 * (x[1] * (g2 - 1) + 1)
        c(g1, g2)
    }
    ineq_fn <- function(x) {
        c(x[1] * x[2] - 1, x[2]^2 + x[1])
    }
    ineq_jac <- function(x) {
        matrix(c(x[2], x[1], 1,2 * x[2]), nrow = 2, byrow = TRUE)
    }
    lower <- c(-1000, -1000)
    upper <- c(0.5, 1000)
    start <- c(-2, 1)
    eq_fn <- NULL
    eq_jac <- NULL
    eq_b <- NULL
    ineq_lower <- c(0, 0)
    ineq_upper <- c(1e8, 1e8)
    list(
        name = "hs15",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = 3.065,
        best_par = c(0.5, 2.0)
    )
}

hs16_problem <- function()
{
    fn <- function(x) {
        100 * (x[2] - x[1]^2)^2 + (1 - x[1])^2
    }
    gr <- function(x) {
        g2 <- 200 * (x[2] - x[1]^2)
        g1 <- -2 * (x[1] * (g2 - 1) + 1)
        c(g1, g2)
    }
    ineq_fn <- function(x) {
        c(x[2]^2 + x[1], x[1]^2 + x[2])
    }
    ineq_jac <- function(x) {
        matrix(c(1, 2 * x[2], 2 * x[1], 1), nrow = 2, byrow = TRUE)
    }
    lower <- c(-0.5, -1000)
    upper <- c(0.5, 1)
    start <- c(0, 0.5)
    eq_fn <- NULL
    eq_jac <- NULL
    eq_b <- NULL
    ineq_lower <- c(0, 0)
    ineq_upper <- c(1e8, 1e8)
    list(
        name = "hs16",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = 0.25,
        best_par = c(0.5, 0.25)
    )
}

hs17_problem <- function()
{
    fn <- function(x) {
        100 * (x[2] - x[1]^2)^2 + (1 - x[1])^2
    }
    gr <- function(x) {
        g2 <- 200 * (x[2] - x[1]^2)
        g1 <- -2 * (x[1] * (g2 - 1) + 1)
        c(g1, g2)
    }
    ineq_fn <- function(x) {
        c(x[2]^2 - x[1], x[1]^2 - x[2])
    }
    ineq_jac <- function(x) {
        matrix(c(-1, 2 * x[2], 2 * x[1], -1), nrow = 2, byrow = TRUE)
    }

    lower <- c(-2, -1000)
    upper <- c(0.5, 1)
    start <- c(-1, 0)
    eq_fn <- NULL
    eq_jac <- NULL
    eq_b <- NULL
    ineq_lower <- c(0, 0)
    ineq_upper <- c(1e8, 1e8)
    list(
        name = "hs17",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = 1,
        best_par = c(0, 0)
    )
}

hs18_problem <- function()
{
    fn <- function(x) {
        0.01 * x[1]^2 + x[2]^2
    }
    gr <- function(x) {
        c(0.02 * x[1], 2 * x[2])
    }
    ineq_fn <- function(x) {
        c(x[1] * x[2] - 25, x[1]^2 + x[2]^2 - 25)
    }
    ineq_jac <- function(x) {
        matrix(c(x[2], x[1], 2 * x[1], 2 * x[2]), nrow = 2, byrow = TRUE)
    }
    lower <- c(2, 0)
    upper <- c(50, 50)
    start <- c(10, 2)
    eq_fn <- NULL
    eq_jac <- NULL
    eq_b <- NULL
    ineq_lower <- c(0, 0)
    ineq_upper <- c(1e8, 1e8)
    x1sol <- sqrt(250)
    x2sol <- 0.1 * x1sol
    best_par <- c(x1sol, x2sol)
    best_fn <- 5
    list(
        name = "hs18",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = best_fn,
        best_par = best_par
    )
}


hs19_problem <- function()
{
    fn <- function(x) {
        (x[1] - 10)^3 + (x[2] - 20)^3
    }
    gr <- function(x) {
        c(3 * (x[1] - 10)^2, 3 * (x[2] - 20)^2)
    }
    ineq_fn <- function(x) {
        c((x[1] - 5)^2 + (x[2] - 5)^2 - 100, 82.81 - (x[1] - 6)^2 - (x[2] - 5)^2)
    }
    ineq_jac <- function(x) {
        matrix(c(2 * (x[1] - 5), 2 * (x[2] - 5), -2 * (x[1] - 6), -2 * (x[2] - 5)), nrow = 2, byrow = TRUE)
    }
    lower <- c(13, 0)
    upper <- c(100, 100)
    start <- c(20.1, 5.84)
    eq_fn <- NULL
    eq_jac <- NULL
    eq_b <- NULL
    ineq_lower <- c(0, 0)
    ineq_upper <- c(1e8, 1e8)
    saex <- 17.280975
    aex <- sqrt(saex)
    x1sol <- 14.095
    x2sol <- 5 - aex
    best_par <- c(x1sol, x2sol)
    best_fn <- (4.095)^3 - (15 + aex)^3
    list(
        name = "hs19",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = best_fn,
        best_par = best_par
    )
}

hs20_problem <- function()
{
    fn <- function(x) {
        100 * (x[2] - x[1]^2)^2 + (1 - x[1])^2
    }
    gr <- function(x) {
        g2 <- 200 * (x[2] - x[1]^2)
        g1 <- -2 * (x[1] * (g2 - 1) + 1)
        c(g1, g2)
    }
    ineq_fn <- function(x) {
        c(x[2]^2 + x[1], x[1]^2 + x[2], x[1]^2 + x[2]^2 - 1)
    }
    ineq_jac <- function(x) {
        matrix(c(1, 2 * x[2], 2 * x[1], 1, 2 * x[1], 2 * x[2]), nrow = 3, byrow = TRUE)
    }
    lower <- c(-0.5, -1000)
    upper <- c(0.5, 1000)
    start <- c(0.1, 1)
    eq_fn <- NULL
    eq_jac <- NULL
    eq_b <- NULL
    ineq_lower <- c(0, 0, 0)
    ineq_upper <- c(1e8, 1e8, 1e8)
    x1sol <- 0.5
    x2sol <- sqrt(3) * 0.5
    best_par <- c(x1sol, x2sol)
    best_fn <- 81.5 - 25 * sqrt(3)
    list(
        name = "hs20",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = best_fn,
        best_par = best_par
    )
}

hs21_problem <- function()
{
    fn <- function(x) {
        0.01 * x[1]^2 + x[2]^2 - 100
    }
    gr <- function(x) {
        c(0.02 * x[1], 2 * x[2])
    }
    ineq_fn <- function(x) {
        10 * x[1] - x[2] - 10
    }
    ineq_jac <- function(x) {
        matrix(c(10, -1), nrow = 1)
    }
    lower <- c(2, -50)
    upper <- c(50, 50)
    start <- c(4, -1)
    eq_fn <- NULL
    eq_jac <- NULL
    eq_b <- NULL
    ineq_lower <- 0
    ineq_upper <- 1e8
    best_fn <- -99.96
    best_par <- c(2, 0)
    list(
        name = "hs21",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = best_fn,
        best_par = best_par
    )
}

hs22_problem <- function()
{
    fn <- function(x) {
        (x[1] - 2)^2 + (x[2] - 1)^2
    }
    gr <- function(x) {
        c(2 * (x[1] - 2), 2 * (x[2] - 1))
    }
    ineq_fn <- function(x) {
        2 - x[1] - x[2]
    }
    ineq_jac <- function(x) {
        matrix(c(-1, -1), nrow = 1)
    }
    eq_fn <- function(x) {
        x[2] - x[1]^2
    }
    eq_jac <- function(x) {
        matrix(c(-2 * x[1], 1), nrow = 1)
    }
    lower <- c(-1000, -1000)
    upper <- c(1000, 1000)
    start <- c(2, 2)
    eq_b <- 0
    ineq_lower <- 0
    ineq_upper <- 1e8
    list(
        name = "hs22",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = 1,
        best_par = c(1, 1)
    )
}

hs23_problem <- function()
{
    fn <- function(x) {
        x[1]^2 + x[2]^2
    }
    gr <- function(x) {
        c(2 * x[1], 2 * x[2])
    }
    ineq_fn <- function(x) {
        c(x[1] + x[2] - 1, x[1]^2 + x[2]^2 - 1, 9 * x[1]^2 + x[2]^2 - 9, x[1]^2 - x[2], x[2]^2 - x[1])
    }
    ineq_jac <- function(x) {
        matrix(c(1, 1, 2 * x[1],  2 * x[2], 18 * x[1], 2 * x[2], 2 * x[1], -1, -1, 2 * x[2]), nrow = 5, byrow = TRUE)
    }
    lower <- c(-50, -50)
    upper <- c(50, 50)
    start <- c(3, 1)
    eq_fn <- NULL
    eq_jac <- NULL
    eq_b <- NULL
    ineq_lower <- rep(0, 5)
    ineq_upper <- rep(1e8, 5)
    list(
        name = "hs23",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = 2,
        best_par = c(1, 1)
    )
}

hs24_problem <- function()
{
    fn <- function(x) {
        A <- sqrt(3)
        ((x[1] - 3)^2 - 9) * x[2]^3 / (27 * A)
    }
    gr <- function(x) {
        A <- sqrt(3)
        g1 <- 2 * (x[1] - 3) * x[2]^3 / (27 * A)
        g2 <- ((x[1] - 3)^2 - 9) * x[2]^2 / (9 * A)
        c(g1, g2)
    }
    ineq_fn <- function(x) {
        A <- sqrt(3)
        c(x[1] / A - x[2], x[1] + x[2] * A, 6 - x[2] * A - x[1])
    }
    ineq_jac <- function(x) {
        A <- sqrt(3)
        matrix(c(1/A, -1, 1, A, -1, -A), nrow = 3, byrow = TRUE)
    }
    lower <- c(0, 0)
    upper <- c(1000, 1000)
    start <- c(500, 500)
    eq_fn <- NULL
    eq_jac <- NULL
    eq_b <- NULL
    ineq_lower <- rep(0, 3)
    ineq_upper <- rep(1e8, 3)
    A <- sqrt(3)
    best_par <- c(3, A)
    best_fn <- -1
    list(
        name = "hs24",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = best_fn,
        best_par = best_par
    )
}

hs25_problem <- function()
{
    fn <- function(x) {
        s <- 1:99
        u <- 25 + (-50 * log(s / 100))^(2/3)
        z <- u - x[2]
        if (any(z <= 0)) {
            return(sum((x - 5)^2))
        }
        exp_term <- exp(-(z^x[3]) / x[1])
        sum((-s / 100 + exp_term)^2)
    }
    gr <- function(x) {
        s <- 1:99
        u <- 25 + (-50 * log(s / 100))^(2/3)
        z <- u - x[2]
        if (any(z <= 0)) {
            return(2 * (x - 5))
        }
        e <- exp(-(z^x[3]) / x[1])
        r <- -s / 100 + e

        de_dx1 <-  (z^x[3]) / (x[1]^2) * e
        de_dx2 <- -x[3] * z^(x[3] - 1) / x[1] * e
        de_dx3 <- -log(z) * z^x[3] / x[1] * e

        grad1 <- sum(2 * r * de_dx1)
        grad2 <- sum(2 * r * de_dx2)
        grad3 <- sum(2 * r * de_dx3)

        c(grad1, grad2, grad3)
    }
    lower <- c(0.1, 1e-5, 1e-5)
    upper <- c(100, 25.6, 5)
    start <- c(100, 12.5, 3)
    eq_fn <- NULL
    eq_jac <- NULL
    ineq_fn <- NULL
    ineq_jac <- NULL
    eq_b <- NULL
    ineq_lower <- NULL
    ineq_upper <- NULL
    list(
        name = "hs25",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = 0,
        best_par = c(50, 25, 1.5)
    )
}

hs26_problem <- function()
{
    fn <- function(x) {
        (x[1] - x[2])^2 + (x[2] - x[3])^4
    }
    gr <- function(x) {
        g1 <- 2 * (x[1] - x[2])
        g3 <- -4 * (x[2] - x[3])^3
        g2 <- -g1 - g3
        c(g1, g2, g3)
    }
    eq_fn <- function(x) {
        x[1] * (1 + x[2]^2) + x[3]^4 - 3
    }
    eq_jac <- function(x) {
        matrix(c(1 + x[2]^2, 2 * x[1] * x[2], 4 * x[3]^3), nrow = 1)
    }
    lower <- rep(-1000, 3)
    upper <- rep(1000, 3)
    start <- c(0, 0, 0)
    ineq_fn <- NULL
    ineq_jac <- NULL
    eq_b <- 0
    ineq_lower <- NULL
    ineq_upper <- NULL
    list(
        name = "hs26",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = 0,
        best_par = c(1, 1, 1)
    )
}

hs27_problem <- function()
{
    fn <- function(x) {
        (x[1] - 1)^2 + 100 * (x[2] - x[1]^2)^2
    }
    gr <- function(x) {
        g1 <- 2 * (x[1] - 1) - 400 * (x[2] - x[1]^2) * x[1]
        g2 <- 200 * (x[2] - x[1]^2)
        g3 <- 0
        c(g1, g2, g3)
    }
    eq_fn <- function(x) {
        x[1] + x[3]^2 + 1
    }
    eq_jac <- function(x) {
        matrix(c(1, 0, 2 * x[3]), nrow = 1)
    }
    lower <- rep(-1000, 3)
    upper <- rep(1000, 3)
    start <- c(2, 2, 2)
    ineq_fn <- NULL
    ineq_jac <- NULL
    eq_b <- 0
    ineq_lower <- NULL
    ineq_upper <- NULL
    list(
        name = "hs27",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = 4,
        best_par = c(-1, 1, 0)
    )
}

hs28_problem <- function()
{
    fn <- function(x) {
        (x[1] + x[2])^2 + (x[2] + x[3])^2
    }
    gr <- function(x) {
        g1 <- 2 * (x[1] + x[2])
        g3 <- 2 * (x[2] + x[3])
        g2 <- g1 + g3
        c(g1, g2, g3)
    }
    eq_fn <- function(x) {
        x[1] + 2 * x[2] + 3 * x[3] - 1
    }
    eq_jac <- function(x) {
        matrix(c(1, 2, 3), nrow = 1)
    }
    lower <- rep(-1000, 3)
    upper <- rep(1000, 3)
    start <- c(-4, 1, 1)
    eq_b <- 0
    ineq_fn <- NULL
    ineq_jac <- NULL
    ineq_lower <- NULL
    ineq_upper <- NULL
    list(
        name = "hs28",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = 0,
        best_par = c(0.5, -0.5, 0.5)
    )
}


hs29_problem <- function()
{
    fn <- function(x) {
        -x[1] * x[2] * x[3]
    }
    gr <- function(x) {
        c(-x[2] * x[3],
          -x[1] * x[3],
          -x[1] * x[2])
    }
    ineq_fn <- function(x) {
        48 - x[1]^2 - 2 * x[2]^2 - 4 * x[3]^2
    }
    ineq_jac <- function(x) {
        matrix(c(-2 * x[1], -4 * x[2], -8 * x[3]), nrow = 1)
    }
    lower <- rep(-1000, 3)
    upper <- rep(1000, 3)
    start <- c(1, 1, 1)
    eq_fn <- NULL
    eq_jac <- NULL
    eq_b <- NULL
    ineq_lower <- 0
    ineq_upper <- 1e8
    best_fn <- -16 * sqrt(2)
    best_par <- c(4, 2 * sqrt(2), 2)
    list(
        name = "hs29",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = best_fn,
        best_par = best_par
    )
}

hs30_problem <- function()
{
    fn <- function(x) {
        x[1]^2 + x[2]^2 + x[3]^2
    }
    gr <- function(x) {
        c(2 * x[1], 2 * x[2], 2 * x[3])
    }
    ineq_fn <- function(x) {
        x[1]^2 + x[2]^2 - 1
    }
    ineq_jac <- function(x) {
        matrix(c(2 * x[1], 2 * x[2], 0), nrow = 1)
    }
    lower <- c(1, -10, -10)
    upper <- c(10, 10, 10)
    start <- c(1, 1, 1)
    eq_fn <- NULL
    eq_jac <- NULL
    eq_b <- NULL
    ineq_lower <- 0
    ineq_upper <- 1e8
    best_fn <- 1
    best_par <- c(1, 0, 0)
    list(
        name = "hs30",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = best_fn,
        best_par = best_par
    )
}


hs31_problem <- function()
{
    fn <- function(x) {
        9 * x[1]^2 + x[2]^2 + 9 * x[3]^2
    }
    gr <- function(x) {
        c(18 * x[1], 2 * x[2], 18 * x[3])
    }
    ineq_fn <- function(x) {
        x[1] * x[2] - 1
    }
    ineq_jac <- function(x) {
        matrix(c(x[2], x[1], 0), nrow = 1)
    }
    lower <- c(-10, 1, -10)
    upper <- c(10, 10, 1)
    start <- c(1, 1, 1)
    eq_fn <- NULL
    eq_jac <- NULL
    eq_b <- NULL
    ineq_lower <- 0
    ineq_upper <- 1e8
    best_fn <- 6
    best_par <- c(1 / sqrt(3), sqrt(3), 0)
    list(
        name = "hs31",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = best_fn,
        best_par = best_par
    )
}

hs32_problem <- function()
{
    fn <- function(x) {
        (x[1] + 3 * x[2] + x[3])^2 + 4 * (x[1] - x[2])^2
    }
    gr <- function(x) {
        g1 <- 10 * x[1] - 2 * x[2] + 2 * x[3]
        g2 <- -2 * x[1] + 26 * x[2] + 6 * x[3]
        g3 <- 2 * (x[1] + 3 * x[2] + x[3])
        c(g1, g2, g3)
    }
    ineq_fn <- function(x) {
        -x[1]^3 + 6 * x[2] + 4 * x[3] - 3
    }
    ineq_jac <- function(x) {
        matrix(c(-3 * x[1]^2, 6, 4), nrow = 1)
    }
    eq_fn <- function(x) {
        1 - x[1] - x[2] - x[3]
    }
    eq_jac <- function(x) {
        matrix(c(-1, -1, -1), nrow = 1)
    }
    lower <- rep(0, 3)
    upper <- rep(1000, 3)
    start <- c(0.1, 0.7, 0.2)
    eq_b <- 0
    ineq_lower <- 0
    ineq_upper <- 1e8
    best_fn <- 1
    best_par <- c(0, 0, 1)
    list(
        name = "hs32",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = best_fn,
        best_par = best_par
    )
}

hs33_problem <- function()
{
    fn <- function(x) {
        (x[1] - 1) * (x[1] - 2) * (x[1] - 3) + x[3]
    }
    gr <- function(x) {
        g1 <- 3 * x[1]^2 - 12 * x[1] + 11
        g2 <- 0
        g3 <- 1
        c(g1, g2, g3)
    }
    ineq_fn <- function(x) {
        c(x[3]^2 - x[1]^2 - x[2]^2, x[1]^2 + x[2]^2 + x[3]^2 - 4)
    }

    ineq_jac <- function(x) {
        matrix(c(-2 * x[1], -2 * x[2],  2 * x[3], 2 * x[1],  2 * x[2],  2 * x[3]), nrow = 2, byrow = TRUE)
    }
    lower <- c(0, 0, 0)
    upper <- c(1000, 1000, 5)
    start <- c(1, 1, 3)
    eq_fn <- NULL
    eq_jac <- NULL
    eq_b <- NULL
    ineq_lower <- c(0, 0)
    ineq_upper <- c(1e8, 1e8)
    best_fn <- sqrt(2) - 6
    best_par <- c(0, sqrt(2), sqrt(2))
    list(
        name = "hs33",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = best_fn,
        best_par = best_par
    )
}

hs34_problem <- function()
{
    fn <- function(x) {
        -x[1]
    }
    gr <- function(x) {
        c(-1, 0, 0)
    }
    ineq_fn <- function(x) {
        c(x[2] - exp(x[1]), x[3] - exp(x[2]))
    }
    ineq_jac <- function(x) {
        matrix(c(-exp(x[1]), 1, 0, 0, -exp(x[2]), 1), nrow = 2, byrow = TRUE)
    }
    lower <- c(0, 0, 0)
    upper <- c(100, 100, 10)
    start <- c(1, 1.05, 2.9)
    eq_fn <- NULL
    eq_jac <- NULL
    eq_b <- NULL
    ineq_lower <- c(0, 0)
    ineq_upper <- c(1e8, 1e8)
    best_par <- c(log(log(10)), log(10), 10)
    best_fn <- -log(log(10))
    list(
        name = "hs34",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = best_fn,
        best_par = best_par
    )
}

hs35_problem <- function()
{
    fn <- function(x) {
        9 - 8 * x[1] - 6 * x[2] - 4 * x[3] + 2 * x[1]^2 + 2 * x[2]^2 + x[3]^2 + 2 * x[1] * x[2] + 2 * x[1] * x[3]
    }
    gr <- function(x) {
        c(-8 + 4 * x[1] + 2 * x[2] + 2 * x[3], -6 + 4 * x[2] + 2 * x[1], -4 + 2 * x[3] + 2 * x[1])
    }
    ineq_fn <- function(x) {
        -x[1] - x[2] - 2 * x[3] + 3
    }
    ineq_jac <- function(x) {
        matrix(c(-1, -1, -2), nrow = 1)
    }
    lower <- c(0, 0, 0)
    upper <- rep(1000, 3)
    start <- rep(0.5, 3)
    eq_fn <- NULL
    eq_jac <- NULL
    eq_b <- NULL
    ineq_lower <- 0
    ineq_upper <- 1e8
    best_fn <- 1/9
    best_par <- c(4/3, 7/9, 4/9)
    list(
        name = "hs35",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = best_fn,
        best_par = best_par
    )
}

hs36_problem <- function()
{
    fn <- function(x) {
        -x[1] * x[2] * x[3]
    }
    gr <- function(x) {
        c(-x[2] * x[3], -x[1] * x[3], -x[1] * x[2])
    }
    ineq_fn <- function(x) {
        72 - x[1] - 2 * x[2] - 2 * x[3]
    }
    ineq_jac <- function(x) {
        matrix(c(-1, -2, -2), nrow = 1)
    }
    lower <- c(0, 0, 0)
    upper <- c(20, 11, 42)
    start <- c(10, 10, 10)
    eq_fn <- NULL
    eq_jac <- NULL
    eq_b <- NULL
    ineq_lower <- 0
    ineq_upper <- 1e8
    best_fn <- -3300
    best_par <- c(20, 11, 15)
    list(
        name = "hs36",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = best_fn,
        best_par = best_par
    )
}

hs37_problem <- function()
{
    fn <- function(x) {
        -x[1] * x[2] * x[3]
    }
    gr <- function(x) {
        c(-x[2] * x[3], -x[1] * x[3], -x[1] * x[2])
    }
    ineq_fn <- function(x) {
        c(72 - x[1] - 2 * x[2] - 2 * x[3], x[1] + 2 * x[2] + 2 * x[3])
    }
    ineq_jac <- function(x) {
        matrix(c(-1, -2, -2, 1,  2,  2), nrow = 2, byrow = TRUE)
    }
    lower <- rep(0, 3)
    upper <- rep(42, 3)
    start <- rep(10, 3)
    eq_fn <- NULL
    eq_jac <- NULL
    eq_b <- NULL
    ineq_lower <- c(0, 0)
    ineq_upper <- c(1e8, 1e8)
    best_fn <- -3456
    best_par <- c(24, 12, 12)
    list(
        name = "hs37",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = best_fn,
        best_par = best_par
    )
}

hs38_problem <- function()
{
    fn <- function(x) {
        (100 * (x[2] - x[1]^2)^2 + (1 - x[1])^2 +
             90 * (x[4] - x[3]^2)^2 + (1 - x[3])^2 +
             10.1 * ((x[2] - 1)^2 + (x[4] - 1)^2) +
             19.8 * (x[2] - 1) * (x[4] - 1))
    }
    gr <- function(x) {
        c(
            -400 * x[1] * (x[2] - x[1]^2) - 2 * (1 - x[1]),
            200 * (x[2] - x[1]^2) + 20.2 * (x[2] - 1) + 19.8 * (x[4] - 1),
            -360 * x[3] * (x[4] - x[3]^2) - 2 * (1 - x[3]),
            180 * (x[4] - x[3]^2) + 20.2 * (x[4] - 1) + 19.8 * (x[2] - 1)
        )
    }
    lower <- rep(-10, 4)
    upper <- rep(10, 4)
    start <- c(-3, -1, -3, -1)
    eq_fn <- NULL
    eq_jac <- NULL
    eq_b <- NULL
    ineq_fn <- NULL
    ineq_jac <- NULL
    ineq_lower <- NULL
    ineq_upper <- NULL
    best_fn <- 0
    best_par <- c(1, 1, 1, 1)
    list(
        name = "hs38",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = best_fn,
        best_par = best_par
    )
}

hs39_problem <- function()
{
    fn <- function(x) {
        -x[1]
    }
    gr <- function(x) {
        c(-1, 0, 0, 0)
    }
    eq_fn <- function(x) {
        c(x[2] - x[1]^3 - x[3]^2, x[1]^2 - x[2] - x[4]^2)
    }
    eq_jac <- function(x) {
        matrix(c(-3 * x[1]^2, 1, -2 * x[3], 0, 2 * x[1], -1, 0, -2 * x[4]), nrow = 2, byrow = TRUE)
    }
    lower <- rep(-1000, 4)
    upper <- rep(1000, 4)
    start <- rep(2, 4)
    eq_b <- c(0, 0)
    ineq_fn <- NULL
    ineq_jac <- NULL
    ineq_lower <- NULL
    ineq_upper <- NULL
    best_fn <- -1
    best_par <- c(1, 1, 0, 0)
    list(
        name = "hs39",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = best_fn,
        best_par = best_par
    )
}

hs40_problem <- function()
{
    fn <- function(x) {
        -x[1] * x[2] * x[3] * x[4]
    }
    gr <- function(x) {
        c(-x[2] * x[3] * x[4], -x[1] * x[3] * x[4], -x[1] * x[2] * x[4], -x[1] * x[2] * x[3])
    }
    eq_fn <- function(x) {
        c(x[1]^3 + x[2]^2 - 1, x[1]^2 * x[4] - x[3], x[4]^2 - x[2])
    }
    eq_jac <- function(x) {
        matrix(c(3 * x[1]^2, 2 * x[2], 0, 0,
                 2 * x[1] * x[4], 0, -1, x[1]^2,
                 0, -1, 0, 2 * x[4]), nrow = 3, byrow = TRUE)
    }
    lower <- rep(-1000, 4)
    upper <- rep(1000, 4)
    start <- rep(0.8, 4)
    eq_b <- c(0, 0, 0)
    ineq_fn <- NULL
    ineq_jac <- NULL
    ineq_lower <- NULL
    ineq_upper <- NULL
    best_fn <- -0.25
    best_par <- c(2^(-1/3), 2^(-0.5), 2^(-11/12), 2^(-0.25))
    list(
        name = "hs40",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = best_fn,
        best_par = best_par
    )
}

hs41_problem <- function()
{
    fn <- function(x) {
        2 - x[1] * x[2] * x[3]
    }
    gr <- function(x) {
        c(-x[2] * x[3], -x[1] * x[3], -x[1] * x[2], 0)
    }
    eq_fn <- function(x) {
        x[1] + 2 * x[2] + 2 * x[3] - x[4]
    }
    eq_jac <- function(x) {
        matrix(c(1, 2, 2, -1), nrow = 1)
    }
    lower <- c(0, 0, 0, 0)
    upper <- c(1, 1, 1, 2)
    start <- rep(1, 4)
    eq_b <- 0
    ineq_fn <- NULL
    ineq_jac <- NULL
    ineq_lower <- NULL
    ineq_upper <- NULL
    best_fn <- 52 / 27
    best_par <- c(2/3, 1/3, 1/3, 2)
    list(
        name = "hs41",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = best_fn,
        best_par = best_par
    )
}

hs42_problem <- function()
{
    fn <- function(x) {
        (x[1] - 1)^2 + (x[2] - 2)^2 + (x[3] - 3)^2 + (x[4] - 4)^2
    }
    gr <- function(x) {
        2 * (x - 1:4)
    }
    eq_fn <- function(x) {
        c(x[1] - 2,  x[3]^2 + x[4]^2 - 2)
    }
    eq_jac <- function(x) {
        matrix(c(1, 0, 0, 0, 0, 0, 2 * x[3], 2 * x[4]), nrow = 2, byrow = TRUE)
    }
    lower <- rep(-1000, 4)
    upper <- rep(1000, 4)
    start <- rep(1, 4)
    eq_b <- c(0, 0)
    ineq_fn <- NULL
    ineq_jac <- NULL
    ineq_lower <- NULL
    ineq_upper <- NULL
    best_fn <- 28 - 10 * sqrt(2)
    best_par <- c(2, 2, sqrt(0.72), sqrt(1.28))
    list(
        name = "hs42",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = best_fn,
        best_par = best_par
    )
}

hs43_problem <- function()
{
    fn <- function(x) {
        x[1]^2 + x[2]^2 + 2 * x[3]^2 + x[4]^2 -
            5 * x[1] - 5 * x[2] - 21 * x[3] + 7 * x[4]
    }
    gr <- function(x) {
        c(2 * x[1] - 5, 2 * x[2] - 5, 4 * x[3] - 21, 2 * x[4] + 7)
    }
    ineq_fn <- function(x) {
        c(-x[1]^2 - x[2]^2 - x[3]^2 - x[4]^2 - x[1] + x[2] - x[3] + x[4] + 8,
          -x[1]^2 - 2 * x[2]^2 - x[3]^2 - 2 * x[4]^2 + x[1] + x[4] + 10,
          -2 * x[1]^2 - x[2]^2 - x[3]^2 - 2 * x[1] + x[2] + x[4] + 5)
    }
    ineq_jac <- function(x) {
        matrix(c(-2 * x[1] - 1, -2 * x[2] + 1, -2 * x[3] - 1, -2 * x[4] + 1,
                 -2 * x[1] + 1, -4 * x[2], -2 * x[3], -4 * x[4] + 1,
                 -4 * x[1] - 2, -2 * x[2] + 1, -2 * x[3], 1), nrow = 3, byrow = TRUE)
    }
    lower <- rep(-1000, 4)
    upper <- rep(1000, 4)
    start <- c(0, 0, 0, 0)
    eq_fn <- NULL
    eq_jac <- NULL
    eq_b <- NULL
    ineq_lower <- rep(0, 3)
    ineq_upper <- rep(1e8, 3)
    best_fn <- -44
    best_par <- c(0, 1, 2, -1)
    list(
        name = "hs43",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = best_fn,
        best_par = best_par
    )
}

hs44_problem <- function()
{
    fn <- function(x) {
        x[1] - x[2] - x[3] - x[1] * x[3] + x[1] * x[4] + x[2] * x[3] - x[2] * x[4]
    }
    gr <- function(x) {
        c(1 - x[3] + x[4], -1 + x[3] - x[4], -1 - x[1] + x[2], x[1] - x[2])
    }
    ineq_fn <- function(x) {
        c(8 - x[1] - 2 * x[2], 12 - 4 * x[1] - x[2], 12 - 3 * x[1] - 4 * x[2], 8 - 2 * x[3] - x[4],
          8 - x[3] - 2 * x[4], 5 - x[3] - x[4])
    }
    ineq_jac <- function(x) {
        matrix(c(-1, -2,  0,  0,
                 -4, -1,  0,  0,
                 -3, -4,  0,  0,
                 0,  0, -2, -1,
                 0,  0, -1, -2,
                 0,  0, -1, -1), nrow = 6, byrow = TRUE)
    }
    lower <- rep(0, 4)
    upper <- rep(1000, 4)
    start <- c(1, 1, 1, 1)
    eq_fn <- NULL
    eq_jac <- NULL
    eq_b <- NULL
    ineq_lower <- rep(0, 6)
    ineq_upper <- rep(1e8, 6)
    best_fn <- -15
    best_par <- c(0, 3, 0, 4)
    list(
        name = "hs44",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = best_fn,
        best_par = best_par
    )
}

hs45_problem <- function()
{
    fn <- function(x) {
        2 - x[1] * x[2] * x[3] * x[4] * x[5] / 120
    }
    gr <- function(x) {
        c(-x[2] * x[3] * x[4] * x[5] / 120,
          -x[1] * x[3] * x[4] * x[5] / 120,
          -x[1] * x[2] * x[4] * x[5] / 120,
          -x[1] * x[2] * x[3] * x[5] / 120,
          -x[1] * x[2] * x[3] * x[4] / 120)
    }
    lower <- rep(0, 5)
    upper <- 1:5
    start <- rep(1, 5)
    eq_fn <- NULL
    eq_jac <- NULL
    eq_b <- NULL
    ineq_fn <- NULL
    ineq_jac <- NULL
    ineq_lower <- NULL
    ineq_upper <- NULL
    best_fn <- 1
    best_par <- 1:5
    list(
        name = "hs45",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = best_fn,
        best_par = best_par
    )
}

hs46_problem <- function()
{
    fn <- function(x) {
        (x[1] - x[2])^2 + (x[3] - 1)^2 + (x[4] - 1)^4 + (x[5] - 1)^6
    }
    gr <- function(x) {
        g1 <- 2 * (x[1] - x[2])
        g2 <- -g1
        g3 <- 2 * (x[3] - 1)
        g4 <- 4 * (x[4] - 1)^3
        g5 <- 6 * (x[5] - 1)^5
        c(g1, g2, g3, g4, g5)
    }
    eq_fn <- function(x) {
        c(x[1]^2 * x[4] + sin(x[4] - x[5]) - 1, x[2] + x[3]^4 * x[4]^2 - 2)
    }
    eq_jac <- function(x) {
        matrix(c(2 * x[1] * x[4], 0, 0, x[1]^2 + cos(x[4] - x[5]), -cos(x[4] - x[5]),
                 0, 1, 4 * x[3]^3 * x[4]^2, 2 * x[3]^4 * x[4], 0), nrow = 2, byrow = TRUE)
    }
    lower <- rep(-1000, 5)
    upper <- rep(1000, 5)
    start <- c(0.5 * sqrt(2), 1.75, 0.5, 2, 2)
    eq_b <- c(0, 0)
    ineq_fn <- NULL
    ineq_jac <- NULL
    ineq_lower <- NULL
    ineq_upper <- NULL
    best_fn <- 0
    best_par <- rep(1, 5)
    list(
        name = "hs46",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = best_fn,
        best_par = best_par
    )
}

hs47_problem <- function()
{
    fn <- function(x) {
        (x[1] - x[2])^2 + (x[2] - x[3])^2 + (x[3] - x[4])^4 + (x[4] - x[5])^4
    }
    gr <- function(x) {
        v1 <- 2 * (x[1] - x[2])
        v2 <- 2 * (x[2] - x[3])
        v3 <- 4 * (x[3] - x[4])^3
        v4 <- 4 * (x[4] - x[5])^3
        c(v1, -v1 + v2, -v2 + v3, -v3 + v4, -v4)
    }
    eq_fn <- function(x) {
        c(x[1] + x[2]^2 + x[3]^3 - 3, x[2] - x[3]^2 + x[4] - 1, x[1] * x[5] - 1)
    }
    eq_jac <- function(x) {
        matrix(c(1, 2 * x[2], 3 * x[3]^2, 0, 0,
                 0, 1, -2 * x[3], 1, 0,
                x[5], 0, 0, 0, x[1]), nrow = 3, byrow = TRUE)
    }
    lower <- rep(-1000, 5)
    upper <- rep(1000, 5)
    start <- c(2, sqrt(2), -1, 2 - sqrt(2), 0.5)
    eq_b <- c(0, 0, 0)
    ineq_fn <- NULL
    ineq_jac <- NULL
    ineq_lower <- NULL
    ineq_upper <- NULL
    best_fn <- 0
    best_par <- rep(1, 5)
    list(
        name = "hs47",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = best_fn,
        best_par = best_par
    )
}

hs48_problem <- function()
{
    fn <- function(x) {
        (x[1] - 1)^2 + (x[2] - x[3])^2 + (x[4] - x[5])^2
    }
    gr <- function(x) {
        g1 <- 2 * (x[1] - 1)
        g2 <- 2 * (x[2] - x[3])
        g3 <- -g2
        g4 <- 2 * (x[4] - x[5])
        g5 <- -g4
        c(g1, g2, g3, g4, g5)
    }
    eq_fn <- function(x) {
        c(sum(x) - 5, x[3] - 2 * (x[4] + x[5]) + 3)
    }
    eq_jac <- function(x) {
        matrix(c(1, 1, 1, 1, 1,
                 0, 0, 1, -2, -2), nrow = 2, byrow = TRUE)
    }
    lower <- rep(-1000, 5)
    upper <- rep(1000, 5)
    start <- c(3, 5, -3, 2, -2)
    eq_b <- c(0, 0)
    ineq_fn <- NULL
    ineq_jac <- NULL
    ineq_lower <- NULL
    ineq_upper <- NULL
    best_fn <- 0
    best_par <- rep(1, 5)
    list(
        name = "hs48",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = best_fn,
        best_par = best_par
    )
}

hs49_problem <- function()
{
    fn <- function(x) {
        (x[1] - x[2])^2 + (x[3] - 1)^2 + (x[4] - 1)^4 + (x[5] - 1)^6
    }
    gr <- function(x) {
        g1 <- 2 * (x[1] - x[2])
        g2 <- -g1
        g3 <- 2 * (x[3] - 1)
        g4 <- 4 * (x[4] - 1)^3
        g5 <- 6 * (x[5] - 1)^5
        c(g1, g2, g3, g4, g5)
    }
    eq_fn <- function(x) {
        c(x[1] + x[2] + x[3] + 4 * x[4] - 7, x[3] + 5 * x[5] - 6)
    }
    eq_jac <- function(x) {
        matrix(c(1, 1, 1, 4, 0,
                 0, 0, 1, 0, 5), nrow = 2, byrow = TRUE)
    }
    lower <- rep(-1000, 5)
    upper <- rep(1000, 5)
    start <- c(10, 7, 2, -3, 0.8)
    eq_b <- c(0, 0)
    ineq_fn <- NULL
    ineq_jac <- NULL
    ineq_lower <- NULL
    ineq_upper <- NULL
    best_fn <- 0
    best_par <- rep(1, 5)
    list(
        name = "hs49",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = best_fn,
        best_par = best_par
    )
}

hs50_problem <- function()
{
    fn <- function(x) {
        (x[1] - x[2])^2 + (x[2] - x[3])^2 + (x[3] - x[4])^4 + (x[4] - x[5])^4
    }
    gr <- function(x) {
        v1 <- 2 * (x[1] - x[2])
        v2 <- 2 * (x[2] - x[3])
        v3 <- 4 * (x[3] - x[4])^3
        v4 <- 4 * (x[4] - x[5])^3
        c(v1, -v1 + v2, -v2 + v3, -v3 + v4, -v4)
    }
    eq_fn <- function(x) {
        c(x[1] + 2 * x[2] + 3 * x[3] - 6,
          x[2] + 2 * x[3] + 3 * x[4] - 6,
          x[3] + 2 * x[4] + 3 * x[5] - 6
        )
    }
    eq_jac <- function(x) {
        matrix(c(1, 2, 3, 0, 0,
                 0, 1, 2, 3, 0,
                 0, 0, 1, 2, 3), nrow = 3, byrow = TRUE)
    }
    lower <- rep(-1000, 5)
    upper <- rep(1000, 5)
    start <- c(35, -31, 11, 5, -5)
    eq_b <- c(0, 0, 0)
    ineq_fn <- NULL
    ineq_jac <- NULL
    ineq_lower <- NULL
    ineq_upper <- NULL
    best_fn <- 0
    best_par <- rep(1, 5)
    list(
        name = "hs50",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = best_fn,
        best_par = best_par
    )
}

hs51_problem <- function()
{
    fn <- function(x) {
        (x[1] - x[2])^2 + (x[2] + x[3] - 2)^2 + (x[4] - 1)^2 + (x[5] - 1)^2
    }
    gr <- function(x) {
        g1 <- 2 * (x[1] - x[2])
        g3 <- 2 * (x[2] + x[3] - 2)
        g2 <- g3 - g1
        g4 <- 2 * (x[4] - 1)
        g5 <- 2 * (x[5] - 1)
        c(g1, g2, g3, g4, g5)
    }
    eq_fn <- function(x) {
        c(x[1] + 3 * x[2] - 4, x[3] + x[4] - 2 * x[5], x[2] - x[5])
    }
    eq_jac <- function(x) {
        matrix(c(1, 3, 0, 0, 0,
                 0, 0, 1, 1, -2,
                 0, 1, 0, 0, -1), nrow = 3, byrow = TRUE)
    }
    lower <- rep(-1000, 5)
    upper <- rep(1000, 5)
    start <- c(2.5, 0.5, 2, -1, 0.5)
    eq_b <- c(0, 0, 0)
    ineq_fn <- NULL
    ineq_jac <- NULL
    ineq_lower <- NULL
    ineq_upper <- NULL
    best_fn <- 0
    best_par <- rep(1, 5)
    list(
        name = "hs51",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = best_fn,
        best_par = best_par
    )
}

hs52_problem <- function()
{
    fn <- function(x) {
        (4 * x[1] - x[2])^2 + (x[2] + x[3] - 2)^2 + (x[4] - 1)^2 + (x[5] - 1)^2
    }
    gr <- function(x) {
        g1 <- 8 * (4 * x[1] - x[2])
        g3 <- 2 * (x[2] + x[3] - 2)
        g2 <- -0.25 * g1 + g3
        g4 <- 2 * (x[4] - 1)
        g5 <- 2 * (x[5] - 1)
        c(g1, g2, g3, g4, g5)
    }
    eq_fn <- function(x) {
        c(x[1] + 3 * x[2], x[3] + x[4] - 2 * x[5], x[2] - x[5])
    }
    eq_jac <- function(x) {
        matrix(c(1, 3, 0, 0, 0,
                 0, 0, 1, 1, -2,
                 0, 1, 0, 0, -1), nrow = 3, byrow = TRUE)
    }
    lower <- rep(-1000, 5)
    upper <- rep(1000, 5)
    start <- rep(2, 5)
    eq_b <- c(0, 0, 0)
    ineq_fn <- NULL
    ineq_jac <- NULL
    ineq_lower <- NULL
    ineq_upper <- NULL
    best_fn <- 1859/349
    best_par <- c(-33/349, 11/349, 180/349, -158/349, 11/349)
    list(
        name = "hs52",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = best_fn,
        best_par = best_par
    )
}

hs53_problem <- function()
{
    fn <- function(x) {
        (x[1] - x[2])^2 + (x[2] + x[3] - 2)^2 + (x[4] - 1)^2 + (x[5] - 1)^2
    }
    gr <- function(x) {
        g1 <- 2 * (x[1] - x[2])
        g3 <- 2 * (x[2] + x[3] - 2)
        g2 <- g3 - g1
        g4 <- 2 * (x[4] - 1)
        g5 <- 2 * (x[5] - 1)
        c(g1, g2, g3, g4, g5)
    }
    eq_fn <- function(x) {
        c(x[1] + 3 * x[2], x[3] + x[4] - 2 * x[5], x[2] - x[5])
    }
    eq_jac <- function(x) {
        matrix(c(1, 3, 0, 0, 0,
                 0, 0, 1, 1, -2,
                 0, 1, 0, 0, -1), nrow = 3, byrow = TRUE)
    }
    lower <- rep(-10, 5)
    upper <- rep(10, 5)
    start <- rep(2, 5)
    eq_b <- c(0, 0, 0)
    ineq_fn <- NULL
    ineq_jac <- NULL
    ineq_lower <- NULL
    ineq_upper <- NULL
    best_fn <- 176/43
    best_par <- c(-33/43, 11/43, 27/43, -5/43, 11/43)
    list(
        name = "hs53",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = best_fn,
        best_par = best_par
    )
}

hs54_problem <- function()
{
    fn <- function(x) {
        v1 <- x[1] - 1e4
        v2 <- x[2] - 1
        v3 <- x[3] - 2e6
        v4 <- x[4] - 10
        v5 <- x[5] - 1e-3
        v6 <- x[6] - 1e8
        v7 <- 1 / 0.96
        v8 <- 1 / 4.9e13
        v9 <- 1 / 2.45e13
        q <- (1.5625e-8 * v1^2 + 5e-5 * v1 * v2 + v2^2) * v7 + v3^2 * v8 +
            4e-4 * v4^2 + 4e2 * v5^2 + 4e-18 * v6^2
        -exp(-0.5 * q)
    }
    gr <- function(x) {
        v1 <- x[1] - 1e4
        v2 <- x[2] - 1
        v3 <- x[3] - 2e6
        v4 <- x[4] - 10
        v5 <- x[5] - 1e-3
        v6 <- x[6] - 1e8
        v7 <- 1 / 0.96
        v8 <- 1 / 4.9e13
        v9 <- 1 / 2.45e13
        q <- (1.5625e-8 * v1^2 + 5e-5 * v1 * v2 + v2^2) * v7 + v3^2 * v8 +
            4e-4 * v4^2 + 4e2 * v5^2 + 4e-18 * v6^2
        dq1 <- (3.125e-8 * v1 + 5e-5 * v2) * v7
        dq2 <- (5e-5 * v1 + 2 * v2) * v7
        dq3 <- v3 * v9
        dq4 <- 8e-4 * v4
        dq5 <- 800 * v5
        dq6 <- 8e-18 * v6
        factor <- 0.5 * exp(-0.5 * q)
        c(factor * dq1, factor * dq2, factor * dq3, factor * dq4, factor * dq5, factor * dq6)
    }
    eq_fn <- function(x) {
        x[1] + 4e3 * x[2] - 1.76e4
    }
    eq_jac <- function(x) {
        matrix(c(1, 4e3, 0, 0, 0, 0), nrow = 1)
    }
    lower <- c(0, -10, 0, 0, 0, 0)
    upper <- c(2e4, 10, 1e7, 20, 1, 2e8)
    start <- c(6e3, 1.5, 4e6, 2, 3e-3, 5e7)
    eq_b <- 0
    ineq_fn <- NULL
    ineq_jac <- NULL
    ineq_lower <- NULL
    ineq_upper <- NULL
    best_fn <- -exp(-27 / 280)
    best_par <- c(9.16e4 / 7, 79 / 70, 2e6, 10, 1e-3, 1e8)
    list(
        name = "hs54",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = best_fn,
        best_par = best_par
    )
}

hs55_problem <- function()
{
    # solver killer problem
    fn <- function(x) {
        x14 <- x[1] * x[4]
        if (x14 > 10) x14 <- 10
        x[1] + 2 * x[2] + 4 * x[5] + exp(x14)
    }
    gr <- function(x) {
        x14 <- x[1] * x[4]
        if (x14 > 10) {
            v1 <- exp(10)
            dfdx1 <- 1 + x[4] * v1
            dfdx4 <- x[1] * v1
        } else {
            v1 <- exp(x14)
            dfdx1 <- 1 + x[4] * v1
            dfdx4 <- x[1] * v1
        }
        c(dfdx1, 2, 0, dfdx4, 4, 0)
    }
    eq_fn <- function(x) {
        c(x[1] + 2 * x[2] + 5 * x[5] - 6, x[1] + x[2] + x[3] - 3,
            x[4] + x[5] + x[6] - 2, x[1] + x[4] - 1,
            x[2] + x[5] - 2, x[3] + x[6] - 2)
    }
    eq_jac <- function(x) {
        matrix(c(1, 2, 0, 0, 5, 0,
                 1, 1, 1, 0, 0, 0,
                 0, 0, 0, 1, 1, 1,
                 1, 0, 0, 1, 0, 0,
                 0, 1, 0, 0, 1, 0,
                 0, 0, 1, 0, 0, 1), nrow = 6, byrow = TRUE)
    }
    lower <- c(0, 0, 0, 0, 0, 0)
    upper <- c(1, 1000, 1000, 1, 1000, 1000)
    start <- c(1, 2, 0, 0, 0, 2)
    eq_b <- rep(0, 6)
    ineq_fn <- NULL
    ineq_jac <- NULL
    ineq_lower <- NULL
    ineq_upper <- NULL
    best_fn <- 19/3
    best_par <- c(0, 4/3, 5/3, 1, 2/3, 1/3)
    list(
        name = "hs55",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = best_fn,
        best_par = best_par
    )
}

hs56_problem <- function()
{
    fn <- function(x) {
        -x[1] * x[2] * x[3]
    }
    gr <- function(x) {
        c(-x[2] * x[3], -x[1] * x[3], -x[1] * x[2], 0, 0, 0, 0)
    }
    eq_fn <- function(x) {
        c(x[1] - 4.2 * sin(x[4])^2, x[2] - 4.2 * sin(x[5])^2,
          x[3] - 4.2 * sin(x[6])^2, x[1] + 2 * x[2] + 2 * x[3] - 7.2 * sin(x[7])^2)
    }
    eq_jac <- function(x) {
        matrix(c(1, 0, 0, -8.4 * sin(x[4]) * cos(x[4]), 0, 0, 0,
                 0, 1, 0, 0, -8.4 * sin(x[5]) * cos(x[5]), 0, 0,
                 0, 0, 1, 0, 0, -8.4 * sin(x[6]) * cos(x[6]), 0,
                 1, 2, 2, 0, 0, 0, -14.4 * sin(x[7]) * cos(x[7])), nrow = 4, byrow = TRUE)
    }
    lower <- rep(-1000, 7)
    upper <- rep(1000, 7)
    start <- c(
        1, 1, 1,
        asin(sqrt(1 / 4.2)),
        asin(sqrt(1 / 4.2)),
        asin(sqrt(1 / 4.2)),
        asin(sqrt(5 / 7.2))
    )
    eq_b <- rep(0, 4)
    ineq_fn <- NULL
    ineq_jac <- NULL
    ineq_lower <- NULL
    ineq_upper <- NULL
    best_fn <- -3.456
    best_par <- c(
        2.4, 1.2, 1.2,
        asin(sqrt(4 / 7)),
        asin(sqrt(2 / 7)),
        asin(sqrt(2 / 7)),
        2 * atan(1)
    )
    list(
        name = "hs56",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = best_fn,
        best_par = best_par
    )
}

hs57_problem <- function()
{
    A <- rep(0, 44); B <- rep(0, 44)

    # Step 2: Fortran style assignments (copy the block structure)
    A[1:2] <- 8; A[16+1:2] <- 18; A[30+1:2] <- 28; A[35+1:2] <- 32; A[38+1:2] <- 36; A[40+1:2] <- 38
    B[1:2] <- 0.49; B[6+1:2] <- 0.46; B[11+1:2] <- 0.43; B[14+1:2] <- 0.43; B[18+1:2] <- 0.42
    B[21+1:2] <- 0.41; B[25+1:2] <- 0.40; B[29+1:2] <- 0.41; B[36+1:2] <- 0.40; B[40+1:2] <- 0.40
    B[42+1:2] <- 0.39
    for (i in 1:3) {
        A[10+i] <- 14;  A[13+i] <- 16;  A[18+i] <- 20
        A[21+i] <- 22;  A[24+i] <- 24;  A[27+i] <- 26;  A[32+i] <- 30
        B[31+i] <- 0.40
    }
    for (i in 1:4) {
        A[2+i] <- 10; A[6+i] <- 12
    }
    A[38] <- 34; A[43] <- 40; A[44] <- 42
    B[3] <- 0.48; B[4] <- 0.47; B[5] <- 0.48; B[6] <- 0.47; B[9] <- 0.45
    B[10] <- 0.43; B[11] <- 0.45; B[14] <- 0.44; B[17] <- 0.46; B[18] <- 0.45
    B[21] <- 0.43; B[24] <- 0.40; B[25] <- 0.42; B[28] <- 0.41; B[29] <- 0.40
    B[35] <- 0.38; B[36] <- 0.41; B[39] <- 0.41; B[40] <- 0.38
    fn <- function(x) {
        F <- B - x[1] - (0.49 - x[1]) * exp(-x[2] * (A - 8))
        sum(F^2)
    }
    gr <- function(x) {
        F <- B - x[1] - (0.49 - x[1]) * exp(-x[2] * (A - 8))
        V1 <- exp(-x[2] * (A - 8))
        DF1 <- -1 + V1
        DF2 <- (A - 8) * (0.49 - x[1]) * V1
        g1 <- 2 * sum(F * DF1)
        g2 <- 2 * sum(F * DF2)
        c(g1, g2)
    }
    ineq_fn <- function(x) {
        -x[1] * x[2] + 0.49 * x[2] - 0.09
    }
    ineq_jac <- function(x) {
        matrix(c(-x[2], -x[1] + 0.49), nrow = 1)
    }
    lower <- c(0.4, -4)
    upper <- c(1000, 1000)
    start <- c(0.42, 5)
    ineq_lower <- 0
    ineq_upper <- 1e8
    eq_fn <- NULL
    eq_jac <- NULL
    eq_b <- NULL
    best_fn <- 0.0284596697213
    best_par <- c(0.419952674511, 1.28484562930)
    list(
        name = "hs57",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = best_fn,
        best_par = best_par
    )
}

hs58_problem <- function()
{
    fn <- function(x) {
        100 * (x[2] - x[1]^2)^2 + (1 - x[1])^2
    }
    gr <- function(x) {
        g2 <- 200 * (x[2] - x[1]^2)
        g1 <- -2 * (x[1] * (g2 - 1) + 1)
        c(g1, g2)
    }
    ineq_fn <- function(x) {
        c(x[2]^2 - x[1], x[1]^2 - x[2], x[1]^2 + x[2]^2 - 1)
    }
    ineq_jac <- function(x) {
        matrix(c(-1, 2 * x[2],
                 2 * x[1], -1,
                 2 * x[1], 2 * x[2]), nrow = 3, byrow = TRUE)
    }
    lower <- c(-2, -1000)
    upper <- c(0.5, 1000)
    start <- c(-2, 1)
    eq_fn <- NULL
    eq_jac <- NULL
    eq_b <- NULL
    ineq_lower <- rep(0, 3)
    ineq_upper <- rep(1e8, 3)
    best_fn <- 3.19033354957
    best_par <- c(-0.786150483331, 0.618034533851)
    list(
        name = "hs58",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = best_fn,
        best_par = best_par
    )
}

hs59_problem <- function()
{
    fn <- function(x) {
        x1 <- x[1]; x2 <- x[2]
        x11 <- x1; x12 <- x1^2; x13 <- x1^3; x14 <- x1^4
        x21 <- x2; x22 <- x2^2; x23 <- x2^3; x24 <- x2^4
        xx11 <- x1 * x2; xx12 <- x1 * x2^2; xx21 <- x1^2 * x2
        xx31 <- x1^3 * x2
        -75.196 + 3.8112 * x11 - 0.12694 * x12 + 2.0567e-3 * x13 - 1.0345e-5 * x14 +
            6.8306 * x21 - 3.0234e-2 * x11 * x21 + 1.28134e-3 * x12 * x21 -
            3.5256e-5 * x13 * x21 + 2.266e-7 * x14 * x21 -
            0.25645 * x22 + 3.4604e-3 * x23 - 1.3514e-5 * x24 +
            28.106 / (x21 + 1) + 5.2375e-6 * x12 * x22 + 6.3e-8 * x13 * x22 -
            7e-10 * x13 * x23 - 3.4054e-4 * x11 * x22 + 1.6638e-6 * x11 * x23 +
            2.8673 * exp(5e-4 * x11 * x21)
    }
    gr <- function(x) {
        x1 <- x[1]; x2 <- x[2]
        x11 <- x1; x12 <- x1^2; x13 <- x1^3; x14 <- x1^4
        x21 <- x2; x22 <- x2^2; x23 <- x2^3
        xx11 <- x1 * x2; xx12 <- x1 * x2^2; xx21 <- x1^2 * x2; xx31 <- x1^3 * x2

        expterm <- exp(5e-4 * xx11)
        g1 <- 3.8112 - 0.25388 * x1 + 6.1701e-3 * x12 - 4.138e-5 * x13 -
            3.0234e-2 * x2 + 2.56268e-3 * xx11 - 1.05768e-4 * xx21 + 9.064e-7 * xx31 +
            1.0475e-5 * xx12 + 1.89e-7 * x12 * x22 - 2.1e-9 * x12 * x23 -
            3.4054e-4 * x22 + 1.6638e-6 * x23 + 1.43365e-3 * x2 * expterm
        g2 <- 6.8306 - 3.0234e-2 * x1 + 1.28134e-3 * x12 - 3.5256e-5 * x13 +
            2.266e-7 * x14 - 0.5129 * x2 + 1.03812e-2 * x22 - 5.4056e-5 * x23 - 28.106 / (x2 + 1)^2 +
            1.0475e-5 * xx21 + 1.26e-7 * xx31 - 2.1e-9 * x13 * x22 - 6.8108e-4 * xx11 + 4.9914e-6 * xx12 +
            1.43365e-3 * x1 * expterm
        c(g1, g2)
    }
    ineq_fn <- function(x) {
        c(x[1] * x[2] - 700, x[2] - 8e-3 * x[1]^2, (x[2] - 50)^2 - 5 * (x[1] - 55))
    }
    ineq_jac <- function(x) {
        matrix(c(x[2], x[1],
                 -1.6e-2 * x[1], 1,
                 -5, 2 * (x[2] - 50)), nrow = 3, byrow = TRUE)
    }
    lower <- c(0, 0)
    upper <- c(75, 65)
    start <- c(90, 10)
    ineq_lower <- rep(0, 3)
    ineq_upper <- rep(1e8, 3)
    eq_fn <- NULL
    eq_jac <- NULL
    eq_b <- NULL
    best_fn <- -7.80422632408
    best_par <- c(13.5501042366, 51.6601812877)
    list(
        name = "hs59",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = best_fn,
        best_par = best_par
    )
}

hs60_problem <- function()
{
    fn <- function(x) {
        (x[1] - 1)^2 + (x[1] - x[2])^2 + (x[2] - x[3])^4
    }
    gr <- function(x) {
        v1 <- 2 * (x[1] - x[2])
        g1 <- 2 * (x[1] - 1) + v1
        g3 <- -4 * (x[2] - x[3])^3
        g2 <- -g3 - v1
        c(g1, g2, g3)
    }
    eq_fn <- function(x) {
        x[1] * (1 + x[2]^2) + x[3]^4 - 4 - 3 * sqrt(2)
    }
    eq_jac <- function(x) {
        matrix(c(1 + x[2]^2, 2 * x[1] * x[2], 4 * x[3]^3), nrow = 1)
    }
    lower <- rep(-10, 3)
    upper <- rep(10, 3)
    start <- rep(2, 3)
    eq_b <- 0
    ineq_fn <- NULL
    ineq_jac <- NULL
    ineq_lower <- NULL
    ineq_upper <- NULL
    best_fn <- 0.0325682002513
    best_par <- c(1.10485902423, 1.19667419413, 1.53526225739)
    list(
        name = "hs60",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = best_fn,
        best_par = best_par
    )
}

hs61_problem <- function()
{
    fn <- function(x) {
        4 * x[1]^2 + 2 * x[2]^2 + 2 * x[3]^2 - 33 * x[1] + 16 * x[2] - 24 * x[3]
    }
    gr <- function(x) {
        c(8 * x[1] - 33, 4 * x[2] + 16, 4 * x[3] - 24)
    }
    eq_fn <- function(x) {
        c(3 * x[1] - 2 * x[2]^2 - 7, 4 * x[1] - x[3]^2 - 11)
    }
    eq_jac <- function(x) {
        matrix(c(3, -4 * x[2], 0, 4, 0, -2 * x[3]), nrow = 2, byrow = TRUE)
    }
    lower <- rep(-1000, 3)
    upper <- rep(1000, 3)
    start <- c(0, 0, 0)
    eq_b <- c(0, 0)
    ineq_fn <- NULL
    ineq_jac <- NULL
    ineq_lower <- NULL
    ineq_upper <- NULL
    best_fn <- -143.646142201
    best_par <- c(5.32677015744, -2.11899863998, 3.21046423906)
    list(
        name = "hs61",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = best_fn,
        best_par = best_par
    )
}

hs62_problem <- function()
{
    fn <- function(x) {
        # Variable substitutions
        B3 <- x[3] + 0.03
        C3 <- 0.13 * x[3] + 0.03
        B2 <- B3 + x[2]
        C2 <- B3 + 0.07 * x[2]
        B1 <- B2 + x[1]
        C1 <- B2 + 0.09 * x[1]
        V5 <- B1 / C1
        V6 <- B2 / C2
        V7 <- B3 / C3

        if (V5 <= 0 || V6 <= 0 || V7 <= 0) {
            # Penalty for infeasible region
            S <- sum((x - 5)^2)
            return(S + 1e3 - 2.67e4)
        }

        -32.174 * (255 * log(V5) + 280 * log(V6) + 290 * log(V7))
    }
    gr <- function(x) {
        # Variable substitutions
        B3 <- x[3] + 0.03
        C3 <- 0.13 * x[3] + 0.03
        B2 <- B3 + x[2]
        C2 <- B3 + 0.07 * x[2]
        B1 <- B2 + x[1]
        C1 <- B2 + 0.09 * x[1]
        # Reciprocals
        RB1 <- 1 / B1
        RB2 <- 1 / B2
        RB3 <- 1 / B3
        RC1 <- 1 / C1
        RC2 <- 1 / C2
        RC3 <- 1 / C3
        V1 <- -32.174 * 255
        V2 <- -32.174 * 280
        V3 <- -32.174 * 290
        V4 <- V1 * (RB1 - RC1)

        # Gradient (as per Fortran)
        g1 <- V1 * (RB1 - 0.09 * RC1)
        g2 <- V4 + V2 * (RB2 - 0.07 * RC2)
        g3 <- V4 + V2 * (RB2 - RC2) + V3 * (RB3 - 0.13 * RC3)
        c(g1, g2, g3)
    }
    eq_fn <- function(x) {
        sum(x) - 1
    }
    eq_jac <- function(x) {
        matrix(rep(1, 3), nrow = 1)
    }
    lower <- rep(0, 3)
    upper <- rep(1, 3)
    start <- c(0.7, 0.2, 0.1)
    eq_b <- 0
    ineq_fn <- NULL
    ineq_jac <- NULL
    ineq_lower <- NULL
    ineq_upper <- NULL
    best_fn <- -26272.5144873
    best_par <- c(0.617813298210, 0.328202155786, 0.0539845460119)
    list(
        name = "hs62",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = best_fn,
        best_par = best_par
    )
}

hs63_problem <- function()
{
    fn <- function(x) {
        1e3 - x[1]^2 - 2 * x[2]^2 - x[3]^2 - x[1] * x[2] - x[1] * x[3]
    }
    gr <- function(x) {
        c(-2 * x[1] - x[2] - x[3], -4 * x[2] - x[1], -2 * x[3] - x[1])
    }
    eq_fn <- function(x) {
        c(8 * x[1] + 14 * x[2] + 7 * x[3] - 56,
          x[1]^2 + x[2]^2 + x[3]^2 - 25)
    }
    eq_jac <- function(x) {
        matrix(c(8, 14, 7, 2 * x[1], 2 * x[2], 2 * x[3]),nrow = 2, byrow = TRUE)
    }
    lower <- rep(0, 3)
    upper <- rep(1000, 3)
    start <- rep(2, 3)
    eq_b <- c(0, 0)
    ineq_fn <- NULL
    ineq_jac <- NULL
    ineq_lower <- NULL
    ineq_upper <- NULL
    best_fn <- 961.715172127
    best_par <- c(3.51211841492, 0.216988174172, 3.55217403459)
    list(
        name = "hs63",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = best_fn,
        best_par = best_par
    )
}

hs64_problem <- function()
{
    fn <- function(x) {
        5 * x[1] + 5e4 / x[1] + 20 * x[2] + 7.2e4 / x[2] + 10 * x[3] + 1.44e5 / x[3]
    }
    gr <- function(x) {
        c(5 - 5e4 / x[1]^2, 20 - 7.2e4 / x[2]^2, 10 - 1.44e5 / x[3]^2)
    }
    ineq_fn <- function(x) {
        1 - 4 / x[1] - 32 / x[2] - 120 / x[3]
    }
    ineq_jac <- function(x) {
        matrix(c(4 / x[1]^2, 32 / x[2]^2, 120 / x[3]^2), nrow = 1)
    }
    lower <- rep(1e-5, 3)
    upper <- rep(1000, 3)
    start <- rep(1, 3)
    eq_fn <- NULL
    eq_jac <- NULL
    eq_b <- NULL
    ineq_lower <- 0
    ineq_upper <- 1e8
    best_fn <- 6299.84242821
    best_par <- c(108.734717597, 85.1261394257, 204.324707858)
    list(
        name = "hs64",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = best_fn,
        best_par = best_par
    )
}

hs65_problem <- function()
{
    fn <- function(x) {
        (x[1] - x[2])^2 + ((x[1] + x[2] - 10)/3)^2 + (x[3] - 5)^2
    }
    gr <- function(x) {
        v1 <- 2 * (x[1] - x[2])
        v2 <- 2 * (x[1] + x[2] - 10) / 9
        c(v1 + v2, -v1 + v2, 2 * (x[3] - 5))
    }
    ineq_fn <- function(x) {
        48 - x[1]^2 - x[2]^2 - x[3]^2
    }
    ineq_jac <- function(x) {
        matrix(c(-2 * x[1], -2 * x[2], -2 * x[3]), nrow = 1)
    }
    lower <- c(-4.5, -4.5, -5)
    upper <- c(4.5, 4.5, 5)
    start <- c(-5, 5, 0)
    eq_fn <- NULL
    eq_jac <- NULL
    eq_b <- NULL
    ineq_lower <- 0
    ineq_upper <- 1e8
    best_fn <- 0.953528856757
    best_par <- c(3.65046182158, 3.65046168940, 4.62041750754)
    list(
        name = "hs65",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = best_fn,
        best_par = best_par
    )
}

alkylation_problem <- function() {
    fn <- function(x) {
        -0.63 * x[4] * x[7] + 50.4 * x[1] + 3.5 * x[2] + x[3] + 33.6 * x[5]
    }

    gr <- NULL  # If you want to implement gradient later

    eq_fn <- function(x) {
        z1 <- 98 * x[3] - 0.1 * x[4] * x[6] * x[9] - x[3] * x[6]
        z2 <- 1000 * x[2] + 100 * x[5] - 100 * x[1] * x[8]
        z3 <- 122 * x[4] - 100 * x[1] - 100 * x[5]
        return(c(z1, z2, z3))
    }

    eq_jac <- NULL  # Optional, if analytic Jacobian is available

    ineq_fn <- function(x) {
        z1 <- (1.12 * x[1] + 0.13167 * x[1] * x[8] - 0.00667 * x[1] * x[8]^2) / x[4]
        z2 <- (1.098 * x[8] - 0.038 * x[8]^2 + 0.325 * x[6] + 57.25) / x[7]
        z3 <- (-0.222 * x[10] + 35.82) / x[9]
        z4 <- (3 * x[7] - 133) / x[10]
        return(c(z1, z2, z3, z4))
    }

    ineq_jac <- NULL  # Optional, if analytic Jacobian is available

    lower <- c(0, 0, 0, 10, 0, 85, 10, 3, 1, 145)
    upper <- c(20, 16, 120, 50, 20, 93, 95, 12, 4, 162)
    ineq_lower <- c(0.99, 0.99, 0.9, 0.99)
    ineq_upper <- c(100 / 99, 100 / 99, 10 / 9, 100 / 99)
    eq_b <- c(0, 0, 0)
    start <- c(17.45, 12, 110, 30, 19.74, 89.2, 92.8, 8, 3.6, 155)
    best_par <- c(16.996427, 16.000000, 57.685751, 30.324940, 20.000000, 90.565147, 95.000000, 10.590461, 1.561636, 153.535354)
    best_fn <- -172.642

    return(list(
        name = "alkylation",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_jac = eq_jac,
        eq_b = eq_b,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_par = best_par,
        best_fn = best_fn
    ))
}

wright4_problem <- function() {
    # Objective function
    fn <- function(x) {
        (x[1] - 1)^2 +
            (x[1] - x[2])^2 +
            (x[2] - x[3])^3 +
            (x[3] - x[4])^4 +
            (x[4] - x[5])^4
    }

    # Gradient of the objective
    gr <- function(x) {
        d1 <- 2 * (x[1] - 1) + 2 * (x[1] - x[2])
        d2 <- -2 * (x[1] - x[2]) + 3 * (x[2] - x[3])^2
        d3 <- -3 * (x[2] - x[3])^2 + 4 * (x[3] - x[4])^3
        d4 <- -4 * (x[3] - x[4])^3 + 4 * (x[4] - x[5])^3
        d5 <- -4 * (x[4] - x[5])^3
        return(c(d1, d2, d3, d4, d5))
    }

    # Equality constraints
    eq_fn <- function(x) {
        z1 <- x[1] + x[2]^2 + x[3]^3
        z2 <- x[2] - x[3]^2 + x[4]
        z3 <- x[1] * x[5]
        return(c(z1, z2, z3))
    }

    # Jacobian of equality constraints: 3 x 5 matrix
    eq_jac <- function(x) {
        jac <- matrix(0, nrow = 3, ncol = 5)

        # dz1/dx
        jac[1, 1] <- 1
        jac[1, 2] <- 2 * x[2]
        jac[1, 3] <- 3 * x[3]^2
        jac[1, 4] <- 0
        jac[1, 5] <- 0

        # dz2/dx
        jac[2, 1] <- 0
        jac[2, 2] <- 1
        jac[2, 3] <- -2 * x[3]
        jac[2, 4] <- 1
        jac[2, 5] <- 0

        # dz3/dx
        jac[3, 1] <- x[5]
        jac[3, 2] <- 0
        jac[3, 3] <- 0
        jac[3, 4] <- 0
        jac[3, 5] <- x[1]

        return(jac)
    }

    ineq_fn <- NULL
    ineq_jac <- NULL
    ineq_lower <- NULL
    ineq_upper <- NULL

    lower <- rep(-10, 5)
    upper <- rep(10, 5)
    eq_b <- c(2 + 3 * sqrt(2), -2 + 2 * sqrt(2), 2)
    start <- c(1, 1, 1, 1, 1)
    best_par <- c(1.116635, 1.220442, 1.537785, 1.972769, 1.791096)
    best_fn <- 0.02931083

    return(list(
        name = "wright4",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_jac = eq_jac,
        eq_b = eq_b,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_par = best_par,
        best_fn = best_fn
    ))
}

wright9_problem <- function() {
    fn <- function(x) {
        10 * x[1] * x[4] -
            6 * x[3] * x[2]^2 +
            x[2] * x[1]^3 +
            9 * sin(x[5] - x[3]) +
            x[5]^4 * x[4]^2 * x[2]^3
    }

    gr <- function(x) {
        df_dx1 <- 10 * x[4] + 3 * x[2] * x[1]^2
        df_dx2 <- -12 * x[3] * x[2] + x[1]^3 + 3 * x[5]^4 * x[4]^2 * x[2]^2
        df_dx3 <- -6 * x[2]^2 - 9 * cos(x[5] - x[3])
        df_dx4 <- 10 * x[1] + 2 * x[5]^4 * x[4] * x[2]^3
        df_dx5 <- 9 * cos(x[5] - x[3]) + 4 * x[5]^3 * x[4]^2 * x[2]^3
        return(c(df_dx1, df_dx2, df_dx3, df_dx4, df_dx5))
    }

    eq_fn <- NULL
    eq_jac <- NULL

    ineq_fn <- function(x) {
        z1 <- x[1]^2 + x[2]^2 + x[3]^2 + x[4]^2 + x[5]^2
        z2 <- x[1]^2 * x[3] - x[4] * x[5]
        z3 <- x[2]^2 * x[4] + 10 * x[1] * x[5]
        return(c(z1, z2, z3))
    }

    ineq_jac <- function(x) {
        # Rows: g1, g2, g3; Columns: partials w.r.t. x1 to x5
        jac <- matrix(0, nrow = 3, ncol = 5)

        # g1
        jac[1, ] <- 2 * x

        # g2
        jac[2, 1] <- 2 * x[1] * x[3]
        jac[2, 2] <- 0
        jac[2, 3] <- x[1]^2
        jac[2, 4] <- -x[5]
        jac[2, 5] <- -x[4]

        # g3
        jac[3, 1] <- 10 * x[5]
        jac[3, 2] <- 2 * x[2] * x[4]
        jac[3, 3] <- 0
        jac[3, 4] <- x[2]^2
        jac[3, 5] <- 10 * x[1]

        return(jac)
    }

    lower <- rep(-5, 5)
    upper <- rep(5, 5)
    ineq_lower <- c(-100, -2, 5)
    ineq_upper <- c(20, 100, 100)
    eq_b <- NULL
    start <- c(1, 1, 1, 1, 1)
    best_par <- c(-0.08145219, 3.69237756, 2.48741102, 0.37713392, 0.17398257)
    best_fn <- -210.4078

    return(list(
        name = "wright9",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_jac = eq_jac,
        eq_b = eq_b,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_par = best_par,
        best_fn = best_fn
    ))
}

entropy_problem <- function() {
    # Objective function
    fn <- function(x) {
        m <- length(x)
        f <- -sum(log(x))
        f - log(sqrt(sum((x - 1)^2)) + 0.1)
    }

    # Gradient of the objective function
    gr <- function(x) {
        # Compute vector norm of (x - 1)
        diff <- x - 1
        norm_diff <- sqrt(sum(diff^2))

        # First term: -1 / x
        g1 <- -1 / x

        # Second term: -(x_i - 1) / [norm(x - 1) * (norm(x - 1) + 0.1)]
        g2 <- -diff / (norm_diff * (norm_diff + 0.1))

        return(g1 + g2)
    }

    # Equality constraint: sum(x) == 10
    eq_fn <- function(x) sum(x)

    # Jacobian of equality constraint
    eq_jac <- function(x) matrix(1, nrow = 1, ncol = length(x))

    # Inequality: none
    ineq_fn <- NULL
    ineq_jac <- NULL
    ineq_lower <- NULL
    ineq_upper <- NULL

    # Bounds and other values
    n <- 10
    lower <- rep(0, n)
    upper <- rep(1000, n)
    eq_b <- 10
    start <- runif(n, 0, 1000)
    best_fn <- 0.1854782
    best_par <- c(
        2.2801555, 0.8577605, 0.8577605, 0.8577605, 0.8577605,
        0.8577605, 0.8577605, 0.8577605, 0.8577605, 0.8577605
    )

    return(list(
        name = "entropy",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_jac = eq_jac,
        eq_b = eq_b,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_par = best_par,
        best_fn = best_fn
    ))
}

box_problem <- function() {
    # Objective function
    fn <- function(x) {
        -x[1] * x[2] * x[3]
    }

    # Gradient of the objective function
    gr <- function(x) {
        c(
            -x[2] * x[3],
            -x[1] * x[3],
            -x[1] * x[2]
        )
    }

    # Equality constraint function: 4*x1*x2 + 2*x2*x3 + 2*x3*x1 == 100
    eq_fn <- function(x) {
        4 * x[1] * x[2] + 2 * x[2] * x[3] + 2 * x[3] * x[1]
    }

    # Jacobian of equality constraint: 1 x 3 matrix
    eq_jac <- function(x) {
        matrix(c(
            4 * x[2] + 2 * x[3],
            4 * x[1] + 2 * x[3],
            2 * x[2] + 2 * x[1]
        ), nrow = 1)
    }

    # No inequality constraints
    ineq_fn <- NULL
    ineq_jac <- NULL
    ineq_lower <- NULL
    ineq_upper <- NULL

    lower <- rep(1, 3)
    upper <- rep(10, 3)
    eq_b <- 100
    start <- c(1.1, 1.1, 9)
    best_par <- c(2.886751, 2.886751, 5.773503)
    best_fn <- -48.11252

    return(list(
        name = "box",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_jac = eq_jac,
        eq_b = eq_b,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_par = best_par,
        best_fn = best_fn
    ))
}

rosen_suzuki_problem <- function() {
    # Objective function
    fn <- function(x) {
        x[1]^2 + x[2]^2 + 2 * x[3]^2 + x[4]^2 - 5 * x[1] - 5 * x[2] - 21 * x[3] + 7 * x[4]
    }

    # Gradient of the objective function
    gr <- function(x) {
        c(
            2 * x[1] - 5,
            2 * x[2] - 5,
            4 * x[3] - 21,
            2 * x[4] + 7
        )
    }

    # Inequality constraint function
    ineq_fn <- function(x) {
        z1 <- 8 - x[1]^2 - x[2]^2 - x[3]^2 - x[4]^2 - x[1] + x[2] - x[3] + x[4]
        z2 <- 10 - x[1]^2 - 2 * x[2]^2 - x[3]^2 - 2 * x[4]^2 + x[1] + x[4]
        z3 <- 5 - 2 * x[1]^2 - x[2]^2 - x[3]^2 - 2 * x[1] + x[2] + x[4]
        return(c(z1, z2, z3))
    }

    # Jacobian of the inequality constraints: 3 x 4 matrix
    ineq_jac <- function(x) {
        jac <- matrix(0, nrow = 3, ncol = 4)

        # dz1/dx
        jac[1, 1] <- -2 * x[1] - 1
        jac[1, 2] <- -2 * x[2] + 1
        jac[1, 3] <- -2 * x[3] - 1
        jac[1, 4] <- -2 * x[4] + 1

        # dz2/dx
        jac[2, 1] <- -2 * x[1] + 1
        jac[2, 2] <- -4 * x[2]
        jac[2, 3] <- -2 * x[3]
        jac[2, 4] <- -4 * x[4] + 1

        # dz3/dx
        jac[3, 1] <- -4 * x[1] - 2
        jac[3, 2] <- -2 * x[2] + 1
        jac[3, 3] <- -2 * x[3]
        jac[3, 4] <- 1

        return(jac)
    }

    eq_fn <- NULL
    eq_jac <- NULL
    eq_b <- NULL

    lower <- rep(-10, 4)
    upper <- rep(10, 4)
    ineq_lower <- rep(0, 3)
    ineq_upper <- rep(1000, 3)
    start <- c(1, 1, 1, 1)
    best_par <- c(2.502771e-07, 9.999997e-01, 2.000000e+00, -1.000000e+00)
    best_fn <- -44

    return(list(
        name = "rosen_suzuki",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_jac = eq_jac,
        eq_b = eq_b,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_par = best_par,
        best_fn = best_fn
    ))
}


powell_problem <- function() {
    # Objective function
    fn <- function(x) {
        exp(x[1] * x[2] * x[3] * x[4] * x[5])
    }

    # Gradient of the objective function
    gr <- function(x) {
        prod_term <- x[1] * x[2] * x[3] * x[4] * x[5]
        fval <- exp(prod_term)

        c(
            fval * x[2] * x[3] * x[4] * x[5],
            fval * x[1] * x[3] * x[4] * x[5],
            fval * x[1] * x[2] * x[4] * x[5],
            fval * x[1] * x[2] * x[3] * x[5],
            fval * x[1] * x[2] * x[3] * x[4]
        )
    }

    # Equality constraint function (3 constraints)
    eq_fn <- function(x) {
        z1 <- x[1]^2 + x[2]^2 + x[3]^2 + x[4]^2 + x[5]^2
        z2 <- x[2] * x[3] - 5 * x[4] * x[5]
        z3 <- x[1]^3 + x[2]^3
        return(c(z1, z2, z3))
    }

    # Jacobian of equality constraints: 3 x 5 matrix
    eq_jac <- function(x) {
        jac <- matrix(0, nrow = 3, ncol = 5)

        # dz1/dx
        jac[1, ] <- 2 * x

        # dz2/dx
        jac[2, 1] <- 0
        jac[2, 2] <- x[3]
        jac[2, 3] <- x[2]
        jac[2, 4] <- -5 * x[5]
        jac[2, 5] <- -5 * x[4]

        # dz3/dx
        jac[3, 1] <- 3 * x[1]^2
        jac[3, 2] <- 3 * x[2]^2
        jac[3, 3] <- 0
        jac[3, 4] <- 0
        jac[3, 5] <- 0

        return(jac)
    }

    ineq_fn <- NULL
    ineq_jac <- NULL
    ineq_lower <- NULL
    ineq_upper <- NULL

    lower <- rep(-10, 5)
    upper <- rep(10, 5)
    eq_b <- c(10, 0, -1)
    start <- c(-2, 2, 2, -1, -1)
    best_par <- c(-1.717144, 1.595710, 1.827245, 0.763643, 0.763643)
    best_fn <- 0.05394985

    return(list(
        name = "powell",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_jac = eq_jac,
        eq_b = eq_b,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_par = best_par,
        best_fn = best_fn
    ))
}

himmelblau5_problem <- function() {
    # Objective function
    fn <- function(x) {
        (x[1] - 5)^2 + (x[2] - 5)^2
    }

    # Gradient of objective
    gr <- function(x) {
        c(2 * (x[1] - 5), 2 * (x[2] - 5))
    }

    # Equality constraint: x1 - 3 * x2 = 0
    eq_fn <- function(x) {
        x[1] - 3 * x[2]
    }

    eq_jac <- function(x) {
        matrix(c(1, -3), nrow = 1)
    }

    # Inequality constraint: x1^2 + x2^2 <= 25
    ineq_fn <- function(x) {
        x[1]^2 + x[2]^2
    }

    ineq_jac <- function(x) {
        matrix(c(2 * x[1], 2 * x[2]), nrow = 1)
    }

    lower <- rep(-5, 2)
    upper <- rep(5, 2)
    eq_b <- 0
    ineq_lower <- -1000
    ineq_upper <- 25
    start <- c(1, 1)
    best_par <- c(4.74342, 1.58114)
    best_fn <- 11.7544

    return(list(
        name = "himmelblau5",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_jac = eq_jac,
        eq_b = eq_b,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_par = best_par,
        best_fn = best_fn
    ))
}


hs118_problem <- function() {
    # Objective function
    fn <- function(x) {
        obj <- 0
        for (i in 0:4) {
            xi <- x[(3 * i + 1):(3 * i + 3)]
            obj <- obj + (2.3 * xi[1] + 0.0001 * xi[1]^2 +
                              1.7 * xi[2] + 0.0001 * xi[2]^2 +
                              2.2 * xi[3] + 0.00015 * xi[3]^2)
        }
        return(obj)
    }

    gr <- function(x) {
        g <- numeric(15)
        for (k in 0:4) {
            i <- 3 * k + 1
            g[i]     <- 2.3 + 0.0002 * x[i]
            g[i + 1] <- 1.7 + 0.0002 * x[i + 1]
            g[i + 2] <- 2.2 + 0.0003 * x[i + 2]
        }
        return(g)
    }

    # Inequality constraints
    ineq_fn <- function(x) {
        c(
            # Temporal order constraints (lower and upper bounds handled separately)
            x[4]  - x[1]  + 7,
            x[7]  - x[4]  + 7,
            x[10] - x[7]  + 7,
            x[13] - x[10] + 7,

            x[5]  - x[2]  + 7,
            x[8]  - x[5]  + 7,
            x[11] - x[8]  + 7,
            x[14] - x[11] + 7,

            x[6]  - x[3]  + 7,
            x[9]  - x[6]  + 7,
            x[12] - x[9]  + 7,
            x[15] - x[12] + 7,

            # Demand constraints (LHS ≥ threshold)
            60 - (x[1] + x[2] + x[3]),
            50 - (x[4] + x[5] + x[6]),
            70 - (x[7] + x[8] + x[9]),
            85 - (x[10] + x[11] + x[12]),
            100 - (x[13] + x[14] + x[15])
        )
    }

    ineq_jac <- function(x) {
        J <- matrix(0, nrow = 17, ncol = 15)

        # Time progression constraints (12 rows)
        idx_pairs <- list(
            c(1, 4), c(4, 7), c(7, 10), c(10, 13),
            c(2, 5), c(5, 8), c(8, 11), c(11, 14),
            c(3, 6), c(6, 9), c(9, 12), c(12, 15)
        )

        for (i in seq_along(idx_pairs)) {
            idx1 <- idx_pairs[[i]][1]
            idx2 <- idx_pairs[[i]][2]
            J[i, idx2] <- 1
            J[i, idx1] <- -1
        }

        # Demand constraints (5 rows)
        for (k in 0:4) {
            row <- 12 + k + 1
            cols <- 3 * k + 1:3
            J[row, cols] <- -1
        }

        return(J)
    }

    # Bounds
    lower <- c(8, 43, 3,
               rep(0, 12))
    upper <- c(21, 57, 16,
               rep(c(90, 120, 60), 4))

    # Initial values
    start <- c(20, 55, 15,
               20, 60, 20,
               20, 60, 20,
               20, 60, 20,
               20, 60, 20)

    # Known best solution
    best_fn <- 664.82045
    best_par <- NULL  # Not given exactly — optional

    return(list(
        name = "hs118",
        fn = fn,
        gr = gr,
        eq_fn = NULL,
        eq_jac = NULL,
        eq_b = NULL,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = c(rep(0, 12), rep(-1000, 5)),  # >= 0 for time constraints, >= threshold for demand
        ineq_upper = c(rep(13, 12), rep(0, 5)),    # <= 13 for time constraints, demand constraint expressed as ≤ 0
        lower = lower,
        upper = upper,
        start = start,
        best_par = best_par,
        best_fn = best_fn
    ))
}


hs119_problem <- function() {
    # Constants
    c_vals <- c(2.5, 1.1, -3.1, -3.5, 1.3, 2.1, 2.3, -1.5)

    # Objective: sum of 46 nonlinear terms
    fn <- function(x) {
        square_expr <- function(z) (z^2 + z + 1)
        t <- numeric(46)
        t[1] <- square_expr(x[1]) * square_expr(x[1])
        t[2] <- square_expr(x[1]) * square_expr(x[4])
        t[3] <- square_expr(x[1]) * square_expr(x[7])
        t[4] <- square_expr(x[1]) * square_expr(x[8])
        t[5] <- square_expr(x[1]) * square_expr(x[16])
        t[6] <- square_expr(x[2]) * square_expr(x[2])
        t[7] <- square_expr(x[2]) * square_expr(x[3])
        t[8] <- square_expr(x[2]) * square_expr(x[7])
        t[9] <- square_expr(x[2]) * square_expr(x[10])
        t[10] <- square_expr(x[3]) * square_expr(x[3])
        t[11] <- square_expr(x[3]) * square_expr(x[7])
        t[12] <- square_expr(x[3]) * square_expr(x[9])
        t[13] <- square_expr(x[3]) * square_expr(x[10])
        t[14] <- square_expr(x[3]) * square_expr(x[14])
        t[15] <- square_expr(x[4]) * square_expr(x[4])
        t[16] <- square_expr(x[4]) * square_expr(x[7])
        t[17] <- square_expr(x[4]) * square_expr(x[11])
        t[18] <- square_expr(x[4]) * square_expr(x[15])
        t[19] <- square_expr(x[5]) * square_expr(x[5])
        t[20] <- square_expr(x[5]) * square_expr(x[6])
        t[21] <- square_expr(x[5]) * square_expr(x[10])
        t[22] <- square_expr(x[5]) * square_expr(x[12])
        t[23] <- square_expr(x[5]) * square_expr(x[16])
        t[24] <- square_expr(x[6]) * square_expr(x[6])
        t[25] <- square_expr(x[6]) * square_expr(x[8])
        t[26] <- square_expr(x[6]) * square_expr(x[15])
        t[27] <- square_expr(x[7]) * square_expr(x[7])
        t[28] <- square_expr(x[7]) * square_expr(x[11])
        t[29] <- square_expr(x[7]) * square_expr(x[13])
        t[30] <- square_expr(x[8]) * square_expr(x[8])
        t[31] <- square_expr(x[8]) * square_expr(x[10])
        t[32] <- square_expr(x[8]) * square_expr(x[15])
        t[33] <- square_expr(x[9]) * square_expr(x[9])
        t[34] <- square_expr(x[9]) * square_expr(x[12])
        t[35] <- square_expr(x[9]) * square_expr(x[16])
        t[36] <- square_expr(x[10]) * square_expr(x[10])
        t[37] <- square_expr(x[10]) * square_expr(x[14])
        t[38] <- square_expr(x[11]) * square_expr(x[11])
        t[39] <- square_expr(x[11]) * square_expr(x[13])
        t[40] <- square_expr(x[11]) * square_expr(x[12])
        t[41] <- square_expr(x[12]) * square_expr(x[14])
        t[42] <- square_expr(x[13]) * square_expr(x[13])
        t[43] <- square_expr(x[13]) * square_expr(x[14])
        t[44] <- square_expr(x[14]) * square_expr(x[14])
        t[45] <- square_expr(x[15]) * square_expr(x[15])
        t[46] <- square_expr(x[16]) * square_expr(x[16])
        return(sum(t))
    }

    gr <- function(x) {
        fi  <- x^2 + x + 1          # f_i(x)
        dfi <- 2*x + 1              # f'_i(x)

        ## list the 46 (a,b) index pairs exactly as used in `fn`
        idx <- matrix(c(
            1,1,  1,4, 1,7, 1,8, 1,16,
            2,2,  2,3, 2,7, 2,10,
            3,3,  3,7, 3,9, 3,10, 3,14,
            4,4,  4,7, 4,11, 4,15,
            5,5,  5,6, 5,10, 5,12, 5,16,
            6,6,  6,8, 6,15,
            7,7,  7,11, 7,13,
            8,8,  8,10, 8,15,
            9,9,  9,12, 9,16,
            10,10, 10,14,
            11,11, 11,13, 11,12,
            12,14,
            13,13, 13,14,
            14,14,
            15,15,
            16,16
        ), byrow = TRUE, ncol = 2)

        g <- numeric(16)

        for (k in seq_len(nrow(idx))) {
            i <- idx[k, 1]
            j <- idx[k, 2]

            if (i == j) {
                ## diagonal term f_i^2
                g[i] <- g[i] + 2 * fi[i] * dfi[i]
            } else {
                ## off-diagonal product f_i * f_j
                g[i] <- g[i] + fi[j] * dfi[i]
                g[j] <- g[j] + fi[i] * dfi[j]
            }
        }
        g
    }

    # Equality constraints: s[1:8] = c[1:8]
    eq_fn <- function(x) {
        c(
            0.22*x[1] + 0.2*x[2] + 0.19*x[3] + 0.25*x[4] + 0.15*x[5] + 0.11*x[6] + 0.12*x[7] + 0.13*x[8] + x[9],
            -1.46*x[1] -1.3*x[3] + 1.82*x[4] -1.15*x[5] + 0.8*x[7] + x[10],
            1.29*x[1] -0.89*x[2] -1.16*x[5] -0.96*x[6] -0.49*x[8] + x[11],
            -1.1*x[1] -1.06*x[2] + 0.95*x[3] -0.54*x[4] -1.78*x[6] -0.41*x[7] + x[12],
            -1.43*x[4] + 1.51*x[5] + 0.59*x[6] -0.33*x[7] -0.43*x[8] + x[13],
            -1.72*x[2] -0.33*x[3] + 1.62*x[5] + 1.24*x[6] + 0.21*x[7] -0.26*x[8] + x[14],
            1.12*x[1] + 0.31*x[4] + 1.12*x[7] -0.36*x[9] + x[15],
            0.45*x[2] + 0.26*x[3] -1.1*x[4] + 0.58*x[5] -1.03*x[7] + 0.1*x[8] + x[16]
        )
    }

    # Jacobian of the equality constraints (8x16)
    eq_jac <- function(x) {
        J <- matrix(0, nrow = 8, ncol = 16)
        # Fill based on coefficients from model
        J[1, 1:9] <- c(0.22, 0.2, 0.19, 0.25, 0.15, 0.11, 0.12, 0.13, 1)
        J[2, c(1, 3, 4, 5, 7, 10)] <- c(-1.46, -1.3, 1.82, -1.15, 0.8, 1)
        J[3, c(1, 2, 5, 6, 8, 11)] <- c(1.29, -0.89, -1.16, -0.96, -0.49, 1)
        J[4, c(1, 2, 3, 4, 6, 7, 12)] <- c(-1.1, -1.06, 0.95, -0.54, -1.78, -0.41, 1)
        J[5, c(4, 5, 6, 7, 8, 13)] <- c(-1.43, 1.51, 0.59, -0.33, -0.43, 1)
        J[6, c(2, 3, 5, 6, 7, 8, 14)] <- c(-1.72, -0.33, 1.62, 1.24, 0.21, -0.26, 1)
        J[7, c(1, 4, 7, 9, 15)] <- c(1.12, 0.31, 1.12, -0.36, 1)
        J[8, c(2, 3, 4, 5, 7, 8, 16)] <- c(0.45, 0.26, -1.1, 0.58, -1.03, 0.1, 1)
        return(J)
    }

    lower <- rep(0, 16)
    upper <- rep(5, 16)
    start <- rep(2, 16)
    eq_b <- c(2.5, 1.1, -3.1, -3.5, 1.3, 2.1, 2.3, -1.5)

    best_par <- c(
        0.03984735, 0.7919832, 0.2028703, 0.8443579,
        1.126991, 0.9347387, 1.681962, 0.1553009,
        1.567870, 0, 0, 0,
        0.6602041, 0, 0.6742559, 0
    )
    best_fn <- 244.899698

    return(list(
        name = "hs119",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_jac = eq_jac,
        eq_b = eq_b,
        ineq_fn = NULL,
        ineq_jac = NULL,
        ineq_lower = NULL,
        ineq_upper = NULL,
        lower = lower,
        upper = upper,
        start = start,
        best_par = best_par,
        best_fn = best_fn
    ))
}



hs110_problem <- function()
{
    # Objective function: sum of obj[1:11]
    fn <- function(x) {
        # x is length 10
        # Intermediates
        prod_vec <- numeric(11)
        prod_vec[1] <- 1
        for (i in 2:11) {
            prod_vec[i] <- x[i-1] * prod_vec[i-1]
        }
        obj1_10 <- (log(x - 2))^2 + (log(10 - x))^2
        obj11 <- -prod_vec[11]^0.2
        sum(obj1_10) + obj11
    }
    # Analytic gradient
    gr <- function(x) {
        grad <- numeric(10)
        # Compute intermediates for obj11
        prod_vec <- numeric(11)
        prod_vec[1] <- 1
        for (i in 2:11) {
            prod_vec[i] <- x[i-1] * prod_vec[i-1]
        }
        log_term1 <- log(x - 2)
        log_term2 <- log(10 - x)
        # Derivatives of obj[1:10]
        grad_obj <- 2 * log_term1 / (x - 2) - 2 * log_term2 / (10 - x)
        # Derivative of obj11
        for (j in 1:10) {
            if (prod_vec[11] > 0) {
                dprodj <- prod_vec[11] / x[j]
                grad[j] <- grad_obj[j] + (-0.2) * prod_vec[11]^(-0.8) * dprodj
            } else {
                # Outside feasible region, set NaN
                grad[j] <- NaN
            }
        }
        grad
    }
    lower <- rep(2.001, 10)
    upper <- rep(9.999, 10)
    start <- rep(9, 10)
    eq_fn <- NULL
    eq_jac <- NULL
    ineq_fn <- NULL
    ineq_jac <- NULL
    eq_b <- NULL
    ineq_lower <- NULL
    ineq_upper <- NULL
    list(
        name = "hs110",
        fn = fn,
        gr = gr,
        eq_fn = eq_fn,
        eq_b = eq_b,
        eq_jac = eq_jac,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = -45.77846971,
        best_par = NULL
    )
}

garch_problem <- function()
{
    n <- 1500
    set.seed(100)
    z <- rnorm(n)
    x <- numeric(n)
    sigma2 <- numeric(n)
    x[1] <- 0.1315172
    sigma2[1] <- 0.2211
    mu <- -0.006184353
    omega <- 0.010760430
    alpha <- 0.153408326
    beta <- 0.805877422
    for (t in 2:n) {
        sigma2[t] <- omega + alpha * (x[t - 1] - mu)^2 + beta * sigma2[t - 1]
        x[t] <- mu + z[t] * sqrt(sigma2[t])
    }

    returns <- x
    garch11_negloglik <- function(par) {
        # Parameterization: par = c(mu, omega, alpha, beta)
        mu    <- par[1]
        omega <- par[2]
        alpha <- par[3]
        beta  <- par[4]

        T <- length(returns)
        eps <- returns - mu
        sig2 <- numeric(T)

        # Initialize sigma^2 with unconditional variance
        sig2[1] <- mean(eps^2)

        for (t in 2:T) {
            sig2[t] <- omega + alpha * eps[t - 1]^2 + beta * sig2[t - 1]
        }

        nll <- 0.5 * sum(log(2*pi) + log(sig2) + eps^2 / sig2)
        return(nll)
    }

    garch11_grad <- function(par) {
        mu    <- par[1]
        omega <- par[2]
        alpha <- par[3]
        beta  <- par[4]
        T <- length(returns)
        eps <- returns - mu
        sig2 <- numeric(T)
        dsig2_domega <- numeric(T)
        dsig2_dalpha <- numeric(T)
        dsig2_dbeta  <- numeric(T)
        dsig2_dmu    <- numeric(T)

        # Initialize with sample variance of residuals
        sig2[1] <- mean(eps^2)
        dsig2_domega[1] <- 0
        dsig2_dalpha[1] <- 0
        dsig2_dbeta[1]  <- 0
        dsig2_dmu[1]    <- -2 * mean(eps)

        # Recursion
        for (t in 2:T) {
            dsig2_domega[t] <- 1 + beta * dsig2_domega[t - 1]
            dsig2_dalpha[t] <- eps[t - 1]^2 + beta * dsig2_dalpha[t - 1]
            dsig2_dbeta[t]  <- sig2[t - 1] + beta * dsig2_dbeta[t - 1]
            dsig2_dmu[t]    <- -2 * alpha * eps[t - 1] + beta * dsig2_dmu[t - 1]
            sig2[t] <- omega + alpha * eps[t - 1]^2 + beta * sig2[t - 1]
        }

        # Gradients
        dmu <- -sum(eps / sig2) + sum(0.5 * (1 / sig2 - (eps^2) / (sig2^2)) * dsig2_dmu)
        domega <- sum((1 / (2 * sig2) - eps^2 / (2 * sig2^2)) * dsig2_domega)
        dalpha <- sum((1 / (2 * sig2) - eps^2 / (2 * sig2^2)) * dsig2_dalpha)
        dbeta  <- sum((1 / (2 * sig2) - eps^2 / (2 * sig2^2)) * dsig2_dbeta)

        return(c(dmu, domega, dalpha, dbeta))
    }
    garch11_ineq_fn <- function(par) {
        # par = (mu, omega, alpha, beta)
        alpha <- par[3]
        beta  <- par[4]
        return(alpha + beta - 1)
    }

    garch11_ineq_jac <- function(par) {
        # Jacobian: returns a 1x4 row vector
        jac <- numeric(4)
        jac[3] <- 1  # derivative wrt alpha
        jac[4] <- 1  # derivative wrt beta
        return(matrix(jac, nrow = 1))
    }

    lower <- c(-1, 1e-12, 1e-12, 1e-12)
    upper <- c(1, 2, 1, 1)
    start <- c(mu = 0, omega = 0.1, alpha = 0.05, beta = 0.9)
    eq_fn <- NULL
    eq_jac <- NULL
    ineq_fn <- garch11_ineq_fn
    ineq_jac <- garch11_ineq_jac
    eq_b <- NULL
    ineq_lower <- -1
    ineq_upper <- 0
    list(
        name = "garch",
        fn = garch11_negloglik,
        gr = garch11_grad,
        eq_fn = NULL,
        eq_b = NULL,
        eq_jac = NULL,
        ineq_fn = ineq_fn,
        ineq_jac = ineq_jac,
        ineq_lower = ineq_lower,
        ineq_upper = ineq_upper,
        lower = lower,
        upper = upper,
        start = start,
        best_fn = 1074.36,
        best_par = c(-0.006184353, 0.010760430, 0.153408326, 0.805877422)
    )
}

