## more sensible to reparametersise as baserate, timehr? then trt eff on baserate gives a PH model
## or is trt eff only in first period a sensible model?
## either way might want utilities to describe various kinds of fitted curves 
hpexp <- function(t, rate, knots){
    check_pexp(rate, knots)
    knots <- sort(knots)
    rate[findInterval(t, knots) + 1]
}

## TODO check t is valid.  return zero for zero, other special values? 
Hpexp <- function(t, rate, knots){
    check_pexp(rate, knots)
    k0 <- c(0,knots)
    int <- findInterval(t, k0)
    dk <- diff(k0)
    cumrate <- c(0, cumsum(dk*rate[-length(rate)]))    # cumulative rate up to each knot
    cumrate[int] + rate[int]*(t - k0[int])             # add on the remainder for each t
}

## TODO handle special values for rate 
check_pexp <- function(rate, knots){
    if (!is.numeric(rate)) stop("`rate` should be numeric")
    if (!is.numeric(knots)) stop("`knots` should be numeric")
    if (any(knots <= 0)) stop("`knots` should all be positive")
    nr <- length(rate)
    nk <- length(knots)
    if (nr != nk + 1)
        stop(sprintf("found nr=%s rates and nk=%s knots, should have nr = nk+1", nr, nk))
}

hpw <- function(knots){
    nk <- length(knots)
    unroll.function(hpexp, rate=0:nk)
}

Hpw <- function(knots){
    nk <- length(knots)
    unroll.function(Hpexp, rate=0:nk)
}

custom_pw <- function(knots){
    nk <- length(knots)
    list(name = "pw", 
         pars = paste0("rate", 0:nk),
         location = "rate0",
         transforms = rep(list(log), nk+1),
         inv.transforms = rep(list(exp), nk+1),
         inits =  function(t,aux) {lam <- sr.exp.inits(t,aux); rep(lam, nk+1)})
}

## todo d, p functions?  q? r even? copy from msm? 
