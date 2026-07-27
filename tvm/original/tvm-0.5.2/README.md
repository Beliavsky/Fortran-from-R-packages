#  Package Description

The tvm package aims to simplify financial calculations, involving loan payments and rates, and the transformation from discount factors to different rate types.

It has two sections.

The first one deals with fixed equal payment loans. There you have functions similar to PMT, RATE, etc from Excel.

The second one deals with rate curves and different rates for different loan structures (zero coupon, bullet, french, german, etc).

# Quick Examples

```
library(tvm)
# Present values and internal rate of return calculations
npv(i = 0.01, cf = c(-1, 0.5, 0.9), ts = c(0, 1, 3))
xnpv(i = 0.01, cf = c(-1, 0.5, 0.9), d = as.Date(c("2015-01-01", "2015-02-15", "2015-04-10")))
irr(cf = c(-1, 0.5, 0.9), ts = c(0, 1, 3))
xirr(cf = c(-1, 1.5), d = Sys.Date() + c(0, 365))
# Typical loan calculations
pmt(amt = 100, maturity = 10, rate = 0.05)
rate(amt = 100, maturity = 10, pmt = 15)
loan(rate = 0.05, maturity = 10, amt = 100, type = "bullet")
# Get the cashflow for a loan
l <- loan(rate = 0.05, maturity = 10, amt = 100, type = "bullet")
cashflow(l)
# Build a rate curve from different inputs
rate_curve(rates = c(0.1, 0.2, 0.3), rate_type = "zero_eff")
rate_curve(fun_r = function(x) rep_len(0.1, length(x)), rate_type = "swap", knots = 1:10)
rate_curve(fun_d = function(x) 1 / (1 + x), knots = 1:10)
# Subset a rate curve, maybe transforming it to another rate type
r <- rate_curve(rates = c(0.1, 0.2, 0.3), rate_type = "zero_eff")
r["zero_eff"]
r["swap",c(1.5, 2)]
# Plot a rate curve
plot(r)
plot(r, rate_type = "german")
plot(r, rate_type = c("french", "german"))
```

# Installation instructions

`tvm` lives on CRAN, so installation is easy with `install.packages("tvm"")`

# More Details

Please read the introductory vignette