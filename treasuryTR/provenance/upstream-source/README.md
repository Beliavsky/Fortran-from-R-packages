# treasuryTR

R Package for generating Total Returns (TR) from bond yield data with fixed maturity, e.g. reported treasury yields.

## Intro

While Treasury yields are easy to come by (see [FRED](https://fred.stlouisfed.org/)), total return (TR) indices are not. The latter is earned by investors, and is therefore of paramount importance e.g. when simulating a treasury-stock diversified portfolio. A supplier for proprietary TR Treasury index data is [CRSP](https://www.crsp.org/). Their data can be purchased or accesses trough a handful of commercial research platforms.

Swinkels (2019) compute returns from publicly available yield-to-maturity data using *standard (fixed-income) textbook formulas*. See also the following post on quant.stackexchange: https://quant.stackexchange.com/a/57403

## Installation

From CRAN:

```
install.packages("treasuryTR")
```

From Github:

```
# install.packages("devtools")
devtools::install_github("mgei/treasuryTR")
```

## Use

Vignette: https://cran.r-project.org/web/packages/treasuryTR/vignettes/treasuryTR.html

## References

Swinkels, Laurens. 2019. “Treasury Bond Return Data Starting in 1962.” Data 4 (3): 91.


