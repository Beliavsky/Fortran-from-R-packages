# Reference generation

The permanent tests use fixed values generated independently from the Fortran
implementation.

## Day counts

The interval 2011-08-31 to 2012-02-29 was evaluated separately for each of the
sixteen source conventions, including leap-year and month-end rules. Accrued
interest was then computed from a 5.25% coupon and redemption value 10,000.
The Business/252 count was checked against the exact extracted
`NonBusDays.Brazil.rda` serial-day table.

## Regular bond

A 5% semiannual bond issued 2020-01-15 and maturing 2030-01-15 was valued on
2024-04-15 at a 4% yield. Direct discounted-cash-flow references are:

- dirty price: 106.33533492789661
- accrued interest: 1.25
- clean price: 105.08533492789661

## Odd-coupon bond

The package documentation example with issue date 2013-11-30, first payment
2015-02-28, last regular payment 2020-02-29, and maturity 2021-04-21 was
reconstructed independently. At settlement 2014-10-15 and 5% yield:

- first coupon: 6.555248618784531
- final coupon: 5.991847826086961
- accrued interest: 4.582872928176794
- dirty price: 105.81073224263731
- clean price: 101.22785931446052
- modified duration in coupon periods: 10.534694027524306
- Macaulay duration in coupon periods: 10.798061378212413
- convexity in coupon periods: 68.40044885201478

## Low-level derivatives

Fixed analytical references for `dm_MyPriceEqn`, `ModDUR`, and `CONV` were
computed directly from the source equations using independent gamma-function
and scalar arithmetic.
