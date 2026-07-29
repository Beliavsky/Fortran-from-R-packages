# API reference

All real-valued arguments and results use `real(dp)`, where
`dp = kind(1.0d0)`. Array arguments are assumed-shape arrays. Each function has
an optional final integer `istat` output.

## Valuation

| Original R export | Fortran signature before optional `istat` |
|---|---|
| `disCouponPmtsBond` | `discouponpmtsbond(couponpmt(:), times(:), r)` |
| `disMaturityValBond` | `dismaturityvalbond(bondmaturityval, n, r)` |
| `bondPriceYearlyCoupons` | `bondpriceyearlycoupons(couponpmt(:), times(:), bondmaturityval, n, r)` |
| `pvCouponDeficiency` | `pvcoupondeficiency(coupondeficiency(:), times(:), r)` |
| `bondPriceDefCoupon` | `bondpricedefcoupon(parvalue, coupondeficiency(:), times(:), r)` |
| `pvExcessCoupon` | `pvexcesscoupon(couponexcess(:), times(:), r)` |
| `bondPriceExcessCoupon` | `bondpriceexcesscoupon(couponexcess(:), times(:), r)` |
| `pricingZeroCouponBond` | `pricingzerocouponbond(maturityval, n, r)` |
| `pricingSaCpnBond` | `pricingsacpnbond(sacoupons(:), times(:), maturityval, n, r)` |
| `pricingQtrlyCpnBond` | `pricingqtrlycpnbond(qcoupons(:), times(:), mv, n, r)` |
| `pricingWithSpots` | `pricingwithspots(coupons(:), spots(:), times(:), mv, n)` |
| `pricingWithSptSeq` | `pricingwithsptseq(cpns(:), sp(:), t(:), mv, n)` |
| `matrixMethod` | `matrixmethod(couponpmt(:), times(:), maturityval, n, r1, r2)` |
| `pricingFRN` | `pricingfrn(estrtrn(:), t(:), mv, maturityperiod, estdisc)` |
| `frPricing` | `frpricing(cpns(:), fri(:), mv, n)` |
| `pricingWithZspread` | `pricingwithzspread(cpns(:), spots(:), t(:), mv, n, zsprd)` |
| `pricingWithGspread` | `pricingwithgspread(coupons(:), t(:), mv, n, ytmbenchgovtbond, gspread)` |

## Yields, spreads, and curves

| Original R export | Fortran signature before optional `istat` |
|---|---|
| `ytmZeroCouponBond` | `ytmzerocouponbond(maturityval, n, price)` |
| `computingBondYtmRateFiveDecimalPlaces` | `computingbondytmratefivedecimalplaces(couponpmt, mv, bondpv, period)` |
| `computingBondYtmRateSixDecimalPlaces` | `computingbondytmratesixdecimalplaces(couponpmt, mv, bondpv, period)` |
| `convertAPRtoDifferentPeriodcity` | `convertaprtodifferentperiodcity(givenapr, givenperiodicity, desiredperiodicity)` |
| `extraCompensationForHigherRisk` | `extracompensationforhigherrisk(aprofriskybond, aprofcomparablebond)` |
| `annualYtmZcbForPeriodicity` | `annualytmzcbforperiodicity(maturityval, yearstomaturity, zcbprice, desiredperiodicity)` |
| `earZcbVariousPeriodicity` | `earzcbvariousperiodicity(maturityval, yearstomaturity, zcbprice, desiredperiodicity)` |
| `returnIncomeFRN` | `returnincomefrn(index, qtdmargin, maturityval, periodicity)` |
| `periodicDiscRateFRN` | `periodicdiscratefrn(estrtrn, mvfrn, pricefrn, maturityyears, periodicity)` |
| `discMarginFRN` | `discmarginfrn(index, estrtrn, mvfrn, pricefrn, maturityyears, periodicity)` |
| `computingYTC` | `computingytc(couponpmt, callableval, bondpv, maturityyears, ytcyears)` |
| `computingParRate` | `computingparrate(spotrates(:), times(:), mv, pv, n)` |
| `forwards` | `forwards(spots(:), yrsfrbegins, yrsfrapplies, t(:), n)` |
| `saForwards` | `saforwards(spots(:), bgn, aply, times(:), n)` |
| `computingGspread` | `computinggspread(ytmcorpbond, ytmbenchgovtbond)` |
| `computingZspread` | `computingzspread(coupons, mv, bondpv, n, spots(:))` |

Added conventional helper:

- `effectiveannualratezcb(maturityval, yearstomaturity, zcbprice)`

## Money-market instruments

| Original R export | Fortran signature before optional `istat` |
|---|---|
| `pricingTbill` | `pricingtbill(maturityval, daystomaturity, daysinyear, mmquoteddiscrate)` |
| `pricingMoneyMarketInstrUsingAOR` | `pricingmoneymarketinstrusingaor(maturityval, daystomaturity, daysinyear, aor)` |
| `fvMoneyMarketInstrUsingAOR` | `fvmoneymarketinstrusingaor(pvmmi, daystomaturity, daysinyear, aor)` |
| `computingAORMoneyMarketInstr` | `computingaormoneymarketinstr(pvmmi, fvmmi, daystomaturity, daysinyear)` |
| `pricingCommercialPaper` | `pricingcommercialpaper(maturityval, daystomaturity, daysinyear, mmquoteddiscrate)` |
| `computingQuotedDiscRateMMI` | `computingquoteddiscratemmi(pvmmi, fvmmi, daystomaturity, daysinyear)` |
| `fvMmiUsingQuotedDiscRate` | `fvmmiusingquoteddiscrate(pvmmi, daystomaturity, daysinyear, mmquoteddiscrate)` |

The R functions that format rounded whole-dollar values as character strings
return numeric `real(dp)` values in Fortran.

## Accrued interest and duration

| Original R export | Fortran signature before optional `istat` |
|---|---|
| `aiActDtCon` | `aiactdtcon(cpmt, dt1, dt2, stdt)` |
| `aiRoundedDaysConv` | `airoundeddaysconv(cpmt, bfrstldt, stldt, elpsmnths, daysbtwncpns)` |
| `macDuration` | `macduration(n, ytm, coupon, maturityval, dayscpntosettle, dayscouponperiod)` |
| `pvFullPrice` | `pvfullprice(n, ytm, coupon, maturityval, dayscpntosettle, dayscouponperiod)` |
| `macDurationOnFP` | `macdurationonfp(fp, n, ytm, cpn, mv, dayscpntosettle, dayscouponperiod)` |
| `macDurationOnCouponRate` | `macdurationoncouponrate(couponrate, n, ytm, tt)` |
| `modifDuration` | `modifduration(n, ytm, coupon, maturityval, dayscpntosettle, dayscouponperiod)` |
| `modifDurationUsingMacDuration` | `modifdurationusingmacduration(macduration_value, ytm)` |
| `estimatedPercentChangePVFullPrice` | `estimatedpercentchangepvfullprice(annualmodifduration, changeinannualytm)` |
| `approxModifDuration` | `approxmodifduration(pvbase, pvplus, pvminus, percentchangeytm)` |
| `approxMacDurationUsingApprModifDuration` | `approxmacdurationusingapprmodifduration(approxmodifduration_value, periodicytm)` |
| `effDurtnCallableBond` | `effdurtncallablebond(pvbase, pvplus, pvminus, perchangebenchytm)` |
| `moneyDuration` | `moneyduration(macduration_value, ytm, pvfullbondprice)` |
| `changePvFullBondPrice` | `changepvfullbondprice(moneyduration_value, changeytm)` |
| `computingBondPVBP` | `computingbondpvbp(pvplus, pvminus)` |

Date arguments in `aiactdtcon` and `airoundeddaysconv` are integer serial-day
values. Only date differences are used, so any consistent day origin is valid.

Added conventional helper:

- `conventionalpercentpricechange(modified_duration, yield_change)`

## Status constants

```fortran
ba_success          = 0
ba_invalid_argument = 1
ba_size_mismatch    = 2
ba_no_root          = 3
ba_out_of_range     = 4
```

On failure, numerical procedures return an IEEE quiet NaN.
