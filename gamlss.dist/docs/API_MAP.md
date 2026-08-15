# API map

The main entry point is:

```fortran
use gamlss_dist
```

## Distribution modules

### v0.1.0 families retained

| R family | Fortran routines | Module |
|---|---|---|
| NO/NO2 | d/p/q/rNO* | gamlss_continuous |
| EXP/GA | d/p/q/rEXP, d/p/q/rGA | gamlss_continuous |
| WEI/WEI2/WEI3 | d/p/q/rWEI* | gamlss_continuous |
| BE/BEo | d/p/q/rBE* | gamlss_continuous |
| BEINF/BEINF0/BEINF1 | d/p/q/rBEINF* | gamlss_continuous |
| LO/GU | d/p/q/rLO, d/p/q/rGU | gamlss_continuous |
| LOGNO/LOGNO2/LNO | d/p/q/rLOGNO*, d/p/q/rLNO | gamlss_continuous |
| TF | dTF, pTF, qTF, rTF | gamlss_continuous |
| PE/PE2 | d/p/q/rPE* | gamlss_continuous |
| EGB2 | d/p/q/rEGB2 | gamlss_continuous |
| GB1/GB2 | d/p/q/rGB1, d/p/q/rGB2 | gamlss_continuous |
| GG | d/p/q/rGG | gamlss_continuous |
| IG/IGAMMA | d/p/q/rIG, d/p/q/rIGAMMA | gamlss_continuous |
| GP | d/p/q/rGP | gamlss_continuous |
| JSUo/JSU | d/p/q/rJSUo, d/p/q/rJSU | gamlss_continuous |
| BCCG/BCT/BCPE | d/p/q/rBCCG, d/p/q/rBCT, d/p/q/rBCPE | gamlss_boxcox |
| PO/BI | d/p/q/rPO, d/p/q/rBI | gamlss_discrete |
| GEOM/GEOMo | d/p/q/rGEOM* | gamlss_discrete |
| NBI/NBII | d/p/q/rNBI, d/p/q/rNBII | gamlss_discrete |
| ZIP/ZIP2/ZAP | d/p/q/rZIP*, d/p/q/rZAP | gamlss_discrete |
| ZINBI/ZANBI | d/p/q/rZINBI, d/p/q/rZANBI | gamlss_discrete |
| ZIBI/ZABI | d/p/q/rZIBI, d/p/q/rZABI | gamlss_discrete |
| BB/BNB | d/p/q/rBB, d/p/q/rBNB | gamlss_discrete |
| PIG | d/p/q/rPIG | gamlss_discrete |

### v0.2.0 additions

| R family | Fortran routines | Module |
|---|---|---|
| GIG | dGIG, pGIG, qGIG, rGIG | gamlss_continuous_v02 |
| SHASHo | dSHASHo, pSHASHo, qSHASHo, rSHASHo | gamlss_continuous_v02 |
| SHASH | dSHASH, pSHASH, qSHASH, rSHASH | gamlss_continuous_v02 |
| SIMPLEX | dSIMPLEX, pSIMPLEX, qSIMPLEX, rSIMPLEX | gamlss_continuous_v02 |
| SEP/SEP1/SEP2 | d/p/q/rSEP* | gamlss_continuous_v02 |
| SEP3/SEP4 | d/p/q/rSEP3, d/p/q/rSEP4 | gamlss_continuous_v02 |
| ST1/ST2 | d/p/q/rST1, d/p/q/rST2 | gamlss_continuous_v02 |
| ST3/ST4 | d/p/q/rST3, d/p/q/rST4 | gamlss_continuous_v02 |
| ST5 | dST5, pST5, qST5, rST5 | gamlss_continuous_v02 |
| NET | dNET, pNET, qNET, rNET | gamlss_continuous_v02 |
| GPO | dGPO, pGPO, qGPO, rGPO | gamlss_discrete_v02 |
| DPO | dDPO, pDPO, qDPO, rDPO | gamlss_discrete_v02 |
| DEL | dDEL, pDEL, qDEL, rDEL | gamlss_discrete_v02 |
| SI | dSI, pSI, qSI, rSI | gamlss_discrete_v02 |
| SICHEL | dSICHEL, pSICHEL, qSICHEL, rSICHEL | gamlss_discrete_v02 |
| YULE | dYULE, pYULE, qYULE, rYULE | gamlss_discrete_v02 |
| WARING | dWARING, pWARING, qWARING, rWARING | gamlss_discrete_v02 |
| ZIPF | dZIPF, pZIPF, qZIPF, rZIPF | gamlss_discrete_v02 |

The v0.2 numerical helper module `gamlss_v02_numerics` exposes adaptive
integration, bracketed root finding, and real-order modified Bessel K support.

Optional `lower_tail`, `log_p`, and `log_density` arguments are provided where
applicable. Scalar routines are `elemental` when their implementation permits it.

Because `gamlss_base` already had a generic Zipf helper named `dzipf`, the
umbrella module exposes those old low-level helpers as `dzipf_base`,
`pzipf_base`, `qzipf_base`, and `rzipf_base`; `dZIPF` etc. refer to the GAMLSS
family parameterization.

## Links

`gamlss_links` supplies common GAMLSS links through `linkfun`, `linkinv`, and
`mu_eta`, with constants for identity, log, logit, probit, cauchit, cloglog,
square-root, inverse, `1/mu^2`, `mu^2`, shifted-log, `[-1,1]`, `(0,2]`, and
`(0,5]` transformations.

## Fitting

`fit_gamlss(y, x_mu, family, result, ...)` accepts independent design matrices
for `mu`, `sigma`, `nu`, and `tau`.

The 62 supported family constants are:

```text
GAMLSS_NO GAMLSS_GA GAMLSS_BE GAMLSS_NBI GAMLSS_NBII GAMLSS_ZIP
GAMLSS_GG GAMLSS_EGB2 GAMLSS_GB2 GAMLSS_JSUO GAMLSS_TF GAMLSS_PIG
GAMLSS_BEINF GAMLSS_WEI GAMLSS_LNO GAMLSS_BCCG GAMLSS_BCT GAMLSS_BCPE
GAMLSS_GIG GAMLSS_SHASHO GAMLSS_SHASH GAMLSS_SIMPLEX
GAMLSS_SEP GAMLSS_SEP1 GAMLSS_SEP2 GAMLSS_SEP3 GAMLSS_SEP4
GAMLSS_ST1 GAMLSS_ST2 GAMLSS_ST3 GAMLSS_ST4 GAMLSS_ST5 GAMLSS_NET
GAMLSS_GPO GAMLSS_DPO GAMLSS_DEL GAMLSS_SI GAMLSS_SICHEL
GAMLSS_YULE GAMLSS_WARING GAMLSS_ZIPF
GAMLSS_ST3C GAMLSS_SN1 GAMLSS_SN2 GAMLSS_SST GAMLSS_GT GAMLSS_EXGAUS
GAMLSS_PARETO GAMLSS_PARETO1 GAMLSS_PARETO2 GAMLSS_PARETO2O
GAMLSS_PIG2 GAMLSS_ZIPIG GAMLSS_ZAPIG GAMLSS_ZISICHEL GAMLSS_ZASICHEL
GAMLSS_ZIBNB GAMLSS_ZABNB GAMLSS_ZAZIPF GAMLSS_GAF GAMLSS_NBF GAMLSS_ZINBF
```

`gamlss_fit_result_t` stores the coefficients, fitted parameter vectors,
log-likelihood, AIC, convergence state, iteration count, and numerical-Hessian
covariance matrix.

## v0.3.0 additions

| R family | Fortran routines | Module |
|---|---|---|
| ST3C | dST3C, pST3C, qST3C, rST3C | gamlss_continuous_v03 |
| SN1/SN2 | d/p/q/rSN1, d/p/q/rSN2 | gamlss_continuous_v03 |
| SST | dSST, pSST, qSST, rSST | gamlss_continuous_v03 |
| GT | dGT, pGT, qGT, rGT | gamlss_continuous_v03 |
| exGAUS | dexGAUS, pexGAUS, qexGAUS, rexGAUS | gamlss_continuous_v03 |
| PARETO/PARETO1 | d/p/q/rPARETO* | gamlss_continuous_v03 |
| PARETO2/PARETO2o | d/p/q/rPARETO2* | gamlss_continuous_v03 |
| GAF | dGAF, pGAF, qGAF, rGAF | gamlss_flexible_v03 |
| DBI | dDBI, pDBI, qDBI, rDBI | gamlss_discrete_v03 |
| PIG2 | dPIG2, pPIG2, qPIG2, rPIG2 | gamlss_discrete_v03 |
| NBF | dNBF, pNBF, qNBF, rNBF | gamlss_flexible_v03 |
| ZINBF | dZINBF, pZINBF, qZINBF, rZINBF | gamlss_flexible_v03 |
| ZIPIG/ZAPIG | d/p/q/rZIPIG, d/p/q/rZAPIG | gamlss_discrete_v03 |
| ZISICHEL/ZASICHEL | d/p/q/rZISICHEL, d/p/q/rZASICHEL | gamlss_discrete_v03 |
| ZIBB/ZABB | d/p/q/rZIBB, d/p/q/rZABB | gamlss_discrete_v03 |
| ZIBNB/ZABNB | d/p/q/rZIBNB, d/p/q/rZABNB | gamlss_discrete_v03 |
| ZAZIPF | dZAZIPF, pZAZIPF, qZAZIPF, rZAZIPF | gamlss_discrete_v03 |

`fit_gamlss` now supports 62 family constants. In addition, `gamlss_fit_v03`
provides `fit_dbi`, `fit_zibb`, and `fit_zabb` for models with a supplied
binomial denominator vector.

The low-level Pareto helpers originally present in `gamlss_base` are exposed
through the umbrella module as `dpareto1_base`, `ppareto1_base`,
`qpareto1_base`, `dpareto2_base`, `ppareto2_base`, and `qpareto2_base` so the
GAMLSS `PARETO1`/`PARETO2` names remain unambiguous in case-insensitive Fortran.
