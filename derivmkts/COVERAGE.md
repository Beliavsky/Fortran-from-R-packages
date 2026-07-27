# Computational coverage

Source package: `derivmkts` 0.2.5.1.

## Black-Scholes and implied quantities

| Upstream export | Fortran interface | Status |
|---|---|---|
| `bscall`, `bsput` | `bscall`, `bsput` | Implemented |
| `assetcall`, `assetput` | Same names | Implemented |
| `cashcall`, `cashput` | Same names | Implemented |
| `bscallimpvol`, `bsputimpvol` | Same names | Implemented with expanding positive-volatility bracket |
| `bscallimps`, `bsputimps` | Same names | Implemented |
| `bsopt` | `bsopt`, `bs_call_greeks`, `bs_put_greeks` | Implemented |
| `greeks`, `greeks2` | `greeks`, `greeks2`, `numerical_greeks` | Implemented through a typed procedure callback |

## Bonds

| Upstream export | Fortran interface | Status |
|---|---|---|
| `bondpv` | `bondpv` | Implemented |
| `bondyield` | `bondyield` | Implemented |
| `duration` | `duration` | Implemented |
| `convexity` | `convexity` | Implemented |

## Asian options

| Upstream export | Fortran interface | Status |
|---|---|---|
| `geomavgprice`, `geomavgpricecall`, `geomavgpriceput` | Same names | Implemented |
| `geomavgstrike`, `geomavgstrikecall`, `geomavgstrikeput` | Same names | Implemented |
| `arithasianmc` | `arithasianmc` | Implemented |
| `geomasianmc` | `geomasianmc` | Implemented |
| `arithavgpricecv` | `arithavgpricecv` | Implemented |

## Barrier and perpetual claims

All standard and binary barrier exports are implemented under the same lower
case names:

- `calldownin`, `calldownout`, `putdownin`, `putdownout`
- `callupin`, `callupout`, `putupin`, `putupout`
- `dicall`, `docall`, `diput`, `doput`
- `uicall`, `uocall`, `uiput`, `uoput`
- `cashdicall`, `cashdocall`, `cashdiput`, `cashdoput`
- `assetdicall`, `assetdocall`, `assetdiput`, `assetdoput`
- `cashuicall`, `cashuocall`, `cashuiput`, `cashuoput`
- `assetuicall`, `assetuocall`, `assetuiput`, `assetuoput`
- `dr`, `ur`, `drdeferred`, `urdeferred`
- `callperpetual`, `putperpetual`

Perpetual functions return `perpetual_result`, containing both price and optimal
exercise barrier.

## Compound options

| Upstream export | Fortran interface | Status |
|---|---|---|
| `binormsdist` | `binormsdist` | Implemented by adaptive one-dimensional Gaussian integration |
| `calloncall`, `putoncall`, `callonput`, `putonput` | Same names | Implemented; return `compound_result` |
| `optionsoncall`, `optionsonput` | Same names | Implemented as output-argument subroutines |

## Jump diffusion

| Upstream export | Fortran interface | Status |
|---|---|---|
| `mertonjump` | `mertonjump` | Implemented |
| `assetjump` | `assetjump` | Implemented |
| `cashjump` | `cashjump` | Implemented |

The Poisson mixture is accumulated recursively until its omitted mass is below
`1e-12`.

## Lattices, simulation, and demonstrations

| Upstream export | Fortran interface | Status |
|---|---|---|
| `binomopt` | `binomopt` | Implemented, including trees, hedge positions, probabilities, and Greeks |
| `binomplot` | Tree arrays returned by `binomopt` | Graphics excluded |
| `simprice` | Generic `simprice` | Implemented for scalar volatility and covariance-matrix inputs, with optional jumps |
| `quincunx` | `quincunx` | Counts and theoretical frequencies implemented; animation excluded |

## R-only infrastructure excluded

- Base-graphics and PDF output
- Animated plotting and interactive delays
- R formula/language introspection
- Vector recycling and data-frame reshaping
- R global random-state save/restore

These exclusions do not remove a numerical pricing, simulation, distribution,
or hedge calculation.
