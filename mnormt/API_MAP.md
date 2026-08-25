# API map

| R `mnormt` | Fortran | Notes |
|---|---|---|
| `dmnorm` | `dmnorm`, `dmnorm_many` | vector / matrix observations |
| `pmnorm` | `pmnorm` | specialized d=2/3, adaptive Genz otherwise |
| `rmnorm` | `rmnorm` | native normal RNG + Cholesky |
| `sadmvn` | `sadmvn_prob` | returns value/error/status |
| `dmt` | `dmt`, `dmt_many` | real-valued df in density |
| `pmt` | `pmt` | real df for d=1; upstream integer rounding for d>1 |
| `rmt` | `rmt` | normal / chi-square mixture |
| `sadmvt` | `sadmvt_prob` | returns value/error/status |
| `biv.nt.prob` | `biv_nt_prob` | Genz bivariate kernel |
| `ptriv.nt` | `ptriv_nt` | Genz TVPACK kernel |
| `dmtruncnorm` | `dmtruncnorm` | normalized truncated density |
| `pmtruncnorm` | `pmtruncnorm` | truncated rectangle CDF |
| `rmtruncnorm` | `rmtruncnorm` | conditional-normal Gibbs sampler |
| `mom.mtruncnorm` | `mom_mtruncnorm` | raw moments + cumulants result type |
| `mom2cum` | `mom2cum` | cumulants through order four |
| `recintab` | `recintab` | flat column-major moment table + shape |
| `dmtrunct` | `dmtrunct` | truncated multivariate t density |
| `pmtrunct` | `pmtrunct` | truncated multivariate t CDF |
| `pd.solve` | `pd_solve` | inverse and optional log determinant |
| `sample_Mardia_measures` | `sample_mardia_measures` | numerical result type |
| `plot_fxy` | omitted | plotting only |
| `.onLoad` | omitted | R runtime only |
