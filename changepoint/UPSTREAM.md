# Upstream provenance

## Source package

- Package: `changepoint`
- Version: `2.3`
- Date: `2024-11-02`
- CRAN publication: `2024-11-04`
- Homepage: <https://github.com/rkillick/changepoint/>
- CRAN: <https://CRAN.R-project.org/package=changepoint>
- Input archive used for this translation: `changepoint-master(20260829-210932).zip`
- Input archive SHA-256: `130044d3ab2c194c9a2429fb0303771e87c6406c59dc66201a516db7c67201ae`

The original `DESCRIPTION`, `NAMESPACE`, `NEWS`, and `inst/CITATION` are retained
under `upstream/` without being compiled into the Fortran library.

## Computational source families used

The translation was derived from the package's computational R sources,
including `PELT_one_func_minseglen.R`, `BinSeg_one_func_minseglen.R`,
`SegNeigh_one_func_minseglen.R`, `single.norm.R`, `multiple.norm.R`,
`single.nonparametric.R`, `multiple.nonparametric.R`, `exp.R`, `gamma.R`,
`poisson.R`, `CROPS.R`, `range_of_penalties.R`, `CptReg.R`, `fit.R`,
`penalty_decision.R`, and `decision.R`, together with the native C kernels in
`src/PELT_one_func_minseglen.c`, `src/BinSeg_one_func_minseglen.c`,
`src/cost_general_functions.c`, and `src/C_cptreg.c`.

`upstream/SOURCE_SHA256S.txt` records hashes of the relevant upstream R and C
source files so later audits can tie this translation to the exact input tree.

## Indexing convention

The R package documents a changepoint as the **last observation of the preceding
segment**. The Fortran API preserves that 1-based convention. For example,
changepoints `[4, 8]` in a 12-observation vector define segments `1:4`, `5:8`,
and `9:12`.
