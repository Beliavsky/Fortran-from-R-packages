# API map

All 194 names exported by the upstream `NAMESPACE` are implemented.

Use `use extra_distr` to access the complete public API.

## Continuous univariate distributions

Fortran module: `extra_distr_continuous`

```text
dbetapr dbhatt dfatigue dfrechet dgev dgompertz dgpd dgumbel
dhcauchy dhnorm dht dhuber dinvchisq dinvgamma dkumar dlaplace
dlomax dlst dnsbeta dpareto dpower dprop drayleigh dsgomp
dslash dtnorm dtriang dwald pbetapr pbhatt pfatigue pfrechet
pgev pgompertz pgpd pgumbel phcauchy phnorm pht phuber
pinvchisq pinvgamma pkumar plaplace plomax plst pnsbeta ppareto
ppower pprop prayleigh psgomp pslash ptnorm ptriang pwald
qbetapr qfatigue qfrechet qgev qgompertz qgpd qgumbel qhcauchy
qhnorm qht qhuber qinvchisq qinvgamma qkumar qlaplace qlomax
qlst qnsbeta qpareto qpower qprop qrayleigh qtlambda qtnorm
qtriang rbetapr rbhatt rfatigue rfrechet rgev rgompertz rgpd
rgumbel rhcauchy rhnorm rht rhuber rinvchisq rinvgamma rkumar
rlaplace rlomax rlst rnsbeta rpareto rpower rprop rrayleigh
rsgomp rslash rtlambda rtnorm rtriang rwald
```

## Discrete univariate distributions

Fortran module: `extra_distr_discrete`

```text
dbbinom dbern dbnbinom dcat ddgamma ddlaplace ddnorm ddunif
ddweibull dgpois dlgser dnhyper dskellam dtbinom dtpois dzib
dzinb dzip pbbinom pbern pbnbinom pcat pdgamma pdlaplace
pdnorm pdunif pdweibull pgpois plgser pnhyper ptbinom ptpois
pzib pzinb pzip qbern qcat qdunif qdweibull qlgser
qnhyper qtbinom qtpois qzib qzinb qzip rbbinom rbern
rbnbinom rcat rcatlp rdgamma rdlaplace rdnorm rdunif rdweibull
rgpois rlgser rnhyper rsign rskellam rtbinom rtpois rzib
rzinb rzip
```

## Mixtures and multivariate distributions

Fortran module: `extra_distr_multivariate`

```text
dbvnorm dbvpois ddirichlet ddirmnom dmixnorm dmixpois dmnom dmvhyper
pmixnorm pmixpois rbvnorm rbvpois rdirichlet rdirmnom rmixnorm rmixpois
rmnom rmvhyper
```

## Shared utilities

- `seed_rng`: deterministic initialization of the intrinsic RNG.
- `dp`: binary64 working kind (`kind(1.0d0)`).
- `extra_distr_math`: internal probability and special-function kernel.
- `extra_distr_rng`: internal scalar random-generation kernel.
