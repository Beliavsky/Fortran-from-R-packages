# Translation coverage

The upstream namespace contains 116 exports. Forty-two are plotting, printing,
summary, density-display, or other presentation-oriented entries. The Fortran
library focuses on the numerical surface and exposes typed equivalents for the
major computational exports.

## Implemented directly or through a typed equivalent

### Parametric and regression mixtures

`normalmixEM`, `normalmixEM2comp`, `normalmixMMlc`,
`tauequivnormalmixEM`, `mvnormalmixEM`, `gammamixEM`, `multmixEM`,
`repnormmixEM`, `regmixEM`, `regmixEM.lambda`, `regmixEM.loc`,
`regmixEM.mixed`, `logisregmixEM`, `poisregmixEM`, `segregmixEM`,
`hmeEM`, `flaremixEM`, and `try.flare`.

### Semiparametric and reliability mixtures

`npEM`, `npEMindrep`, `npEMindrepbw`, `npMSL`, `spEM`, `mvnpEM`,
`spEMsymloc`, `spEMsymlocN01`, `spregmix`, `expRMM_EM`,
`weibullRMM_SEM`, and `spRMM_SEM`.

### Initialization interfaces

`normalmix.init`, `mvnormalmix.init`, `gammamix.init`, `multmix.init`,
`repnormmix.init`, `regmix.init`, `regmix.lambda.init`,
`regmix.mixed.init`, `logisregmix.init`, `poisregmix.init`,
`segregmix.init`, and `flaremix.init` are represented by one-step typed
initialization procedures.

### Selection, bootstrap, Bayesian, and tests

`multmixmodel.sel`, `regmixmodel.sel`, `repnormmixmodel.sel`,
`boot.comp`, `boot.se`, `regmixMH`, `post.beta`, `regcr`,
`test.equality`, and `test.equality.mixed`.

The two bootstrap procedures are currently normal-mixture implementations.

### Utilities and simulation

`aug.x`, `compCDF`, `ddirichlet`, `depth`, `dexpmixt`, `dmvnorm`,
`logdmvnorm`, `ellipse`, `ise.npEM`, `lambda`, `lambda.pert`, `ldc`,
`ldmult`, `makemultdata`, `matsqrt`, `parse.constraints`, `perm`,
`rexpmix`, `rlnormscalemix`, `rmvnorm`, `rmvnormmix`, `rnormmix`,
`rweibullmix`, `normmix.sim`, `normmixrm.sim`, `wIQR`, `wkde`, and
`wquantile`.

`density.npEM`, `density.spEM`, and the numerical part of FDR calculation are
available as grid interpolation and posterior-FDR procedures.

## Omitted from the compiled API

- base-R and Plotly plotting functions
- print and summary methods
- formula/S3 dispatch and object conversion
- `mixturegram`, whose principal purpose is visual model comparison
- bundled R datasets and `.RData` serialization
- interactive display helpers

## Important adapted routines

The following are computationally represented but are not exact line-by-line
ports: `flaremixEM`, `npMSL`, `spRMM_SEM`, `regmixEM.mixed`,
`segregmixEM`, `regmixMH`, and generic multi-family bootstrap dispatch.
Details are in `PORTING.md`.
