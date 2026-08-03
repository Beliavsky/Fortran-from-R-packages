# Translation coverage

The upstream namespace exports 44 names.

Forty-three are represented in the compiled Fortran API:

- Solvers/orchestration: `optimr`, `optimx`, `opm`, `multistart`, `polyopt`,
  `proptimr`, `Rvmmin`, `Rvmminu`, `Rvmminb`, `Rcgmin`, `Rcgminu`, `Rcgminb`,
  `nvm`, `ncg`, `hjn`, `tn`, `tnbc`, `snewton`, `snewtm`
- Bounds/search: `bmchk`, `bmstep`, `axsearch`
- Derivatives: `grfwd`, `grback`, `grcentral`, `grnd`, `grpracma`, `gHgen`,
  `gHgenb`
- Checks: `fnchk`, `grchk`, `hesschk`, `kktchk`, `optchk`, `pd_check`,
  `scalechk`
- Utilities: `ctrldefault`, `dispdefault`, `checksolver`, `checkallsolvers`,
  `opm2optimr`, `optimr2opm`, and `optsp`

`coef<-` is an R S3 replacement method and is omitted.

The internal R setup, S3 summary/subsetting/data-frame methods, package
installation checks, and dynamic access to optional R namespaces are retained
under `original/` but are not compiled.
