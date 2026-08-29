# Notices and attribution

This project is a Fortran translation of computational code from the R package **gmm** version 1.9-1.

Upstream package:
- Author and maintainer: Pierre Chausse
- License: GPL (>= 2), represented here as GPL-2.0-or-later
- Upstream source metadata is retained under `upstream/`.

Please cite the upstream package/methodology as requested by its `CITATION` file, in particular:

Pierre Chausse (2010), *Computing Generalized Method of Moments and Generalized Empirical Likelihood with R*, Journal of Statistical Software 34(11), DOI 10.18637/jss.v034.i11.

The statistical procedures also build on the literature cited by the upstream package, including Hansen's GMM, continuously updated GMM, and generalized empirical-likelihood work by Smith, Kitamura, Newey and Smith, and Anatolyev. The upstream `DESCRIPTION`, `CITATION`, R computational sources, and native sources are retained for provenance.

`src/r_mod.F90` is a formatting-only build copy of the user-supplied MIT-licensed `r_mod.f90`; the exact supplied source is retained as `upstream/r_mod-original.f90`. No copyright-holder line was present in the supplied helper, so none is invented. See `LICENSES/MIT-r_mod.txt`.
