# NOTICE

This project is a computational source translation of the R package `distr` 2.9.7.

Upstream `distr` is credited to the authors and copyright holders named in `upstream/DESCRIPTION` and should be cited as described in `upstream/CITATION`. The upstream package declares LGPL-3.

The Kolmogorov-Smirnov routines in `src/distr_ks.f90` are translated from `distr/src/ks.c`, which was taken from R Core and is copyright the R Core Team under GPL-2.0-or-later. The original C source is retained verbatim in `upstream/ks.c`. `src/distr_qq.f90` translates the non-graphical QQ confidence-band calculations and is GPL-3.0-or-later because it depends on that KS module.

See `LICENSES.md` and the `COPYING*` files for details.
