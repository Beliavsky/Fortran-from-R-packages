# Notices

This work translates the computational code of RcppNumerical 0.7-0 into
modern Fortran.

- RcppNumerical package and wrappers: Copyright 2016-2026 Yixuan Qiu,
  GPL-2.0-or-later for the package and MPL-2.0 for the public callback wrappers.
- One-dimensional integration tables and algorithm: modified from
  NumericalIntegration by Sreekumar T. Balan, Matt Beall, and Mark Sauder,
  MPL-2.0.
- Infinite-interval integration contribution: Ralf Stubner.
- Multidimensional integration: Cuhre from the Cuba library by Thomas Hahn,
  LGPL-3.0. The Fortran module ports the default rule-13, rule-11, and rule-9
  paths used by the public RcppNumerical wrapper.
- L-BFGS interface concepts: RcppNumerical/LBFGSpp by Yixuan Qiu, MIT.
- Bundled L-BFGS implementation: see `dependencies/lbfgs/NOTICE.md` and its
  license directory.
- Bundled L-BFGS-B implementation: see `dependencies/lbfgsb3/NOTICE.md` and
  its license directory.

The complete supplied R package source is retained under `upstream/` for
attribution and source comparison.
