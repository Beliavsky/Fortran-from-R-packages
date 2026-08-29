# Notice and provenance

This project is a translation of the computational core of the R package `lmtest` version 0.9-40.

Upstream authors listed in `DESCRIPTION` include Torsten Hothorn, Achim Zeileis, Richard W. Farebrother, Clint Cummins, Giovanni Millo, and David Mitchell. The upstream license field is:

```text
GPL-2 | GPL-3
```

This translation preserves that license choice. The original `DESCRIPTION`, `NAMESPACE`, `NEWS`, `THANKS`, R computational sources, and `src/pan.f` are retained under `upstream/` for attribution and auditability.

`pan.f` is based on the amended version of Applied Statistics Algorithm AS 153 (AS R52), Farebrother (1984), with revisions attributed in the upstream source to Clint Cummins. Its algorithm has been translated to modern free-format Fortran in `src/lmtest_pan.f90` without changing the licensing of the inherited code.
