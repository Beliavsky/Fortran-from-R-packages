# Origin and licensing

This project is a computational translation of:

- Package: `gogarch`
- Version: 0.7-6
- Date: 2026-01-24
- Title: Generalized Orthogonal GARCH (GO-GARCH) Models
- Author and maintainer: Bernhard Pfaff
- Original license declaration: `GPL (>= 2)`

The original `DESCRIPTION` and `NAMESPACE` files are retained in `reference/`.

The translation is distributed under `GPL-2.0-or-later`. Every Fortran source
file contains:

```text
SPDX-License-Identifier: GPL-2.0-or-later
```

and an explicit GNU GPL version 2-or-later notice with original-author
attribution.

The numerical implementation was written for this translation. It does not
copy R S4 infrastructure, plotting code, bundled data objects, or the `fastICA`
and `fGarch` dependency implementations.
