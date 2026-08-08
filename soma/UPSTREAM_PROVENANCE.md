# Upstream provenance

This project is a Fortran translation of the computational code in:

- Package: `soma`
- Version: 1.2.0
- Date: 2022-05-01
- Author/Maintainer: Jon Clayden
- Upstream URL: `https://github.com/jonclayden/soma/`
- Upstream license declaration: GPL-2

The complete uploaded upstream tree is preserved verbatim under:

```text
original/soma-master/
```

The translated optimizer is based primarily on `R/soma.R`. The original
package has `NeedsCompilation: no`; there is no native C/C++/Fortran numerical
backend to preserve separately.

Algorithm references retained from the upstream documentation:

- I. Zelinka (2004), SOMA - self-organizing migrating algorithm.
- Q.B. Diep (2019), SOMA Team To Team Adaptive (T3A).
- Q.B. Diep, I. Zelinka and S. Das (2019), Pareto-Based SOMA.
