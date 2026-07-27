# Origin and licensing

This project is a modern Fortran translation of the computational routines in:

- R package: `fMultivar`
- Version: 4031.84
- Package date: 2023-07-07
- Source archive used: `fMultivar-master.zip`
- Original package URL: `https://www.rmetrics.org`

The original `DESCRIPTION` file declares:

```text
License: GPL (>= 2)
```

The historical R source headers use older "GNU Library General Public License"
wording. The package-level CRAN license declaration is GPL version 2 or later,
which is preserved here as `GPL-2.0-or-later`.

Original authors listed by the package:

- Diethelm Wuertz
- Tobias Setz
- Stefan Theussl
- Yohan Chalabi
- Martin Maechler
- CRAN team contributors

The translation is independent source code written in modern Fortran from the
published R formulas and documented behavior. It does not copy code from the
external `mvtnorm`, `sn`, or `cubature` packages. Their imported algorithms are
represented by separately implemented numerical analogues described in
`README.md` and `API_MAP.md`.

`ORIGINAL_DESCRIPTION` is retained for provenance. `LICENSE` contains the
complete GNU GPL version 2 text.
