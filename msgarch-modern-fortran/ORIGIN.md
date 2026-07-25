# Origin and Licensing

This project is a modern Fortran computational translation of the R package `MSGARCH`.

The attached source package reported:

```text
Package: MSGARCH
Version: 2.51
License: GPL (>= 2)
```

The original metadata is retained in:

- `reference/DESCRIPTION.original`
- `reference/NAMESPACE.original`
- `reference/COPYRIGHTS.original`
- `reference/README.original.md`

The translation therefore uses the SPDX identifier:

```text
GPL-2.0-or-later
```

`LICENSE` contains the GNU General Public License version 2 text. Every Fortran source, application, example, and test file includes the SPDX identifier and a GPL version 2-or-later notice. `test/check_license.sh` enforces these headers.

The Fortran project does not incorporate the R class infrastructure, plotting methods, or bundled R object representation. It translates and tests numerical algorithms and exposes their results through plain Fortran arrays and derived types.
