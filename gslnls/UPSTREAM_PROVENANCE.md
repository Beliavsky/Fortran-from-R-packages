# Upstream provenance

Translated from the supplied source archive of:

- Package: `gslnls`
- Version: 1.4.2
- Date: 2025-09-27
- Author/Maintainer: Joris Chau
- Upstream package license declaration: `LGPL-3`
- Upstream system requirement: GNU Scientific Library (GSL) >= 2.3

The supplied source archive is preserved verbatim at:

```text
original/gslnls-master/
```

Important upstream computational sources include:

```text
src/nls.c
src/nls_large.c
src/trust.c
src/nls_irls.c
src/nls_mstart.c
src/fdjac.c
src/fdf.c
src/fdfvv.c
src/nls_fit.c
src/nls_utils.c
src/test_nls.f90
```

The modern Fortran source is a standalone translation/reimplementation of the
package-owned numerical behavior and the GSL algorithms it invokes. The original
source remains available in the archive for detailed comparison and attribution.

No upstream license file was present in the supplied package archive. The
translation includes the GNU Lesser General Public License version 3 text in
`LICENSE`, consistent with the package's `License: LGPL-3` declaration.
