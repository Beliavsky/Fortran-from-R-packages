# Notice and provenance

This directory is a modern Fortran translation of the computational portions of
**Rssa 1.1**, supplied as `Rssa-master.zip`.

## Upstream attribution

The R package names Anton Korobeynikov, Alex Shlemov, Konstantin Usevich and
Nina Golyandina as authors. The retained source contains, among others, these
copyright notices:

- Copyright (c) 2008-2016 Anton Korobeynikov.
- Copyright (c) 2009, 2013 Konstantin Usevich.
- Copyright (c) 2012-2018 Alex (Alexander) Shlemov.
- Copyright (c) 2015 Nina Golyandina.

The exact notices, years, author addresses and contributor history remain in the
unmodified files under `upstream/` and are authoritative.

Rssa declares `GPL (>= 2)`, and the source headers permit redistribution under
GNU GPL version 2 or, at the recipient's option, any later version. The
translation is therefore distributed as `GPL-2.0-or-later`. The GPL version 2
text is included as `LICENSE`.

## Supplied source

SHA-256 of the supplied archive:

```text
e68590dad5d5634cc14a613dcde78ee4f60df903b5fcf6d193661ab5b011f4d7  Rssa-master.zip
```

`provenance/upstream-files.sha256` records hashes of the retained upstream
files. The original sources are kept under `upstream/`; they are provenance,
not FPM compilation inputs.

## Reused sibling packages

The translation does not vendor numerical dependencies. In the intended
`Fortran-from-R-packages` repository layout it uses:

- `../svd` for SVD/eigensolver compatibility operations. That package in turn
  reuses the repository's iterative eigensolver and pure-Fortran LAPACK stack.
- `../rfortran-linalg` for QR, least-squares, linear solves and general
  eigenvalue calculations.

No BLAS, LAPACK, ARPACK, FFTW, `r.f90`, `r_mod.f90`, or translated dependency
source is copied into this package.

## Translation approach

Rssa's FFTW/external-pointer acceleration layer is intentionally not copied.
Trajectory, Hankel, block-Hankel and Toeplitz operations are expressed as dense
portable Fortran reference algorithms. This changes performance and iteration
paths but preserves the mapped numerical definitions. See `PORTING_NOTES.md`
and `API_COVERAGE.md` for material compatibility differences.
