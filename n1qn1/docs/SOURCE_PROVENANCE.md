# Source provenance

## Supplied package archive

- Archive: `n1qn1c-main.zip`
- SHA-256: `e5a718be5243e7c1737d6209e97a8a41aa06a17fa4b924a233bd031be4ebc4bc`
- Package version in `DESCRIPTION`: `6.0.1-14`
- Package license: CeCILL version 2.1

The extracted package sources used for translation are retained in
`original/n1qn1c/`.

## Scilab reference sources

The translation was also compared against the Scilab 2026.1.0 source files:

- `scilab/modules/optimization/src/fortran/n1qn1.f`
- `scilab/modules/optimization/src/fortran/n1qn1a.f`
- `scilab/modules/optimization/src/fortran/majour.f`

The first two downloaded source files are retained in `original/scilab/`. The
supplied package's `src/n1qn1_all.c` contains the complete f2c-derived `majour`
routine used for the translation and is retained in `original/n1qn1c/src/`.

The Scilab headers state that these files are available under GNU GPL version 2
pursuant to article 5.3.4 of CeCILL version 2.1 and remain available under
CeCILL version 2.1.
