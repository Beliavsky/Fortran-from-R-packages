# License and provenance

## Rmalschains / librealea translation

The attached `Rmalschains` package declares `License: GPL-3`.  Its documentation
states that Rmalschains is adapted from Daniel Molina Cabrera's librealea and
that both are licensed under GNU GPL version 3.

The translated Fortran source in `src/` is therefore distributed under
**GPL-3.0-only**.  The GPL v3 text is included as `COPYING`.

## Original bundled components

The original package's own detailed provenance file is retained verbatim as:

`ORIGINAL_LICENSES.txt`

and in:

`original/Rmalschains-master/inst/original_licenses_of_code_parts`

It identifies, among other material:

* `RmalschainsEvaluate.h`, adapted from RcppDE 0.1.0 (GPL >= 2);
* `origcmaes.h` / `origcmaes.cc`, Nikolaus Hansen's CMA-ES implementation
  (GPL-3 according to the package notice);
* portions of Robert Davies' newmat library, under its permissive no-restriction
  notice with disclaimer;
* Richard J. Wagner's `ConfigFile.h`, under an MIT-style permission notice.

The complete attached source tree is retained under `original/` so those
notices remain with the corresponding files.
