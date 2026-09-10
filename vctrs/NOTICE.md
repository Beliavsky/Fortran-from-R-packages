# Notice and provenance

This directory contains an independent modern-Fortran translation of portable vector
semantics from **vctrs 0.7.3**, *Vector Helpers*.

Upstream project: https://github.com/r-lib/vctrs

Upstream authors are Hadley Wickham, Lionel Henry, and Davis Vaughan. The data.table
team is credited upstream for radix-sort work, and Posit Software, PBC is an upstream
copyright holder and funder. Relevant upstream package and license metadata are
preserved in `upstream/`.

The Fortran implementation follows the documented vctrs invariants for type and size
stability but uses a new explicit tagged-vector representation. It does not translate
or embed upstream R runtime, C registration, S3 dispatch, ALTREP, or radix-sort code.
Consequently no data.table radix-sort source is present in this directory.

The translated source is distributed under the MIT license, matching upstream vctrs.
