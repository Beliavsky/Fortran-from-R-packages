# Notice and provenance

This directory contains an independent modern-Fortran translation of the portable
table semantics of **tibble 3.3.1**, *Simple Data Frames*.

Upstream project: https://github.com/tidyverse/tibble

Upstream authors and contributors include Kirill Müller, Hadley Wickham, Romain
Francois, Jennifer Bryan, and Posit Software, PBC. Relevant upstream `DESCRIPTION`
metadata, the exported API and S3-method portion of `NAMESPACE`, and the CRAN MIT
license metadata are preserved in `upstream/`.

The homogeneous real/integer tibble API in the MIT-licensed
[`r_mod.f90`](https://github.com/Beliavsky/R-to-Fortran/blob/main/src/r_mod.f90)
was reviewed and used as a design and testing
reference. Heterogeneous tagged-vector storage is supplied by the sibling MIT-licensed
`vctrs` translation; this package does not vendor either that source or the full
R-to-Fortran runtime.

The translated source is distributed under the MIT license, matching upstream tibble
and the reusable runtime material. R-specific S3, tidy-evaluation, vctrs, pillar,
interactive, and arbitrary-object behavior is intentionally excluded.
