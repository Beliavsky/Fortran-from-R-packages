# Licenses

## orthopolynom-derived translation

The upstream R package `orthopolynom` 1.0-6.1 declares:

```text
License: GPL (>= 2)
```

Accordingly, the Fortran files in `src/`, tests, examples, and documentation
that derive from the upstream computational code are distributed under
GPL-2.0-or-later. A copy of GPL version 2 is in `licenses/GPL-2.0.txt`.

The original package metadata and source are retained in
`upstream/orthopolynom-master/` and `upstream/orthopolynom-master.zip`.

## polynom-fortran dependency

`vendor/polynom-fortran` is the user-supplied translation of R package
`polynom`. That translation declares GPL-2.0-only. Its own license,
copyright, notices, upstream source, and provenance files are retained inside
the vendored directory.

Because this release includes and links against that GPL-2.0-only dependency,
a combined executable is subject to the compatible GPL version 2 terms.
