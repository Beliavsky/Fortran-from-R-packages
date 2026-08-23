# Licensing

## Active library

`sgt-fortran` is a translation of the `sgt` R package 2.0, whose DESCRIPTION
states:

```text
License: GPL (>= 3)
```

The translated source is therefore distributed under GPL-3.0-or-later. The
full GPL version 3 text is in `LICENSE`.

## Upstream source snapshot

`upstream/sgt-master.zip` is the supplied original package snapshot and retains
its original notices and GPL (>= 3) licensing.

## Reference dependency snapshots

The following supplied translations are retained for provenance/reference but
are not compiled by this FPM project:

- `reference-dependencies/optimx-fortran.zip`
- `reference-dependencies/numDeriv-fortran.zip`

The `optimx` translation is marked GPL-2.0-only. The `numDeriv` translation is
marked GPL-2.0-or-later, although its bundled original package metadata says
`GPL-2`. To avoid changing or guessing those licenses, neither snapshot is
linked into the active GPL-3.0-or-later library.
