# Notices and provenance

This project is a Fortran translation of the computational parts of the R package
`forecast` version 9.0.2 by Rob J Hyndman and contributors. The complete source snapshot
used for translation is retained in `upstream/forecast-master/`.

Upstream `forecast` declares the GNU General Public License version 3. A copy of GPL-3 is
included as `LICENSE`.

The package also bundles user-supplied Fortran translations under `dependencies/`:

- `fracdiff`, whose bundled port declares GPL-2.0-or-later;
- `urca`, whose bundled port declares GPL-2.0-or-later;
- `nnet`, whose bundled port declares GPL-2.0-only OR GPL-3.0-only.

For this GPL-3 combined distribution, the compatible GPL-3 licensing option is used where
applicable. Original authorship, notices, upstream sources and dependency license material
are retained inside their respective trees.

No plotting code from `forecast` has been translated. Retaining upstream plotting sources
inside `upstream/` is for license/provenance/audit purposes only; they are not built by FPM.
