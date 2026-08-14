# Licensing

The original `mlrMBO` package is distributed under BSD-2-Clause with the
copyright-holder file included as `LICENSE-UPSTREAM`. Files in this port that
are derived from mlrMBO retain that notice; see `LICENSE-MLRMBO-BSD-2-CLAUSE`.

This source distribution also contains the separately translated
DiceKriging numerical package under `vendor/DiceKriging-fortran-v0.1.0`.
DiceKriging permits GPL-2 or GPL-3. This combined distribution elects the
GPL-3 option, so the distribution as a whole is released under GPL-3. The
complete GPL-3 text is in `LICENSE`; DiceKriging's original GPL-2/GPL-3 files
remain in its vendor directory.

The original mlrMBO C hypervolume source (`src/hv.c`) is GPL-2-or-later. It
is not copied into this Fortran port: `mlrmbo_multiobjective` contains an
independently written recursive exact hypervolume implementation. The
original upstream package remains available in the supplied source archive
and is identified in `UPSTREAM.md`.
