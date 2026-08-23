# Licensing

## Active Fortran port

The translated Fortran source is distributed under **GPL-3.0-or-later**, matching the attached R package's declaration:

```text
License: GPL (>= 3)
```

A copy of GNU GPL version 3 is included as `LICENSE`. The "or later" permission comes from the upstream package declaration.

## Upstream Rnanoflann

Rnanoflann 0.0.3 is copyright Manos Papadakis and is declared GPL (>= 3).

## Bundled nanoflann in the upstream snapshot

The attached upstream archive contains `inst/include/nanoflann.hpp`, whose copyright/license notice names Marius Muja, David G. Lowe, and Jose Luis Blanco and permits redistribution under a BSD-style two-clause license. The exact upstream notice is preserved in `upstream/Rnanoflann-master.zip` under `inst/COPYRIGHT` and in the header itself.

The active Fortran implementation does not compile or link the bundled C++ nanoflann header.
