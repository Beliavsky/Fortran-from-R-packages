# Upstream provenance

Translated package: `FNN` 1.1.4.1, CRAN package metadata dated 2023-12-31 and
packaged 2024-09-22.

Upstream authors listed in `DESCRIPTION`:

- Alina Beygelzimer, Sham Kakadet, and John Langford (cover-tree library);
- Sunil Arya and David Mount (ANN library 1.1.2);
- Shengqiao Li.

The original package declares `GPL (>= 2)`.

The embedded ANN code carries this notice:

> ANN Copyright (c) 1997-2010 University of Maryland and Sunil Arya and David
> Mount. All Rights Reserved.

and is distributed under LGPL 2.1 or later.  The upstream `COPYRIGHTS` file and
full original sources are retained under `upstream/FNN-1.1.4.1/`.

The supplied `mvtnorm-fortran` project was inspected only to determine whether
it was a runtime dependency.  FNN uses `mvtnorm` solely in a documentation
example, so it is not incorporated into this distribution.
