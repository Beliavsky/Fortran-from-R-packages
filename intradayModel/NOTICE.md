# Notice

This project is a modern Fortran translation of `intradayModel` 0.0.1 by Shengjie Xiu,
Yifan Yu, and Daniel P. Palomar.

The original package declares the Apache License, Version 2.0. The original source,
metadata, documentation, tests, and data files are retained in `original/` except for
the generated R build artifact directory.

The Fortran translation preserves the numerical model while replacing R lists, S3
objects, and time-series containers with typed array-based interfaces. Plotting code is
not included in the compiled library.
