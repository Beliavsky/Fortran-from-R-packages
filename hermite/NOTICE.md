# Notice

This project translates the computational code of the R package `hermite`
1.2.1 by David Morina Soler, Manuel Higueras, Pedro Puig, and Maria Oliveira
into modern Fortran/FPM.

Distribution recurrences and Edgeworth/Cornish-Fisher approximations are
preserved. R formula parsing, model frames, S3 summary/printing, and the
`maxLik` dependency are replaced by a standalone numerical regression API.
