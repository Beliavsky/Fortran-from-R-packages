# Notice

This project is a modern Fortran translation of the computational portions of
`mfGARCH` version 0.2.2 by Onno Kleen. The upstream package is distributed
under the MIT License. The complete upstream source snapshot used for the
translation is retained in `original/mfGARCH-master/`.

The translation omits plotting, R S3 classes, R data-frame/date handling,
bundled data sets as callable package objects, and R-specific printing. The
underlying source files and data remain available in the retained upstream
snapshot for provenance.

The Fortran optimizer and numerical linear algebra were implemented directly
for this translation and do not copy the GPL-licensed `maxLik` or `numDeriv`
implementations used by the R package. The resulting Fortran project therefore
remains MIT-licensed.
