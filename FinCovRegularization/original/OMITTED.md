# R code not translated

`R/Generic.CovCv.R` contains `plot.CovCv`, `print.CovCv`, and `summary.CovCv`.
These are plotting/display and R S3 infrastructure rather than computational
algorithms, so they are outside the Fortran port's scope.

The binary R dataset in `data/m.excess.c10sp9003.rda` is preserved but is not
loaded by the Fortran library.
