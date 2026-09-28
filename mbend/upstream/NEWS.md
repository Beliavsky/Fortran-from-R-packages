# mbend: Matrix bending

## Version: 1.2.1

* The 1st official version on CRAN

## Version: 1.2.2

* Modified for returning the matrix if it is already positive-definite.
* Skip other checks if the matrix is already positive-definite.

## Version: 1.2.3

* Set VignetteIndexEntry

## Version: 1.2.4

* Giving reports and statistics
* Iteration 1: if `small.positive` is greater than the smallest positive eigenvalue (x), replace it with x/10 and report the user.

## Version: 1.2.5

* Correcting citation Schaeffer (2010) to Schaeffer (2014)
* Reporting the average and the range of deviations
* Statistical reports for correlation matrices do not include diagonal elements.
* Updated README.md for changes in the previous and the current versions.

## Version: 1.3.0

* Different statistics are obtained for bending performance.
* The statistics became a part of a list output.
* `small.positive` is not being overwritten anymore if it is greater than the smallest positive eigenvalue.

## Version: 1.3.1

* Bug fixed: Matrices with given colnames and rownames could not be a correlation matrix.

