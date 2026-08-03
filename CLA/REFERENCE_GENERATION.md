# Reference generation

The fixed three-asset reference uses the example covariance matrix and expected
returns from `man/CLA.Rd`. An independent NumPy implementation followed the
current `R/CLA.R` equations for initialization, matrix solves, lambda events,
and free-weight updates. Double precision produced the values recorded in
`test/test_cla_reference.f90`.

The interpolation reference at target mean `0.04` is obtained by linearly
interpolating the two adjacent turning-point weight vectors and recomputing
portfolio variance from the covariance matrix. It gives:

- sigma: `0.080300141561153`
- weights: `(0.987460815047022, 0, 0.012539184952978)`
