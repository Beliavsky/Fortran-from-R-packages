# Optional mvtnorm companion

The user supplied a modern Fortran translation of `mvtnorm` alongside FNN.
Inspection of FNN 1.1.4.1 shows that `mvtnorm` is listed under `Suggests`, not
`Depends` or `Imports`.  The only reference is an example in `man/knn.dist.Rd`
that uses `rmvnorm` to generate correlated multivariate-normal data.

Therefore `FNN-fortran` deliberately has no `mvtnorm` dependency.

To reproduce that style of example in a larger FPM application, add both
projects as dependencies and generate the sample with `mvtnorm`'s `rmvnorm`,
then pass the resulting `n x p` matrix directly to `get_knn` or `knn_dist`.
No adapter layer is required because both ports represent observations by rows.
