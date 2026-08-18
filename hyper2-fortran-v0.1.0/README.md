# hyper2-fortran

Modern free-format Fortran translation of the computational code in the R
package `hyper2` 3.2-3.

## Highlights

The library implements the package's sparse likelihood algebra for products of
sums of player strengths (`hyper2`) and weighted sums (`hyper3`). It provides:

- canonical sparse hyper2 and hyper3 term storage;
- model addition/subtraction/scaling/equality and term access/update;
- log likelihoods, exact gradients, and simplex Hessians;
- constrained simplex maximum-likelihood fitting and multistart fitting;
- Dirichlet, generalized-Dirichlet, weighted-Dirichlet, rank, pairwise,
  home/away, draw, and Bradley-Terry/Zermelo constructors;
- repeated-player weighted race likelihoods (`ordervec2supp3`);
- player keeping/discarding/substitution, balance, and PWA transformations;
- grouped-rank likelihood sums (`suplist`) using `partitions-fortran`;
- list-of-sum-likelihood (`lsl`) evaluation and optimization;
- simplex normalizing constants, densities, moments, and MGF calculations using
  `cubature-fortran`;
- rank-table/order-table conversion and support construction;
- Metropolis support sampling, Dirichlet/rank/race simulators, and hyper3 race
  simulation;
- equal-probability and known-probability likelihood-ratio tests.

R classes, printing, data-frame conversion, plotting, package datasets, and
other presentation-only machinery are intentionally omitted. See
`TRANSLATION_NOTES.md` for the detailed mapping and differences.

## Build

```console
fpm build
fpm test
fpm run --example bradley_terry
```

All Fortran source is free-format `.f90`; no C or C++ code is compiled or
called. Original R/Rcpp/C++ files under `upstream/` are provenance only.

## Minimal example

```fortran
program demo
    use hyper2
    implicit none
    type(hyper2_model) :: h
    type(fit_result) :: fit
    real(dp) :: wins(2,2)
    character(len=name_len) :: names(2)

    names = ['a                                                               ', &
             'b                                                               ']
    wins = reshape([0.0_dp, 2.0_dp, 8.0_dp, 0.0_dp], [2,2])
    h = pairwise(wins, names)
    fit = maxp(h, startp=[0.8_dp,0.2_dp])
    print *, fit%p
end program demo
```

## Dependencies

The package vendors:

- `cubature-fortran-v0.1.0` for adaptive simplex integration;
- `partitions-fortran-v0.1.0` for grouped-rank permutation enumeration.

Both dependencies are pure Fortran and are declared as local FPM dependencies.
