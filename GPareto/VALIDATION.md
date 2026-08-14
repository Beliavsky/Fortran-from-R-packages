# Validation

Validation was performed with GNU Fortran 14.2 using strict free-form Fortran
2018 checks, explicit-interface warnings as errors, and runtime bounds/allocation
checking.

Permanent tests cover:

1. Pareto filtering, 2D/3D exact hypervolume, distance semantics and benchmarks.
2. Analytical 2D EHI, including deterministic zero-variance behavior.
3. Exact 2D and 3D probability of non-domination.
4. DiceKriging model fitting/prediction, EHI, qEHI, SUR and model update.
5. CPF/Vorob'ev calculations.
6. Differential-evolution optimization.
7. Sequential `gparetoptim` and conditional Pareto-set density workflow.

Independent development differential tests:

- 250 random hypervolume cases (2-4 objectives, up to 6 points) versus an
  independent inclusion-exclusion implementation: maximum absolute difference
  `1.9984e-15`.
- 200 random analytical two-objective EHI cases versus an independent direct
  implementation of the upstream C++ cell formula: zero difference at printed
  double precision.
- 200 random two/three-objective probability-of-nondomination cases versus an
  independent Gaussian inclusion-exclusion calculation: maximum absolute
  difference `4.1633e-16`.

A fixed analytical EHI regression value used in the permanent tests is
`0.209242250192` for front `(1,4),(2,2),(4,1)`, reference `(5,5)`, predictive
mean `(2.5,2.5)` and standard deviations `(0.5,0.7)`.
