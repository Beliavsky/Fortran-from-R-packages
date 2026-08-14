# Algorithm notes

## Objective convention

The translated nondominated sorting treats smaller objective values as better,
matching the active `ecr::doNondominatedSorting` path in rmoo.

## NSGA-II

The implementation uses Deb et al.'s fast pairwise nondominated sorting and
standard front-wise crowding distance.  Boundary points receive IEEE positive
infinity.  Environmental selection takes complete fronts and then the most
crowded members of the split front.

## NSGA-III

The port implements the numerical steps contained in rmoo directly:

1. update ideal and worst points;
2. determine the split front;
3. compute extreme points with achievement scalarizing functions;
4. construct the hyperplane nadir estimate with a pivoted linear solve and
   fallback to the worst first-front point;
5. normalize objectives and associate them with reference directions;
6. perform Deb-Jain niching for the split front.

The package's use of the raw fitness norm in the final perpendicular-distance
scaling is preserved.

## R-NSGA-II

`rnsga2_survivors` implements reference-point preference ranking, weighted
normalized objective distance, epsilon grouping, and split-front truncation.
The low-level routine accepts the three normalization modes (`NORM_EVER`,
`NORM_FRONT`, `NORM_NONE`).  The high-level optimizer uses the upstream default
front normalization.  The optional upstream behavior of permanently appending
new extreme points to the reference-direction set is not enabled by the
high-level v0.1.0 driver.

## Genetic operators

The default high-level choices follow rmoo's control table:

- real: SBX + polynomial mutation;
- binary: single-point crossover + random bit mutation;
- permutation: ordered crossover + inversion mutation;
- discrete integer: uniform crossover + per-gene integer mutation.

HUX, uniform crossover, linear-rank selection, and random real mutation remain
available as public lower-level routines.

## Intentional corrections

A few clear R implementation defects are not reproduced literally:

- the active vectorized `PerformScalarizing` code computes one scalarizing
  value across the whole weight matrix and can select the same point for every
  objective; the Fortran routine implements the intended objective-specific ASF;
- the later-front branch of `sharing()` can reuse an old `nichecount` because it
  is not reset for each individual; the Fortran routine resets it to one;
- random decisions are never embedded in Fortran `.and.`/`.or.` expressions,
  avoiding optimizer-dependent function-elimination and evaluation-order issues.

These corrections affect implementation defects, not the published NSGA
mathematics.
