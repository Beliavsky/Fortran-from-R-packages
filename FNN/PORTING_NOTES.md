# Porting notes

## Scope

This is a numerical/algorithmic translation of FNN 1.1.4.1.  The original
package contains R wrappers plus two historical C++ nearest-neighbor libraries:
ANN 1.1.2 and a cover-tree implementation.  The Fortran port replaces those
pointer-heavy/template implementations with native typed Fortran structures
while preserving exact nearest-neighbor results.

No plotting routines exist in the supplied FNN package.  S3 printing and R
factor/model metadata are interface concerns and were omitted.

## kd-tree

The kd-tree recursively chooses the coordinate with the largest spread,
partitions around a median, and performs exact split-plane branch-and-bound
search.  A deterministic `(coordinate,index)` ordering resolves ties.  Search
uses Euclidean distance and `epsilon = 0`, corresponding to exact ANN queries.

## cover hierarchy

The cover-tree backend uses the upstream scale base 1.3.  Points are inserted
into a scaled covering hierarchy; exact maximum distance from every node to
its descendants is then computed.  Queries use these subtree radii as rigorous
lower bounds, making pruning exact even when the translated hierarchy differs
in layout from the historical C++ template implementation.

## duplicate observations

ANN's original `get_KNN_kd/get_KNN_brute` wrappers ask for `k+1` results and
blindly discard the first result, assuming it is the query observation.  With
duplicated coordinates a different zero-distance observation can be first.
The Fortran implementation excludes the query by its index instead.  This is
an intentional correctness fix.

## correlation distance

`"CR"` preserves the package's literal metric

```text
1 - dot_product(x, y)
```

and does not normalize the vectors.  Callers wanting cosine distance should
normalize rows themselves, just as with the original implementation.

## information estimators

The C++ KL kernels store squared ANN distances and divide log-distance terms
by two.  The Fortran neighbor API stores ordinary Euclidean distances, so the
same formulas are implemented directly with `log(distance)` and no factor of
one half.

A dependency-free asymptotic/recurrence implementation of the digamma function
is included so the estimators do not require a special-function library.

## mutual information

`mutinfo` implements the multidimensional KSG route in `mdmutinfo`: the joint
neighbor radius is the maximum norm over all X and Y coordinates, while the
marginal counts use strict `< radius` comparisons and include the observation
itself when the radius is positive.

## OWNN automatic k

The original R code applies five-fold labels to columns of an `n x n` matrix
whose columns actually represent neighbor ranks, not observation identities.
The translated automatic-k routine implements the intended five-fold
cross-validation directly: training observations from the held-out fold are
removed, predictions are made for that fold, and classification errors are
accumulated across the 21 candidate k values.  A seed controls a deterministic
fold shuffle.

The OWNN `kstar` and BNN probability parameter are also bounded to dimensions
that yield valid finite arrays.  These are robustness corrections to edge
cases; ordinary interior cases use the upstream formulas unchanged.

## mvtnorm

`mvtnorm` is not an imported dependency of FNN.  It appears under `Suggests`
and is used only in documentation examples.  No core source file calls it, so
linking it into this library would introduce an unnecessary dependency.
