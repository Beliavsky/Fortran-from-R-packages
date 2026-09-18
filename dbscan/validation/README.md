# Independent validation

`parity_driver.f90` exports a fixed deterministic data set through the Fortran
implementations of DBSCAN, LOF, HDBSCAN, and OPTICS.  `compare_sklearn.py`
compares the results with scikit-learn where definitions agree.

HDBSCAN cluster labels are compared modulo label numbering.  R `dbscan` uses a
core-distance membership-probability formula that differs from scikit-learn,
so the script verifies that formula independently from nearest-neighbor core
distances.  OPTICS reachability/core distances are compared by original point;
ordering is not compared because the two implementations use different tie
rules.
