# Validation

The maintained Fortran sources, deterministic tests, example, and parity driver
were compiled directly with GNU Fortran 14.2.0 using:

```text
-std=f2018 -pedantic -Wall -Wextra -Wimplicit-interface -Werror -fcheck=all -fbacktrace -O0 -g
```

`test/test_dbscan.f90` passes and covers exact kNN/frNN identities and tie
ordering, DBSCAN and core-point classification, LOF reference values,
point-density, SNN/Jarvis-Patrick, OPTICS with DBSCAN/Xi extraction, core and
mutual-reachability distances, MST construction, HDBSCAN flat extraction,
exact GLOSH scores, membership probabilities, the upstream 118-point
HDBSCAN(e) epsilon-selection regression, constrained/mixed FOSC extraction,
DBCV, and structural conversions. `example/basic.f90` also compiles and runs.

The independent `validation/parity_driver.f90` / `compare_sklearn.py` check
passes with scikit-learn 1.8.0 for the deterministic validation data. It checks
DBSCAN and HDBSCAN partitions modulo label numbering, LOF values, and OPTICS
reachability/core distances. HDBSCAN probabilities are checked against the R
package's core-distance rule because scikit-learn uses a different membership
probability definition.

`python validation/audit_source.py` also passes. It checks the maintained
Fortran source for the requested free-form/source-policy constraints, including
dummy-argument `INTENT`/`VALUE` declarations and FORD comments, a single `dp`
definition, forbidden semicolon statements, legacy real kinds and D-exponents,
self-comparison NaN idioms, duplicate Fortran source content, packaged build
products, disallowed system BLAS/LAPACK linkage, and translation-coverage count
consistency.

FPM and fprettify are not installed in the execution environment used for this
translation. Their executables were searched for directly. An attempt to fetch
the current stable FPM release could not proceed because the container has no
external DNS/download access. The direct strict-gfortran build, test, example,
parity, and static-audit path described above was therefore used for executed
validation. All generated compiler/build products are outside the package tree
and are removed before packaging. The package remains structured for the
required `fpm build`, `fpm test`, and `fpm clean --all` commands on a normal FPM
installation.
