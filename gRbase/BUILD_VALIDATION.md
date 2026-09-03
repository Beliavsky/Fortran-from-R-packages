# Build validation

Validation was performed from clean temporary directories so no `.o`, `.mod`,
executable, or prior build product inside `gRbase/` could satisfy an interface
by accident.

## Compiler validation

GNU Fortran 14.2.0 was used in two independent clean builds.

Checked build:

```text
gfortran -std=f2018 -O0 -g -fcheck=all -Wall -Wextra \
  -Wimplicit-interface -Werror=implicit-interface -pedantic ...
```

Optimized checked build:

```text
gfortran -std=f2018 -O2 -fcheck=all -Wall -Wextra \
  -Wimplicit-interface -Werror=implicit-interface -pedantic ...
```

Both builds compiled without warnings, all deterministic tests passed, and the
example completed successfully.

Test output:

```text
All gRbase deterministic tests passed.
```

Example output:

```text
fill edges added to 4-cycle: 1
p(variable 1 | variable 2):  0.4000  0.6000  0.3000  0.7000
```

The validation shims used outside the package expose the dependency module
interfaces consumed by this translation (`r_kinds`, `r_linalg`,
`igraph_graph`, `igraph_cliques`, and `igraph_components`). They are not part of
the package or ZIP.

## FPM and formatter availability

The execution environment does not contain Fortran Package Manager or
`fprettify`. The required commands were attempted immediately before final
packaging and the shell returned command-not-found status 127 for each:

```text
fpm build
fpm test
fpm clean --all
fprettify src/*.f90 test/*.f90 example/*.f90
```

Therefore those literal commands could not be truthfully reported as
successful. The package contains an FPM manifest with sibling path dependencies
and was instead validated with the strict direct-gfortran builds described
above. No system BLAS/LAPACK link flags are present.

## Deterministic coverage

Tests exercise:

- cell/index conversion, slices, permutations and combinations;
- set operations and subset enumeration;
- table margins, permutations, equality, expansion modes, alignment,
  list folds, arithmetic, normalization, slicing and deterministic sampling;
- row/column reductions, recycled column multiplication and matrix nonzero
  indices;
- DAG/topological/moralization queries;
- MCS, chordality, elimination/MCWH/minimal triangulation;
- separation, simplicial nodes and ancestry/descendancy queries;
- sibling-igraph clique and component bridge semantics;
- RIP/junction-tree and maximal-prime decomposition;
- SPD inverse and partial-correlation reference values.

## Final source/package audit

Before archive creation the maintained sources are checked for:

- duplicate Fortran files;
- lines over 132 columns and tab characters;
- semicolon-separated Fortran statements;
- `double precision`, `real*8`, `kind(0.0d0)`, and D-exponent literals;
- self-comparison NaN tests;
- dummy arguments lacking explicit `INTENT`/`VALUE`;
- dummy declarations sharing one declaration line;
- dummy declarations lacking meaningful trailing `!!` FORD comments;
- copied dependency sources and build/cache/archive products.
