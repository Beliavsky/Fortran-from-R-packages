# API coverage

Coverage basis: distinct exported computational R functions in the upstream
`NAMESPACE`. The five exports are all computational and all have mappings.

## abind

`abind_arrays` accepts an array of `array_value` objects. It implements the core
R dimension rules: inputs may have the final rank or one fewer dimension; a
fractional or outside `along` inserts a singleton dimension; all non-join
shapes must conform; data are copied in native column-major order. Optional
argument names, hierarchical join labels, first/last non-join label selection,
`rev_along`, and dimension-name propagation are supported.

Not translated: R list/data-frame coercion, `force.array=FALSE` return types,
call-expression-derived anonymous names, and the `new.names` postprocessing API.

## asub

`asub_array` applies one-based selectors to arbitrary dimensions in any order.
An unallocated or zero-length selector means the whole dimension. Singleton
results are dropped by default, matching ordinary R array subsetting; callers
can request rank preservation with `drop=.false.`. Labels are subset with data.

The R function can accept arbitrary R index objects. The Fortran routine uses
resolved integer positions. `indices_from_names` and `indices_from_mask` provide
common conversions.

## afill<-

`afill_array` accepts one selector per destination dimension. Empty selectors
correspond to RHS dimensions and are resolved by matching RHS dimension labels
to destination labels. Explicit destination dimensions replicate the RHS while
preserving its multidimensional layout. `excess_ok` discards unmatched RHS
labels when enabled.

R's replacement syntax, unevaluated empty arguments, and `local=FALSE` environment
mutation have no direct Fortran analogue and are not reproduced.

## adrop

`adrop_array` validates that every requested drop dimension has length one,
removes exactly those dimensions, and propagates labels. The optional
`named_vector` and `one_d_array` flags preserve the meaningful metadata aspects
of the R behavior, while `array_value` always uses explicit dimension metadata.

## acorn

`acorn_array` generalizes `n`, `m`, `r`, and `...` to a signed integer vector.
Positive counts select from the beginning; negative counts select from the end.
Unspecified later dimensions default to one slice. With `addrownums` enabled,
previously unlabeled output dimensions receive labels such as `[3]` and `[4]`.
