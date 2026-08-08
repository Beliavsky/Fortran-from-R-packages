# n1qn1 6.0.1-14

* As requested by CRAN, gcc-asan and valgrind

* Make package thread-safe: convert global state to `thread_local`, remove
  `static` from local variables in Fortran-translated C code, add integer
  overflow guards, and fix memory leaks in the R callback wrappers.

* Review all C integer types; confirmed existing types are correct given
  Fortran ABI constraints (`int*` required for Fortran INTEGER) and R's
  `INTEGER()` returning `int*`.

* Fix undefined behavior: add virtual destructor to `EvalBase` so deleting
  derived callback objects through a base pointer is well-defined.

* Fix memory leak: validate `n` and `nzm` before allocating callback objects,
  preventing a leak when `Rf_error` is called on bad inputs.

* Fix `restart=TRUE`: typo caused mode to remain 2 instead of being set to 3.

* Fix `assign=TRUE`: wrong field name meant the compressed Hessian was never
  stored in the supplied environment.

* Add tests for restart, assign, single-variable optimization, repeated calls,
  and input error conditions.

# n1qn1 6.0.1-13

* Make package thread-safe: convert global state to `thread_local`, remove
  `static` from local variables in Fortran-translated C code, add integer
  overflow guards, and fix memory leaks in the R callback wrappers.

* Review all C integer types; confirmed existing types are correct given
  Fortran ABI constraints (`int*` required for Fortran INTEGER) and R's
  `INTEGER()` returning `int*`.

# n1qn1 6.0.1-12

* Add non binary (function pointer) interface/api so that `nlmixr2est`
  will not have to be re submitted when changes to this package occur.

# n1qn1 6.0.1-11

* Added strict prototype fixes as requested by CRAN

* Added a `NEWS.md` file to track changes to the package.
