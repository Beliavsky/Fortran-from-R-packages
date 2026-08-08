# NMOF integration demo

This secondary FPM project demonstrates direct use of a `neighbours`
numeric neighbourhood as the callback expected by the supplied NMOF
Fortran translation's `local_search` routine.

From this directory:

```text
fpm run
```

Its dependencies are local paths to the parent `neighbours` package and the
vendored copy of the user-supplied NMOF translation.
