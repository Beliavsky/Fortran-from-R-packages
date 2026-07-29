# Independent reference generation

The fixed numerical values in `test/test_moments.f90` were generated with
`tools/generate_reference.py`.

The script does not call the Fortran implementation. It independently:

1. centers each return column,
2. computes sample covariance with divisor `T-1`,
3. constructs third- and fourth-order co-moments with explicit loops,
4. performs Kronecker-ordered portfolio contractions,
5. evaluates the source package's decomposition formulas.

Run with:

```sh
python3 tools/generate_reference.py
```

NumPy is needed only for regenerating references, not for building or using the
Fortran library.
