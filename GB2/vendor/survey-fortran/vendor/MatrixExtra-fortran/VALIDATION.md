# Validation

The translated code is compiled with:

```text
-std=f2018 -O0 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all
```

Regression programs compare sparse results with equivalent dense Fortran
calculations for:

1. CSR/CSC/COO/sparse-vector conversions and arbitrary slicing.
2. Sparse assignment plus rbind/cbind.
3. Sparse/dense/sparse matrix and vector multiplication, crossprod/tcrossprod.
4. Elementwise arithmetic, pattern AND/OR, R-style vector recycling and utilities.
5. Norms, diagonals and opposite-storage transpose conversion.
6. Callback-based sparse filtering/mapping.

The upstream R package could not be executed in the validation environment
because R is not installed; therefore validation is against exact dense algebra
semantics rather than an R runtime cross-check.

Final strict validation result: **6/6 test programs passed**, and both examples
compiled and ran successfully. The GNU linker may emit an executable-stack note
for `test_map_filter` because that test passes internal procedures as callbacks;
the translated library itself contains no such internal callback and compiles
without warnings.
