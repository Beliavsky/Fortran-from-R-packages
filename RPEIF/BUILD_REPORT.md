# Build report

Translation version: 0.1.0

Compiler used for validation:

```text
GNU Fortran 14.2.0
```

Checked flags:

```text
-std=f2018 -O0 -g -Wall -Wextra -Wpedantic -Wconversion-extra
-fcheck=all -fbacktrace -ffree-line-length-none
```

Optimized flags:

```text
-std=f2018 -O3 -Wall -Wextra -Wpedantic -ffree-line-length-none
```

Results:

```text
test_dispatch_shape: PASS
test_influence_core: PASS
test_nuisance_stats: PASS
test_robust_prewhiten: PASS
test_tail_ratios: PASS
```

The example `influence_demo` also compiled and ran successfully.

The validation environment did not contain an FPM executable. The package was
therefore compiled directly with GNU Fortran using the same module graph and
flags represented by the FPM manifest. The final package includes an FPM
manifest and path dependency suitable for `fpm build` and `fpm test`.

The clean ZIP archive was independently extracted. Both checked and optimized
Makefile builds, all five tests, and the example executable passed from that
extracted copy.
