# Validation

Validated with GNU Fortran 14.2 using:

```text
-std=f2008
-Wall -Wextra
-Wimplicit-interface -Werror=implicit-interface
-fcheck=all
-ffpe-trap=invalid,zero,overflow
```

The regression suite covers:

1. probability identities and finite-grid normalization;
2. simulation moments;
3. BP and both BZIP EM fits;
4. BNB/BZINB EM fitting, rho calculations, and covariance estimates;
5. direct checks of the translated C++ expectation/score kernel against the
   public likelihood for the five continuous parameters;
6. source-specific mixture score conventions and historical-maximum EM return;
7. exact p4 covariance reconstruction from the 8-parameter information matrix;
8. all-zero fit shortcuts;
9. source-compatible inverse digamma;
10. weighted Pearson correlation;
11. compact and full pairwise BZINB fitting, including pair subsampling.

The transformed-parameter direct likelihood optimizers from v0.1 are retained
and exercised indirectly as reference/fallback APIs, but the public BNB/BZINB
fitters now use the translated specialized EM engine by default.

Simulation is deterministic in the tests through `set_bzinb_seed()`.
