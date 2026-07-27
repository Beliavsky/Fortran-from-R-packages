# Validation record

## Toolchain

- GNU Fortran 14.2.0
- Fortran 2018 mode
- Checked build: bounds checking, backtraces, interface diagnostics,
  conversion diagnostics, and warnings promoted to errors
- Optimized build: `-O2` with the same warning policy

## Test programs

| Test | Coverage |
|---|---|
| `test_distributions_parameters` | Five emission families and parameter round trips |
| `test_hmm` | Fixed likelihood reference, filtering, smoothing, Viterbi, residuals, deterministic simulation |
| `test_hhmm` | Fixed hierarchical likelihood reference, hierarchical smoothing/decoding/simulation, fitting smoke test |
| `test_diagnostics` | Reordering, corrected and upstream forecasts, AIC/BIC, chunk lengths |
| `test_estimation` | Seeded simulation and ordinary-HMM maximum-likelihood fitting |

Fixed HMM and HHMM likelihood references were calculated independently in
NumPy using direct scaled-forward recursions.

Expected checked-build output:

```text
test_diagnostics: PASS
test_distributions_parameters: PASS
test_estimation: PASS
test_hhmm: PASS
test_hmm: PASS
```

Both checked `-O0` and optimized `-O2` builds passed, together with the demo and
both examples.

## Reproduction

On Unix-like systems:

```text
scripts/validate.sh
```

On Windows with GNU Fortran available in `PATH`:

```text
scripts\validate.bat
```

FPM commands:

```text
fpm build
fpm test
fpm run fhmm_demo
fpm run --example basic_hmm
fpm run --example hierarchical_hmm
```

FPM was not installed in the translation environment. The manifest was parsed
as TOML and the project was validated using the same automatic `src`, `test`,
`app`, and `example` directory structure that FPM discovers.
