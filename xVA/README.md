# xVA modern Fortran port

This package translates the computational code of the R package **xVA 1.3**
to modern Fortran 2018. It uses the separately structured modern Fortran ports
of **Trading** and **SACCR** as local FPM dependencies.

Implemented calculations include:

- net/gross ratio and CEM exposure at default;
- credit-spread-implied marginal default probabilities;
- discounted CVA/DVA/FVA valuation adjustments;
- Hull-White one-factor simulated IRS exposure profiles;
- collateralized and uncollateralized EE, NEE, PFE, and EEE;
- CEM, IMM, standard SA-CCR, simplified SA-CCR, and OEM EAD;
- effective maturity and IRB default capital;
- BA-CVA, STD-CVA, and SA-CVA capital charges;
- interest-rate and credit-spread SA-CVA sensitivities;
- KVA and SA-CCR-based CVA, DVA, FCA, FBA, and MVA;
- typed loading of all eleven supervisory CSV tables distributed by xVA.

The original package contains no plotting routines. R reference-class, list,
and `data.table` infrastructure is represented by explicit Fortran derived
types and arrays.

## Requirements

- A Fortran 2018 compiler, tested with GNU Fortran 14.2
- FPM for the preferred build, or GNU Make-style shell/batch scripts included
  in this package

## Build with FPM

```text
fpm test
fpm run --target xva_demo
```

Run commands from the package root so the example can find the `data/`
directory.

## Build with GNU Fortran

On Unix-like systems:

```text
./build_gfortran.sh
./build/test_xva
./build/xva_demo
```

On Windows with GNU Fortran available in `PATH`:

```text
build_gfortran.bat
build\test_xva.exe
build\xva_demo.exe
```

## Main module

Applications normally need only:

```fortran
use xva
```

The main types are `simulation_data_t`, `regulatory_data_t`,
`exposure_profile_t`, `cva_sensitivity_t`, `supervisory_cva_data_t`, and
`xva_result_t`. The end-to-end entry point is `xva_calculator`.

See `docs/API_MAP.md` and `docs/PORTING_NOTES.md` for the R-to-Fortran mapping
and deliberate numerical corrections.

## Licensing

The xVA port is GPL-3.0-only, matching the original package. The bundled
Trading and SACCR translations preserve their own GPL-3.0-only notices and
attribution. See `LICENSE`, `NOTICE`, and the dependency directories.
