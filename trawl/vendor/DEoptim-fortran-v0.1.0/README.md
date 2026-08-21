# DEoptim-fortran

Modern Fortran translation of the computational core of R package `DEoptim`
2.2-8, packaged for FPM.

## Implemented

- Differential Evolution minimization over box constraints.
- Strategy 1: DE/rand/1/bin.
- Strategy 2: DE/local-to-best/1/bin (DEoptim default).
- Strategy 3: DE/best/1/bin with jitter.
- Strategy 4: DE/rand/1/bin with per-vector dither.
- Strategy 5: DE/rand/1/bin with per-generation dither.
- Strategy 6: DE/current-to-p-best/1, including DEoptim's optional adaptive
  crossover/weight update controlled by `c`.
- Random bound repair used by the DEoptim C engine.
- User-supplied initial populations.
- `VTR`, `reltol`, `steptol`, `NP`, `F`, `CR`, `p`, iteration limits, and trace
  controls.
- Best-member/value histories, final population, optional stored populations,
  and function-evaluation counts.
- `fnMap`-style mapping callback. As in the R wrapper, duplicate mapped members
  are detected and regenerated/remapped for up to five rounds.
- Portable standalone pseudo-random generator with uniform, normal, and Cauchy
  variates, and a user-settable seed.

## Intentionally omitted

- R `.Call`, S3 methods, `...` argument marshalling, package loading, and other
  R-specific infrastructure.
- Plotting code.
- R `parallel` / `foreach` orchestration. The numerical engine accepts a
  Fortran objective callback; callers are free to parallelize expensive work
  inside that callback or build a batch wrapper around the library.
- `bs=.true.`. Although old DEoptim documentation describes it, the active
  DEoptim 2.2-8 C engine explicitly stops with `bs = TRUE not currently
  supported`. The Fortran port returns `de_unsupported` for the same option.

## Basic use

```fortran
program example
    use deoptim, only : dp, i8, de_control, de_result, deoptim_solve
    implicit none
    type(de_control) :: control
    type(de_result) :: result
    real(dp) :: lower(2), upper(2)

    lower = -10.0_dp
    upper =  10.0_dp
    control%np = 50
    control%itermax = 400
    control%seed = 21_i8

    call deoptim_solve(rosenbrock, lower, upper, result, control)
    print *, result%bestval, result%bestmem
contains
    function rosenbrock(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp) :: f
        f = 100.0_dp*(x(2)-x(1)**2)**2 + (1.0_dp-x(1))**2
    end function rosenbrock
end program example
```

The objective callback can be an internal procedure, so additional problem
parameters can be supplied naturally through host association.

## Population layout

The public Fortran API uses the natural row-oriented shapes

- population: `(NP, npar)`
- best-member history: `(npar, iter)`
- stored populations: `(NP, npar, nstore)`

This avoids reproducing the R/C column-major marshalling layer.

## Random numbers

The R package uses R's global RNG. A standalone FPM library cannot use that RNG
without depending on R, so this port supplies a portable Park-Miller generator
and Box-Muller/Cauchy transforms. Set `control%seed` to a nonzero integer for
reproducible Fortran runs. The same seed is not expected to reproduce R's exact
population sequence, but the DE algorithms and random decisions are the same.

## Build

With FPM installed:

    fpm build
    fpm test
    fpm run --example rosenbrock

The code is Fortran 2018 and has no external numerical-library dependency.

## License

GPL-2.0-or-later. See `LICENSES.md`, `LICENSE-GPL-2`, and `LICENSE-GPL-3`.
