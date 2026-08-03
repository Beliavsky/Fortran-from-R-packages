# GenSA modern Fortran

A self-contained modern Fortran/FPM translation of the computational core of
R package **GenSA 1.1.15**. It minimizes bounded nonlinear objective functions
with generalized simulated annealing (GSA), optionally followed by bounded
local optimization.

## Implemented functionality

- Tsallis generalized visiting distribution
- Generalized acceptance probability
- Original GenSA temperature schedule and temperature restart
- Full-vector and single-coordinate Markov proposals
- Periodic wraparound at box bounds
- Optional nonlinear feasibility callback
- Projected BFGS local search for smooth objectives
- Bounded coordinate-pattern local search for nonsmooth objectives
- Objective threshold, function-call, CPU-time, iteration, and
  no-improvement stopping controls
- Reproducible stateful `ran2` random generation
- Trace history with step, temperature, current value, and best value
- Defensive handling of nonfinite objective values

The project has no external numerical dependencies.

## Build with FPM

```text
fpm build
fpm test
fpm run --example rastrigin_example
fpm run demo_gensa
```

The package version is the FPM-compatible numeric value `1.1.15`.

## Minimal example

```fortran
program example
   use gensa
   implicit none

   type(gensa_control) :: control
   type(gensa_result) :: result
   real(dp) :: lower(2), upper(2)

   lower = -5.12_dp
   upper = 5.12_dp
   control%maxit = 1500
   control%has_threshold = .true.
   control%threshold_stop = 1.0e-8_dp

   call gensa_minimize(rastrigin, lower, upper, result, control)
   print *, result%value
   print *, result%par

contains

   function rastrigin(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: value
      value = sum(x*x - 10.0_dp*cos(2.0_dp*acos(-1.0_dp)*x)) &
         + 10.0_dp*real(size(x), dp)
   end function rastrigin

end program example
```

Additional data can be supplied to an objective through host association, a
module procedure, or a derived-type object used by module-level state. This is
the Fortran counterpart of R's `...` arguments.

## Public API

The main entry point is:

```fortran
call gensa_minimize(objective, lower, upper, result, control, initial, constraint)
```

`control`, `initial`, and `constraint` are optional. The objective callback has
signature:

```fortran
function objective(x) result(value)
   use gensa, only : dp
   real(dp), intent(in) :: x(:)
   real(dp) :: value
end function objective
```

A constraint callback returns `.true.` for feasible points.

See `docs/API_MAP.md`, `docs/ALGORITHM.md`, and `docs/PORTING_NOTES.md` for
control mappings and implementation details.

## Validation

The regression suite covers the upstream random generator, visiting deviates,
smooth convex optimization, multimodal Rastrigin optimization, nonsmooth
optimization, nonlinear constraints, trace invariants, and exact call-limit
termination.

The supplied GNU Fortran scripts build all tests, examples, and the demo in
checked or optimized modes.

## License

GNU General Public License version 2. See `LICENSE` and `NOTICE.md`.
