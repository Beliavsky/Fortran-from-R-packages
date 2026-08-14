# tabuSearch-fortran

Modern Fortran/FPM translation of the computational code in the R package
`tabuSearch` 1.2.0 by Katarina Domijan.

The package implements tabu search for binary configurations.  Plotting and
R S3 infrastructure are intentionally omitted; the search algorithm and the
computational parts of `summary.tabu` are provided natively in Fortran.

## Implemented functionality

- binary configuration search
- user-supplied objective callback
- random or user-supplied starting configuration
- complete or sampled one-bit neighborhoods
- aspiration criterion
- FIFO tabu list
- preliminary search
- intensification restarts from the best configuration found
- diversification using the most frequent moves
- repeated complete searches
- complete configuration/objective history
- best-value and best-configuration helpers
- computational summary statistics
- explicit reproducible local RNG state

## Basic use

```fortran
program demo
   use tabu_search, only : dp, tabu_control, tabu_result, run_tabu_search
   implicit none

   type(tabu_control) :: control
   type(tabu_result) :: result
   integer :: initial(10)

   initial = 0
   control%iters = 50
   control%neigh = 10
   control%list_size = 3
   control%n_restarts = 4
   control%seed = 12345

   call run_tabu_search(10, objective, result, control, initial)
   print *, result%best_value()
   print *, result%best_configuration()

contains

   function objective(config) result(value)
      integer, intent(in) :: config(:)
      real(dp) :: value

      value = real(sum(config), dp)
   end function objective
end program demo
```

`control%neigh = 0` means use all `size` one-bit neighbors, matching the
upstream default `neigh = size`.

## FPM

```text
fpm build
fpm test
fpm run --example basic_tabu
```

The project intentionally requires modern language semantics:

```toml
[fortran]
implicit-typing = false
implicit-external = false
source-form = "free"
```

There are no external numerical-library dependencies.

## Numerical/semantic corrections

The Fortran code follows the intended algorithm while correcting several
R-storage/evaluation edge cases.  In particular, only sampled neighbors are
eligible for selection, so an unevaluated neighbor cannot acquire an implicit
utility of zero; diversification frequencies use only completed history rows;
and negative objective functions are valid during intensification.  See
`ALGORITHM_NOTES.md` for details.

## License

The upstream package is GPL (>= 2).  This translation is therefore distributed
under GPL-2.0-or-later.  See `LICENSE`, `COPYING.GPL-2`, and `COPYING.GPL-3`.
