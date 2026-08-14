# rngWELL-fortran

Modern Fortran translation of the computational code in R package `rngWELL`
version 0.10-10, organized as an FPM package.

The library implements all 17 WELL generators contained in the supplied
package:

- WELL512a
- WELL521a, WELL521b
- WELL607a, WELL607b
- WELL800a, WELL800b
- WELL1024a, WELL1024b
- WELL19937a, WELL19937b, WELL19937c
- WELL21701a
- WELL23209a, WELL23209b
- WELL44497a, WELL44497b

It also implements the MT2002-style scalar seed expansion used by rngWELL,
raw 32-bit output, uniform real output, vector/matrix filling, canonical state
export/import, and an option mapper corresponding to the R `WELL2test()`
order/version/tempering choices.

## Build

```text
fpm build
fpm test
fpm run --example basic
```

The project has no R or C runtime dependency.

## Basic use

```fortran
program demo
   use, intrinsic :: iso_fortran_env, only : int64, real64
   use rngwell, only : well_rng
   implicit none
   type(well_rng) :: rng
   real(real64) :: x(10)

   call rng%init('19937c', seed=12345_int64)
   call rng%fill(x)
   print *, x
end program demo
```

`next_uint32()` returns the raw signed `int32` bit pattern. It should be
interpreted as an unsigned 32-bit word when comparing hexadecimal reference
streams. `next()` returns the same word multiplied by 2^-32 as `real64`.

## State handling

`get_state()` returns the state in the same canonical cyclic order as the
upstream `GetWELLRNG*()` functions. Passing that state to `init(..., state=)`
or `put_state()` resumes the stream exactly.

## Licensing

The supplied R package declares BSD-3-Clause licensing; its original CRAN
license metadata is retained under `LICENSES/rngWELL-CRAN-LICENSE.txt`.
The bundled WELL C sources additionally carry a source notice restricting
commercial use. That notice is preserved in `LICENSES/WELL-SOURCE-NOTICE.txt`.
Users should review both before redistribution or commercial use.
