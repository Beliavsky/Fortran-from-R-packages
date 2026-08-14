program basic
   use, intrinsic :: iso_fortran_env, only : int64, real64
   use rngwell, only : well_rng
   implicit none
   type(well_rng) :: rng
   real(real64) :: x(5)

   call rng%init('19937c', seed=12345_int64)
   call rng%fill(x)
   print '(5(f12.9,1x))', x
end program basic
