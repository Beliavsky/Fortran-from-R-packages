program basic_randompack
   use randompack, only : dp, randompack_rng, randompack_rng_type
   implicit none
   type(randompack_rng_type) :: rng
   real(dp) :: x(5)
   integer :: info

   rng = randompack_rng('pcg64', seed=12345)
   call rng%normal(x, info=info)
   if (info /= 0) error stop 'normal generation failed'
   print '(a,5(1x,f10.5))', 'normal draws:', x
end program basic_randompack
