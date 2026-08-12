program test_dispatch_seed_max
   use metaheuristic_opt, only : dp, mh_control, mh_result, metaopt, metaopt_many
   implicit none
   type(mh_control) :: c
   type(mh_result) :: r1, r2, r3, many(2)
   real(dp) :: lo(3), hi(3)

   lo = -4.0_dp
   hi = 4.0_dp
   c%num_population = 30
   c%max_iter = 80
   c%seed = 20260810
   c%legacy_quirks = .false.

   call metaopt('GWO', sphere_local, lo, hi, r1, c)
   call metaopt('gwo', sphere_local, lo, hi, r2, c)
   call metaopt_many([character(len=3) :: 'GWO', 'PSO'], sphere_local, lo, hi, many, c)
   if (maxval(abs(r1%par-r2%par)) > 1.0e-14_dp) error stop 'seed reproducibility failed'
   if (abs(r1%value-r2%value) > 1.0e-14_dp) error stop 'dispatcher case handling failed'
   if (abs(r1%value-many(1)%value) > 1.0e-14_dp) error stop 'metaopt_many seed reset failed'

   c%maximize = .true.
   c%seed = 1111
   call metaopt('SCA', peak, lo, hi, r3, c)
   if (r3%value < -5.0e-2_dp) error stop 'maximization failed'
   if (maxval(abs(r3%par-1.0_dp)) > 0.30_dp) error stop 'maximizer location poor'
contains
   real(dp) function sphere_local(x) result(f)
      real(dp), intent(in) :: x(:)
      f = sum(x*x)
   end function sphere_local

   real(dp) function peak(x) result(f)
      real(dp), intent(in) :: x(:)
      f = -sum((x-1.0_dp)**2)
   end function peak
end program test_dispatch_seed_max
