program test_de_strategies
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use metaheuristic_opt, only : dp, mh_control, mh_result, de
   implicit none
   character(len=16), parameter :: s(6) = [character(len=16) :: &
      'clasical','best 1','target to best','best 2','rand 2','rand 2 dir']
   type(mh_control) :: c
   type(mh_result) :: r
   real(dp) :: lo(4), hi(4)
   integer :: i

   lo = -5.0_dp
   hi = 5.0_dp
   c%num_population = 30
   c%max_iter = 50
   c%legacy_quirks = .false.
   do i = 1, size(s)
      c%seed = 8000+i
      c%de_strategy = s(i)
      call de(sphere_local, lo, hi, r, c)
      if (ieee_is_nan(r%value)) error stop 'DE produced NaN'
      if (any(r%par < lo) .or. any(r%par > hi)) error stop 'DE result out of bounds'
      if (r%value > 20.0_dp) error stop 'DE failed gross regression'
   end do
contains
   real(dp) function sphere_local(x) result(f)
      real(dp), intent(in) :: x(:)
      f = sum(x*x)
   end function sphere_local
end program test_de_strategies
