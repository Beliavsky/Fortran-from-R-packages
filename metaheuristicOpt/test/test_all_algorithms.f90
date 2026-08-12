program test_all_algorithms
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use metaheuristic_opt, only : dp, mh_control, mh_result, metaopt
   implicit none
   character(len=8), parameter :: algs(21) = [character(len=8) :: &
      'PSO','ALO','GWO','DA','FFA','GA','GOA','HS','MFO','SCA','WOA', &
      'CLONALG','DE','SFL','CSO','ABC','KH','CS','BA','GBS','BHO']
   type(mh_control) :: c
   type(mh_result) :: r
   real(dp) :: lo(4), hi(4)
   integer :: i

   lo = -5.0_dp
   hi = 5.0_dp
   c%num_population = 24
   c%max_iter = 30
   c%legacy_quirks = .false.
   c%selection_size = 6
   c%num_memeplex = 6
   c%cdc = 2
   c%smp = 5

   do i = 1, size(algs)
      c%seed = 987654 + i
      call metaopt(trim(algs(i)), sphere_local, lo, hi, r, c)
      if (ieee_is_nan(r%value)) error stop 'NaN result'
      if (any(r%par < lo - 1.0e-10_dp)) error stop 'result below lower bound'
      if (any(r%par > hi + 1.0e-10_dp)) error stop 'result above upper bound'
      if (r%value > 100.0_dp) error stop 'unexpectedly bad result'
      if (r%evaluations <= 0) error stop 'evaluation count not updated'
   end do
contains
   real(dp) function sphere_local(x) result(f)
      real(dp), intent(in) :: x(:)
      f = sum(x*x)
   end function sphere_local
end program test_all_algorithms
