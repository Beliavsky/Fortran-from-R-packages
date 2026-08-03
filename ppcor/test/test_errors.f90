program test_errors
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use ppcor, only : dp, ppcor_result, pcor, ppcor_pearson, &
                     ppcor_invalid_shape, ppcor_nonfinite_data, &
                     ppcor_invalid_method, ppcor_constant_variable
   implicit none
   real(dp) :: small(1,2), bad(6,3)
   type(ppcor_result) :: r
   integer :: i

   small = 1.0_dp
   call pcor(small, r, ppcor_pearson)
   call assert_true(r%status == ppcor_invalid_shape, 'invalid shape')

   do i = 1, 6
      bad(i,:) = [real(i,dp), sin(real(i,dp)), cos(real(i,dp))]
   end do
   bad(2,1) = ieee_value(0.0_dp, ieee_quiet_nan)
   call pcor(bad, r, ppcor_pearson)
   call assert_true(r%status == ppcor_nonfinite_data, 'nonfinite rejection')

   call pcor(constant_data(), r, ppcor_pearson)
   call assert_true(r%status == ppcor_constant_variable, 'constant variable')

   call pcor(bad, r, 99)
   call assert_true(r%status == ppcor_invalid_method, 'invalid method')

   print '(a)', 'test_errors: PASS'

contains

   function constant_data() result(a)
      real(dp) :: a(6,3)
      integer :: k
      do k = 1, 6
         a(k,:) = [real(k,dp), 1.0_dp, real(k*k,dp)]
      end do
   end function constant_data

   subroutine assert_true(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) then
         write(*,'(a)') trim(label)//' failed'
         error stop 1
      end if
   end subroutine assert_true

end program test_errors
