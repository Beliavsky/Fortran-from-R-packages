program test_rng
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use gensa_kinds, only : dp, i8
   use gensa_rng, only : ran2_state, gensa_visit
   implicit none

   type(ran2_state) :: a, b
   real(dp), parameter :: expected(5) = [ &
      0.2166913633322185_dp, 0.4993763041901336_dp, &
      0.9153339019060981_dp, 0.1408293410066953_dp, &
      0.7214695212081583_dp ]
   real(dp) :: u, z, step
   integer :: i

   call a%seed(-100377_i8)
   call b%seed(-100377_i8)
   do i = 1, 5
      u = a%uniform()
      call assert_true(abs(u - expected(i)) < 2.0e-15_dp, 'ran2 reference sequence')
      call assert_true(abs(u - b%uniform()) < 1.0e-16_dp, 'ran2 reproducibility')
      call assert_true(u > 0.0_dp .and. u < 1.0_dp, 'uniform range')
   end do

   do i = 1, 1000
      z = a%normal()
      step = gensa_visit(2.62_dp, 10.0_dp, a)
      call assert_true(ieee_is_finite(z), 'normal is finite')
      call assert_true(ieee_is_finite(step), 'visiting deviate is finite')
   end do

   print '(a)', 'PASS test_rng'

contains

   subroutine assert_true(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) then
         write(*, '(a)') 'FAIL: ' // label
         error stop 1
      end if
   end subroutine assert_true

end program test_rng
