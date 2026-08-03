program test_singular
   use ppcor, only : dp, ppcor_result, pcor, ppcor_pearson, ppcor_success
   implicit none
   real(dp) :: x(24,4)
   type(ppcor_result) :: r
   integer :: i
   real(dp) :: t

   do i = 1, size(x,1)
      t = real(i,dp)
      x(i,1) = sin(0.2_dp*t) + 0.02_dp*t
      x(i,2) = cos(0.31_dp*t) - 0.01_dp*t
      x(i,3) = x(i,1)
      x(i,4) = sin(0.73_dp*t) + 0.2_dp*cos(0.11_dp*t)
   end do

   call pcor(x, r, ppcor_pearson)
   call assert_true(r%status == ppcor_success, 'singular status')
   call assert_true(r%used_pseudoinverse, 'pseudoinverse flag')
   call assert_true(r%rank == 3, 'rank of duplicated-column matrix')
   call assert_true(index(r%message, 'pseudoinverse') > 0, 'pseudoinverse message')
   call assert_true(all(abs(r%estimate) <= 1.0_dp), 'bounded estimates')

   print '(a)', 'test_singular: PASS'

contains

   subroutine assert_true(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) then
         write(*,'(a)') trim(label)//' failed'
         error stop 1
      end if
   end subroutine assert_true

end program test_singular
