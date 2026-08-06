program test_wilson_kernel
   use smith_wilson_kinds, only : dp
   use smith_wilson, only : create_kernel_matrix, sw_success, wilson_function
   implicit none

   real(dp), allocatable :: w(:, :)
   real(dp) :: times(2)
   integer :: failures, info
   character(len=256) :: message

   failures = 0
   times = [5.0_dp, 20.0_dp]

   call check_close(wilson_function(0.0_dp, 10.0_dp, 0.042_dp, 0.1_dp), 0.0_dp, &
                    1.0e-15_dp, 'kernel is zero at time zero')
   call check_close(wilson_function(5.0_dp, 20.0_dp, 0.042_dp, 0.1_dp), &
                    wilson_function(20.0_dp, 5.0_dp, 0.042_dp, 0.1_dp), &
                    1.0e-15_dp, 'kernel symmetry')

   call create_kernel_matrix(times, w, 0.042_dp, 0.1_dp, info, message)
   call check(info == sw_success, 'kernel matrix status')
   call check(size(w, 1) == 2 .and. size(w, 2) == 2, 'kernel matrix shape')
   call check_close(w(1, 1), 0.12085700844851822_dp, 2.0e-15_dp, 'W(5,5)')
   call check_close(w(1, 2), 0.15029036138212082_dp, 2.0e-15_dp, 'W(5,20)')
   call check_close(w(2, 2), 0.28126774328081260_dp, 2.0e-15_dp, 'W(20,20)')

   if (failures /= 0) error stop 1
   print '(a)', 'test_wilson_kernel: PASS'

contains

   subroutine check(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label

      if (.not. condition) then
         failures = failures + 1
         print '(a)', 'FAIL: '//trim(label)
      end if
   end subroutine check

   subroutine check_close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual, expected, tolerance
      character(len=*), intent(in) :: label

      call check(abs(actual - expected) <= tolerance, label)
   end subroutine check_close

end program test_wilson_kernel
