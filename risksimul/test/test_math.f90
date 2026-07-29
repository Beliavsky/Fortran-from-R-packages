! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
program test_math
   use risksimul
   use ghyp_special, only : student_cdf
   implicit none
   real(dp) :: x, p
   real(dp), allocatable :: q(:,:), direction(:)
   type(allocation_result) :: alloc
   logical :: ok

   p = 0.975_dp
   x = student_quantile(p,7.0_dp)
   call assert_close(student_cdf(x,7.0_dp),p,2.0e-12_dp)

   x = gamma_quantile(0.9_dp,3.5_dp,1.7_dp)
   call assert_close(gamma_cdf(x,3.5_dp,1.7_dp),0.9_dp,2.0e-12_dp)

   direction = [1.0_dp,2.0_dp,-1.0_dp]
   call orthogonal_completion(direction,q,ok)
   call assert_true(ok)
   call assert_close(maxval(abs(matmul(transpose(q),q)-identity(3))),0.0_dp,2.0e-13_dp)
   call assert_true(dot_product(q(:,1),direction) > 0.0_dp)

   alloc = optimal_allocation_heuristic(reshape([4.0_dp,1.0_dp,1.0_dp,4.0_dp],[2,2]))
   call assert_true(alloc%ok)
   call assert_close(sum(alloc%fractions),1.0_dp,2.0e-14_dp)
   call assert_close(alloc%fractions(1),0.5_dp,1.0e-5_dp)

   print '(a)', 'test_math: PASS'
contains
   pure function identity(n) result(a)
      integer, intent(in) :: n
      real(dp) :: a(n,n)
      integer :: i
      a = 0.0_dp
      do i = 1, n
         a(i,i) = 1.0_dp
      end do
   end function identity

   subroutine assert_true(condition)
      logical, intent(in) :: condition
      if (.not. condition) error stop 1
   end subroutine assert_true

   subroutine assert_close(actual, expected, tolerance)
      real(dp), intent(in) :: actual, expected, tolerance
      if (abs(actual-expected) > tolerance) then
         print '(a,3es25.16)', 'mismatch: ',actual,expected,abs(actual-expected)
         error stop 1
      end if
   end subroutine assert_close
end program test_math
