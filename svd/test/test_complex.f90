! SPDX-License-Identifier: GPL-2.0-or-later
program test_complex
   use svd
   implicit none

   complex(dp) :: a(3, 2)
   type(ztrlan_svd_result) :: result
   real(dp) :: gram(2, 2)

   a = cmplx(0.0_dp, 0.0_dp, dp)
   a(1, 1) = cmplx(3.0_dp, 4.0_dp, dp)
   a(2, 2) = cmplx(2.0_dp, 0.0_dp, dp)

   result = ztrlan_svd_complex(a, 2)
   call check(result%info == 0, "complex status")
   call check(maxval(abs(result%d - [5.0_dp, 2.0_dp])) < 1.0e-10_dp, "complex singular values")
   gram = real(matmul(conjg(transpose(result%u)), result%u), dp)
   call check(maxval(abs(gram - identity2())) < 1.0e-10_dp, "complex left orthogonality")

   print '(a)', 'test_complex: PASS'

contains

   pure function identity2() result(a2)
      real(dp) :: a2(2, 2)

      a2 = 0.0_dp
      a2(1, 1) = 1.0_dp
      a2(2, 2) = 1.0_dp
   end function identity2

   subroutine check(condition, label)
      logical, intent(in) :: condition !! Assertion condition that must be true.
      character(len=*), intent(in) :: label !! Short assertion label printed on failure.

      if (.not. condition) then
         print '(a,1x,a)', 'FAIL:', trim(label)
         error stop 1
      end if
   end subroutine check

end program test_complex
