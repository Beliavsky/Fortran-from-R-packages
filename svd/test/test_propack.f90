! SPDX-License-Identifier: GPL-2.0-or-later
program test_propack
   use rspectra, only : dense_operator, make_dense_operator
   use svd
   implicit none

   real(dp) :: a(4, 3)
   type(propack_svd_result) :: dense_result, operator_result
   type(dense_operator) :: op

   a = 0.0_dp
   a(1, 1) = 4.0_dp
   a(2, 2) = 3.0_dp
   a(3, 3) = 2.0_dp
   op = make_dense_operator(a)

   dense_result = propack_svd_dense(a, 2)
   call check(dense_result%info == 0, "dense status")
   call check(size(dense_result%d) == 2, "dense singular-value count")
   call check(maxval(abs(dense_result%d - [4.0_dp, 3.0_dp])) < 1.0e-9_dp, "dense singular values")
   call check(maxval(abs(matmul(transpose(dense_result%u), dense_result%u) - identity2())) < 1.0e-9_dp, &
      "dense left orthogonality")
   call check(maxval(abs(matmul(transpose(dense_result%v), dense_result%v) - identity2())) < 1.0e-9_dp, &
      "dense right orthogonality")

   operator_result = propack_svd_operator(op, 2)
   call check(operator_result%info == 0, "operator status")
   call check(maxval(abs(operator_result%d - [4.0_dp, 3.0_dp])) < 1.0e-8_dp, "operator singular values")

   print '(a)', 'test_propack: PASS'

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

end program test_propack
