! SPDX-License-Identifier: GPL-2.0-or-later
program test_trlan
   use rspectra, only : dense_operator, make_dense_operator
   use svd
   implicit none

   real(dp) :: a(4, 3), s(4, 4), warm_values(1), warm_vectors(4, 1)
   type(dense_operator) :: aop, sop
   type(trlan_svd_result) :: sdense, sop_result
   type(trlan_eigen_result) :: edense, eop_result

   a = 0.0_dp
   a(1, 1) = 4.0_dp
   a(2, 2) = 3.0_dp
   a(3, 3) = 2.0_dp
   aop = make_dense_operator(a)

   sdense = trlan_svd_dense(a, 2)
   call check(sdense%info == 0, "dense trlan.svd status")
   call check(maxval(abs(sdense%d - [4.0_dp, 3.0_dp])) < 1.0e-8_dp, "dense trlan.svd singular values")

   sop_result = trlan_svd_operator(aop, 2)
   call check(sop_result%info == 0, "operator trlan.svd status")
   call check(maxval(abs(sop_result%d - [4.0_dp, 3.0_dp])) < 1.0e-8_dp, "operator trlan.svd singular values")

   s = 0.0_dp
   s(1, 1) = 5.0_dp
   s(2, 2) = 3.0_dp
   s(3, 3) = 1.0_dp
   s(4, 4) = -2.0_dp
   sop = make_dense_operator(s)

   edense = trlan_eigen_dense(s, 2)
   call check(edense%info == 0, "dense trlan.eigen status")
   call check(maxval(abs(edense%d - [5.0_dp, 3.0_dp])) < 1.0e-9_dp, "dense trlan.eigen values")

   warm_values = [5.0_dp]
   warm_vectors = 0.0_dp
   warm_vectors(1, 1) = 1.0_dp
   eop_result = trlan_eigen_operator(sop, 2, lambda_start=warm_values, u_start=warm_vectors)
   call check(eop_result%info == 0, "warm operator trlan.eigen status")
   call check(maxval(abs(eop_result%d - [5.0_dp, 3.0_dp])) < 1.0e-8_dp, "warm operator trlan.eigen values")

   print '(a)', 'test_trlan: PASS'

contains

   subroutine check(condition, label)
      logical, intent(in) :: condition !! Assertion condition that must be true.
      character(len=*), intent(in) :: label !! Short assertion label printed on failure.

      if (.not. condition) then
         print '(a,1x,a)', 'FAIL:', trim(label)
         error stop 1
      end if
   end subroutine check

end program test_trlan
