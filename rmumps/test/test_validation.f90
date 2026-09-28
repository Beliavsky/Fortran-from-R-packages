program test_validation
   use rmumps_api, only : dp, rmumps_matrix, rmumps_init_triplet, rmumps_solve, rmumps_det, &
      rmumps_determinant, rmumps_set_values, rmumps_set_permutation, rmumps_get_sym
   implicit none

   type(rmumps_matrix) :: a
   type(rmumps_matrix) :: bad
   integer :: info
   integer :: i5(5)
   integer :: j5(5)
   integer :: i4(4)
   integer :: j4(4)
   integer :: i3(3)
   integer :: j3(3)
   integer :: ri(1)
   integer :: rj(1)
   real(dp) :: v5(5)
   real(dp) :: v4(4)
   real(dp) :: v3(3)
   real(dp) :: rv(1)
   real(dp) :: b3(3)
   real(dp) :: x3(3)
   real(dp) :: xbad(3, 2)
   real(dp) :: xshape(2, 2)
   real(dp) :: d
   real(dp) :: modulus
   real(dp) :: sign_det

   i5 = [1, 1, 2, 2, 3]
   j5 = [1, 2, 2, 3, 3]
   v5 = [4.0_dp, 1.0_dp, 5.0_dp, 2.0_dp, 3.0_dp]
   call rmumps_init_triplet(a, i5, j5, v5, 3, sym=1, info=info)
   call assert_true(info == 0, "symmetric init")
   call assert_true(rmumps_get_sym(a) == 1, "symmetric mode query")
   b3 = [6.0_dp, 17.0_dp, 13.0_dp]
   call rmumps_solve(a, b3, x3, info)
   call assert_true(info == 0, "symmetric solve status")
   call assert_close(maxval(abs(x3 - [1.0_dp, 2.0_dp, 3.0_dp])), 0.0_dp, 1.0e-12_dp, "symmetric solve")

   call rmumps_init_triplet(a, i5, j5, v5, 3, sym=2, info=info)
   call assert_true(info == 0, "symmetric positive-definite mode init")
   call assert_true(rmumps_get_sym(a) == 2, "symmetric positive-definite mode query")
   call rmumps_solve(a, b3, x3, info)
   call assert_true(info == 0, "symmetric positive-definite mode solve status")
   call assert_close(maxval(abs(x3 - [1.0_dp, 2.0_dp, 3.0_dp])), 0.0_dp, 1.0e-12_dp, &
      "symmetric positive-definite mode solve")

   i4 = [1, 1, 1, 2]
   j4 = [1, 1, 2, 2]
   v4 = [1.0_dp, 2.0_dp, 1.0_dp, 2.0_dp]
   call rmumps_init_triplet(a, i4, j4, v4, 2, info=info)
   call assert_true(info == 0, "duplicate matrix init")
   call assert_close(a%dense(1, 1), 3.0_dp, 0.0_dp, "duplicate matrix entries summed")

   i3 = [1, 1, 2]
   j3 = [1, 2, 2]
   v3 = [1.0_dp, 2.0_dp, 4.0_dp]
   call rmumps_init_triplet(a, i3, j3, v3, 2, sym=1, info=info)
   call assert_true(info == 0, "singular symmetric init")
   d = rmumps_det(a, info)
   call assert_true(info == -10, "singular factorization status")
   call assert_close(d, 0.0_dp, 0.0_dp, "singular determinant")
   call rmumps_determinant(a, modulus, sign_det, info=info)
   call assert_true(info == -10, "singular determinant decomposition status")
   call assert_close(sign_det, 0.0_dp, 0.0_dp, "singular determinant sign")

   call rmumps_init_triplet(bad, [1], [1], [1.0_dp], -1, info=info)
   call assert_true(info == -1, "negative matrix order rejected")
   call rmumps_init_triplet(bad, [1], [1, 1], [1.0_dp], 1, info=info)
   call assert_true(info == -2, "triplet length mismatch rejected")
   call rmumps_init_triplet(bad, [1], [1], [1.0_dp], 1, sym=3, info=info)
   call assert_true(info == -3, "invalid symmetry rejected")
   call rmumps_init_triplet(bad, [0], [1], [1.0_dp], 1, info=info)
   call assert_true(info == -4, "invalid matrix index rejected")

   call rmumps_init_triplet(a, [1, 2], [1, 2], [1.0_dp, 1.0_dp], 2, info=info)
   call rmumps_set_values(a, [2.0_dp], info)
   call assert_true(info == -1, "value length mismatch rejected")
   call rmumps_set_permutation(a, 1, info)
   call assert_true(info == -1, "invalid permutation code rejected")

   ri = [1]
   rj = [3]
   rv = [1.0_dp]
   call rmumps_solve(a, ri, rj, rv, 2, 2, xbad, info)
   call assert_true(info == -3, "sparse rhs output shape rejected before index check")

   rj = [3]
   call rmumps_solve(a, ri, rj, rv, 2, 2, xshape, info)
   call assert_true(info == -5, "sparse rhs column index rejected")

   print '(a)', 'test_validation: PASS'

contains

   subroutine assert_true(condition, label)
      logical, intent(in) :: condition !! Assertion condition that must evaluate true.
      character(len=*), intent(in) :: label !! Human-readable assertion label emitted on failure.

      if (.not. condition) error stop 'FAIL: ' // trim(label)
   end subroutine assert_true

   subroutine assert_close(actual, expected, tol, label)
      real(dp), intent(in) :: actual !! Computed scalar value being checked.
      real(dp), intent(in) :: expected !! Reference scalar value expected from the deterministic test.
      real(dp), intent(in) :: tol !! Maximum permitted absolute error; must be nonnegative.
      character(len=*), intent(in) :: label !! Human-readable assertion label emitted on failure.

      if (abs(actual - expected) > tol) error stop 'FAIL: ' // trim(label)
   end subroutine assert_close

end program test_validation
