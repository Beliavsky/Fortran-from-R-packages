program test_rhs_modes
   use rmumps_api, only : dp, rmumps_matrix, rmumps_init_triplet, rmumps_solve, rmumps_solve_inplace, &
      rmumps_solvet, rmumps_inverse, rmumps_determinant, rmumps_set_icntl, rmumps_get_icntl
   implicit none

   type(rmumps_matrix) :: a
   integer, parameter :: n = 3
   integer :: info
   integer :: k
   integer :: icntl(33)
   integer :: ai(6)
   integer :: aj(6)
   integer :: ri(7)
   integer :: rj(7)
   integer :: rti(6)
   integer :: rtj(6)
   integer :: iz(0)
   real(dp) :: av(6)
   real(dp) :: rv(7)
   real(dp) :: rtv(6)
   real(dp) :: vz(0)
   real(dp) :: xtrue(n, 2)
   real(dp) :: b(n, 2)
   real(dp) :: bt(n, 2)
   real(dp) :: x(n, 2)
   real(dp) :: xinplace(n, 2)
   real(dp) :: xsparse(n, 2)
   real(dp) :: xt(n, 2)
   real(dp) :: xtsparse(n, 2)
   real(dp) :: ainv(n, n)
   real(dp) :: ainv_sparse(n, n)
   real(dp) :: ainvt_sparse(n, n)
   real(dp) :: eye(n, n)
   real(dp) :: expected_eye(n, n)
   real(dp) :: modulus
   real(dp) :: sign_det

   ai = [1, 1, 2, 2, 3, 3]
   aj = [2, 3, 1, 2, 1, 3]
   av = [2.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 2.0_dp, 1.0_dp]
   call rmumps_init_triplet(a, ai, aj, av, n, info=info)
   call assert_true(info == 0, "pivoting matrix init")

   xtrue(:, 1) = [1.0_dp, 2.0_dp, 3.0_dp]
   xtrue(:, 2) = [-1.0_dp, 0.5_dp, 4.0_dp]
   b = matmul(a%dense, xtrue)
   bt = matmul(transpose(a%dense), xtrue)

   call rmumps_solve(a, b, x, info)
   call assert_true(info == 0, "generic dense matrix solve status")
   call assert_matrix_close(x, xtrue, 1.0e-12_dp, "generic dense matrix solve")

   xinplace = b
   call rmumps_solve_inplace(a, xinplace, info)
   call assert_true(info == 0, "generic dense matrix inplace status")
   call assert_matrix_close(xinplace, xtrue, 1.0e-12_dp, "generic dense matrix inplace")

   ri = [1, 1, 2, 3, 1, 2, 3]
   rj = [1, 1, 1, 1, 2, 2, 2]
   rv = [3.0_dp, 4.0_dp, 3.0_dp, 5.0_dp, 5.0_dp, -0.5_dp, 2.0_dp]
   call rmumps_solve(a, ri, rj, rv, n, 2, xsparse, info)
   call assert_true(info == 0, "sparse rhs solve status")
   call assert_matrix_close(xsparse, xtrue, 1.0e-12_dp, "sparse rhs solve with duplicate")
   icntl = rmumps_get_icntl(a)
   call assert_true(icntl(20) == 1, "sparse rhs leaves ICNTL(20) sparse mode")

   expected_eye = 0.0_dp
   do k = 1, n
      expected_eye(k, k) = 1.0_dp
   end do
   call rmumps_solve(a, ai, aj, av, n, n, eye, info)
   call assert_true(info == 0, "matrix-as-sparse-rhs solve status")
   call assert_matrix_close(eye, expected_eye, 1.0e-12_dp, "matrix-as-sparse-rhs identity")

   call rmumps_set_icntl(a, [2], [9], info)
   call rmumps_solve(a, bt, xt, info)
   call assert_true(info == 0, "ICNTL(9) transpose solve status")
   call assert_matrix_close(xt, xtrue, 1.0e-12_dp, "ICNTL(9) transpose solve")
   icntl = rmumps_get_icntl(a)
   call assert_true(icntl(9) == 2, "ordinary solve preserves ICNTL(9)")

   call rmumps_solvet(a, bt, xt, info)
   call assert_true(info == 0, "generic transpose matrix solve status")
   call assert_matrix_close(xt, xtrue, 1.0e-12_dp, "generic transpose matrix solve")
   icntl = rmumps_get_icntl(a)
   call assert_true(icntl(9) == 1, "solvet restores ICNTL(9) to normal solve mode")

   rti = [1, 2, 3, 1, 2, 3]
   rtj = [1, 1, 1, 2, 2, 2]
   rtv = [8.0_dp, 4.0_dp, 4.0_dp, 8.5_dp, -1.5_dp, 3.0_dp]
   call rmumps_solvet(a, rti, rtj, rtv, n, 2, xtsparse, info)
   call assert_true(info == 0, "sparse rhs transpose solve status")
   call assert_matrix_close(xtsparse, xtrue, 1.0e-12_dp, "sparse rhs transpose solve")
   icntl = rmumps_get_icntl(a)
   call assert_true(icntl(9) == 1 .and. icntl(20) == 1, "sparse transpose control state")

   call rmumps_inverse(a, ainv, info)
   call assert_true(info == 0, "inverse status")
   call rmumps_solve(a, iz, iz, vz, 0, 0, ainv_sparse, info)
   call assert_true(info == 0, "0x0 sparse rhs inverse status")
   call assert_matrix_close(ainv_sparse, ainv, 1.0e-12_dp, "0x0 sparse rhs inverse")

   call rmumps_solvet(a, iz, iz, vz, 0, 0, ainvt_sparse, info)
   call assert_true(info == 0, "0x0 sparse rhs inverse transpose status")
   call assert_matrix_close(ainvt_sparse, transpose(ainv), 1.0e-12_dp, "0x0 sparse rhs inverse transpose")

   call rmumps_determinant(a, modulus, sign_det, info=info)
   call assert_true(info == 0, "determinant log status")
   call assert_close(modulus, log(4.0_dp), 1.0e-12_dp, "determinant log modulus")
   call assert_close(sign_det, -1.0_dp, 0.0_dp, "determinant sign")

   call rmumps_determinant(a, modulus, sign_det, logarithm=.false., info=info)
   call assert_true(info == 0, "determinant absolute status")
   call assert_close(modulus, 4.0_dp, 1.0e-12_dp, "determinant absolute modulus")
   call assert_close(sign_det, -1.0_dp, 0.0_dp, "determinant absolute sign")

   print '(a)', 'test_rhs_modes: PASS'

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

   subroutine assert_matrix_close(actual, expected, tol, label)
      real(dp), intent(in) :: actual(:, :) !! Computed matrix being checked elementwise.
      real(dp), intent(in) :: expected(:, :) !! Reference matrix with the same shape as actual.
      real(dp), intent(in) :: tol !! Maximum permitted elementwise absolute error; must be nonnegative.
      character(len=*), intent(in) :: label !! Human-readable assertion label emitted on failure.

      if (any(shape(actual) /= shape(expected))) error stop 'FAIL shape: ' // trim(label)
      if (maxval(abs(actual - expected)) > tol) error stop 'FAIL: ' // trim(label)
   end subroutine assert_matrix_close

end program test_rhs_modes
