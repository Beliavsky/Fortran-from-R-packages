program test_rmumps
   use rmumps_api, only : dp, rmumps_matrix, rmumps_init_triplet, rmumps_solve, rmumps_solvet, &
      rmumps_inverse, rmumps_det, rmumps_log2det, rmumps_set_values, rmumps_get_triplet, &
      rmumps_set_permutation, rmumps_get_permutation, rmumps_dim, rmumps_nrow, rmumps_ncol, &
      RMUMPS_PERM_AMD
   implicit none

   type(rmumps_matrix) :: a
   integer, parameter :: n = 3
   integer :: info
   integer :: ii(7)
   integer :: jj(7)
   integer :: dims(2)
   integer, allocatable :: io(:)
   integer, allocatable :: jo(:)
   real(dp) :: vv(7)
   real(dp) :: b(n)
   real(dp) :: bt(n)
   real(dp) :: x(n)
   real(dp) :: xt(n)
   real(dp) :: ainv(n, n)
   real(dp) :: ident(n, n)
   real(dp) :: d
   real(dp) :: logd
   real(dp) :: sign_det
   real(dp), allocatable :: vo(:)

   ii = [1, 1, 2, 2, 2, 3, 3]
   jj = [1, 2, 1, 2, 3, 2, 3]
   vv = [4.0_dp, 1.0_dp, 2.0_dp, 5.0_dp, 1.0_dp, 2.0_dp, 3.0_dp]
   call rmumps_init_triplet(a, ii, jj, vv, n, info=info)
   call assert_true(info == 0, "init")

   b = matmul(a%dense, [1.0_dp, 2.0_dp, 3.0_dp])
   call rmumps_solve(a, b, x, info)
   call assert_true(info == 0, "solve status")
   call assert_close(maxval(abs(x - [1.0_dp, 2.0_dp, 3.0_dp])), 0.0_dp, 1.0e-12_dp, "solve")

   bt = matmul(transpose(a%dense), [1.0_dp, 2.0_dp, 3.0_dp])
   call rmumps_solvet(a, bt, xt, info)
   call assert_true(info == 0, "solvet status")
   call assert_close(maxval(abs(xt - [1.0_dp, 2.0_dp, 3.0_dp])), 0.0_dp, 1.0e-12_dp, "solvet")

   call rmumps_inverse(a, ainv, info)
   call assert_true(info == 0, "inverse status")
   ident = matmul(a%dense, ainv)
   call assert_close(maxval(abs(ident - eye(n))), 0.0_dp, 1.0e-12_dp, "inverse")

   d = rmumps_det(a, info)
   call assert_true(info == 0, "det status")
   call assert_close(d, 46.0_dp, 1.0e-12_dp, "det")
   logd = rmumps_log2det(a, sign_det, info)
   call assert_true(info == 0, "log2det status")
   call assert_close(logd, log(46.0_dp) / log(2.0_dp), 1.0e-12_dp, "log2det")
   call assert_close(sign_det, 1.0_dp, 0.0_dp, "det sign")

   call rmumps_get_triplet(a, io, jo, vo)
   call assert_true(all(io == ii), "triplet rows")
   call assert_true(all(jo == jj), "triplet cols")
   call assert_close(maxval(abs(vo - vv)), 0.0_dp, 0.0_dp, "triplet values")

   call rmumps_set_permutation(a, RMUMPS_PERM_AMD, info)
   call assert_true(info == 0, "set permutation")
   call assert_true(rmumps_get_permutation(a) == RMUMPS_PERM_AMD, "get permutation")

   dims = rmumps_dim(a)
   call assert_true(all(dims == [n, n]), "dim")
   call assert_true(rmumps_nrow(a) == n, "nrow")
   call assert_true(rmumps_ncol(a) == n, "ncol")

   vv = 2.0_dp * vv
   call rmumps_set_values(a, vv, info)
   call assert_true(info == 0, "set values")
   b = matmul(a%dense, [1.0_dp, 2.0_dp, 3.0_dp])
   call rmumps_solve(a, b, x, info)
   call assert_close(maxval(abs(x - [1.0_dp, 2.0_dp, 3.0_dp])), 0.0_dp, 1.0e-12_dp, "solve after update")

   print '(a)', 'test_rmumps: PASS'

contains

   pure function eye(m) result(aeye)
      integer, intent(in) :: m !! Order of the identity matrix to construct.
      real(dp) :: aeye(m, m)
      integer :: k

      aeye = 0.0_dp
      do k = 1, m
         aeye(k, k) = 1.0_dp
      end do
   end function eye

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

end program test_rmumps
