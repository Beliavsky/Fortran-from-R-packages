program test_class_state
   use rmumps_api, only : dp, rmumps_matrix, rmumps_init_triplet, rmumps_init_triplet_zero_based, &
      rmumps_solve, rmumps_numeric, &
      rmumps_set_copy, rmumps_get_copy, rmumps_set_ncore, rmumps_get_ncore, rmumps_set_rhs, &
      rmumps_get_rhs, rmumps_set_mrhs, rmumps_get_mrhs, rmumps_set_perm_in, rmumps_get_icntl
   implicit none

   type(rmumps_matrix) :: a
   type(rmumps_matrix) :: a0
   integer :: info
   integer :: icntl(33)
   real(dp), allocatable :: stored_rhs(:)
   real(dp), allocatable :: stored_mrhs(:, :)
   real(dp) :: b(2)
   real(dp) :: x(2)
   real(dp) :: bmat(2, 2)
   real(dp) :: xmat(2, 2)

   call rmumps_init_triplet(a, [1, 1, 2, 2], [1, 2, 1, 2], &
      [3.0_dp, 1.0_dp, 1.0_dp, 2.0_dp], 2, info=info)
   call assert_true(info == 0, "class state init")
   call assert_true(rmumps_get_copy(a), "default copy flag")

   call rmumps_init_triplet_zero_based(a0, [0, 0, 1, 1], [0, 1, 0, 1], &
      [3.0_dp, 1.0_dp, 1.0_dp, 2.0_dp], 2, copy=.false., info=info)
   call assert_true(info == 0, "zero-based constructor")
   call assert_true(.not. rmumps_get_copy(a0), "constructor copy flag")
   call assert_matrix_close(a0%dense, a%dense, 0.0_dp, "zero-based constructor matrix")
   call assert_true(rmumps_get_ncore(a) == 1, "default ncore")

   call rmumps_set_copy(a, .false.)
   call assert_true(.not. rmumps_get_copy(a), "copy flag state")
   call rmumps_set_ncore(a, 4, info)
   call assert_true(info == 0 .and. rmumps_get_ncore(a) == 4, "ncore state")
   call rmumps_set_ncore(a, 0, info)
   call assert_true(info == -1 .and. rmumps_get_ncore(a) == 4, "invalid ncore")

   b = [5.0_dp, 5.0_dp]
   call rmumps_set_rhs(a, b, info)
   call assert_true(info == 0, "set rhs")
   call rmumps_get_rhs(a, stored_rhs)
   call assert_vector_close(stored_rhs, b, 0.0_dp, "stored rhs before solve")

   call rmumps_solve(a, b, x, info)
   call assert_true(info == 0, "solve for rhs state")
   call rmumps_get_rhs(a, stored_rhs)
   call assert_vector_close(stored_rhs, x, 0.0_dp, "stored rhs solution")

   bmat(:, 1) = b
   bmat(:, 2) = [2.0_dp, 3.0_dp]
   call rmumps_set_mrhs(a, bmat, info)
   call assert_true(info == 0, "set mrhs")
   call rmumps_get_mrhs(a, stored_mrhs)
   call assert_matrix_close(stored_mrhs, bmat, 0.0_dp, "stored mrhs before solve")

   call rmumps_solve(a, bmat, xmat, info)
   call assert_true(info == 0, "solve for mrhs state")
   call rmumps_get_mrhs(a, stored_mrhs)
   call assert_matrix_close(stored_mrhs, xmat, 0.0_dp, "stored mrhs solution")

   call rmumps_numeric(a, info)
   call assert_true(info == 0 .and. a%numeric_ready, "numeric before user permutation")
   call rmumps_set_perm_in(a, [2, 1], info)
   call assert_true(info == 0, "set user permutation")
   call assert_true(allocated(a%perm_in), "user permutation stored")
   call assert_true(all(a%perm_in == [2, 1]), "user permutation content")
   call assert_true(.not. a%symbolic_ready .and. .not. a%numeric_ready, "user permutation invalidates factors")
   icntl = rmumps_get_icntl(a)
   call assert_true(icntl(7) == 1 .and. icntl(28) == 1, "user permutation controls")

   call rmumps_set_perm_in(a, [1, 1], info)
   call assert_true(info == -2, "duplicate user permutation rejected")
   call rmumps_set_perm_in(a, [1], info)
   call assert_true(info == -1, "short user permutation rejected")

   print '(a)', 'test_class_state: PASS'

contains

   subroutine assert_true(condition, label)
      logical, intent(in) :: condition !! Assertion condition that must evaluate true.
      character(len=*), intent(in) :: label !! Human-readable assertion label emitted on failure.

      if (.not. condition) error stop 'FAIL: ' // trim(label)
   end subroutine assert_true

   subroutine assert_vector_close(actual, expected, tol, label)
      real(dp), intent(in) :: actual(:) !! Computed vector being checked elementwise.
      real(dp), intent(in) :: expected(:) !! Reference vector with the same shape as actual.
      real(dp), intent(in) :: tol !! Maximum permitted elementwise absolute error; must be nonnegative.
      character(len=*), intent(in) :: label !! Human-readable assertion label emitted on failure.

      if (size(actual) /= size(expected)) error stop 'FAIL size: ' // trim(label)
      if (size(actual) > 0) then
         if (maxval(abs(actual - expected)) > tol) error stop 'FAIL: ' // trim(label)
      end if
   end subroutine assert_vector_close

   subroutine assert_matrix_close(actual, expected, tol, label)
      real(dp), intent(in) :: actual(:, :) !! Computed matrix being checked elementwise.
      real(dp), intent(in) :: expected(:, :) !! Reference matrix with the same shape as actual.
      real(dp), intent(in) :: tol !! Maximum permitted elementwise absolute error; must be nonnegative.
      character(len=*), intent(in) :: label !! Human-readable assertion label emitted on failure.

      if (any(shape(actual) /= shape(expected))) error stop 'FAIL shape: ' // trim(label)
      if (size(actual) > 0) then
         if (maxval(abs(actual - expected)) > tol) error stop 'FAIL: ' // trim(label)
      end if
   end subroutine assert_matrix_close

end program test_class_state
