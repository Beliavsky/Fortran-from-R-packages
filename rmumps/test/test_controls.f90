program test_controls
   use rmumps_api, only : dp, rmumps_matrix, rmumps_infos, rmumps_init_triplet, rmumps_numeric, &
      rmumps_set_icntl, rmumps_get_icntl, rmumps_set_cntl, rmumps_get_cntl, rmumps_set_keep, &
      rmumps_get_keep, rmumps_get_infos, rmumps_get_permutation, rmumps_set_permutation, RMUMPS_PERM_QAMD
   implicit none

   type(rmumps_matrix) :: a
   type(rmumps_infos) :: infos
   integer :: info
   integer :: icntl(33)
   integer :: keep(500)
   real(dp) :: cntl(5)
   real(dp) :: det_from_info

   call rmumps_init_triplet(a, [1, 1, 2, 2, 3, 3], [2, 3, 1, 2, 1, 3], &
      [2.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 2.0_dp, 1.0_dp], 3, info=info)
   call assert_true(info == 0, "control test init")

   icntl = rmumps_get_icntl(a)
   call assert_true(icntl(1) == -1, "default ICNTL(1)")
   call assert_true(icntl(2) == -1, "default ICNTL(2)")
   call assert_true(icntl(3) == -1, "default ICNTL(3)")
   call assert_true(icntl(7) == 7, "default ICNTL(7)")
   call assert_true(icntl(8) == 77, "default ICNTL(8)")
   call assert_true(icntl(14) == 50, "default ICNTL(14)")
   call assert_true(icntl(27) == -8, "default ICNTL(27)")
   call assert_true(icntl(33) == 1, "default ICNTL(33)")

   call rmumps_numeric(a, info)
   call assert_true(info == 0 .and. a%numeric_ready, "numeric before control invalidation")
   call rmumps_set_icntl(a, [0, 999], [7, 34], info)
   call assert_true(info == 0, "set ICNTL status")
   call assert_true(rmumps_get_permutation(a) == 0, "ICNTL(7) synchronizes permutation")
   call assert_true(a%numeric_ready .and. a%symbolic_ready, "direct ICNTL(7) preserves cached factors")
   icntl = rmumps_get_icntl(a)
   call assert_true(icntl(7) == 0, "ICNTL(7) update")

   call rmumps_set_permutation(a, RMUMPS_PERM_QAMD, info)
   call assert_true(info == 0, "set permutation status")
   icntl = rmumps_get_icntl(a)
   call assert_true(icntl(7) == RMUMPS_PERM_QAMD, "set permutation synchronizes ICNTL(7)")

   call rmumps_set_cntl(a, [0.25_dp, 0.5_dp, 99.0_dp], [1, 5, 6], info)
   call assert_true(info == 0, "set CNTL status")
   cntl = rmumps_get_cntl(a)
   call assert_close(cntl(1), 0.25_dp, 0.0_dp, "CNTL(1)")
   call assert_close(cntl(5), 0.5_dp, 0.0_dp, "CNTL(5)")

   call rmumps_set_keep(a, [17, 29, 101], [1, 500, 501], info)
   call assert_true(info == 0, "set KEEP status")
   keep = rmumps_get_keep(a)
   call assert_true(keep(1) == 17, "KEEP(1)")
   call assert_true(keep(500) == 29, "KEEP(500)")

   call rmumps_init_triplet(a, [1], [1], [1.0_dp], 1, info=info)
   keep = rmumps_get_keep(a)
   call assert_true(keep(34) == storage_size(0) / 8, "default KEEP(34) integer byte size")
   call assert_true(keep(10) == 8 / keep(34), "default KEEP(10) real-to-integer size ratio")

   call rmumps_set_icntl(a, [1, 2], [1], info)
   call assert_true(info == -1, "ICNTL length mismatch")
   call rmumps_set_cntl(a, [1.0_dp], [1, 2], info)
   call assert_true(info == -1, "CNTL length mismatch")
   call rmumps_set_keep(a, [1, 2], [1], info)
   call assert_true(info == -1, "KEEP length mismatch")

   call rmumps_init_triplet(a, [1, 1, 2, 2, 3, 3], [2, 3, 1, 2, 1, 3], &
      [2.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 2.0_dp, 1.0_dp], 3, info=info)
   call rmumps_numeric(a, info)
   call assert_true(info == 0, "numeric before get_infos")
   call rmumps_get_infos(a, infos)
   call assert_true(infos%info(1) == 0, "get_infos local status")
   call assert_true(infos%infog(1) == 0, "get_infos global status")
   det_from_info = infos%rinfog(12) * 2.0_dp ** infos%infog(34)
   call assert_close(det_from_info, -4.0_dp, 1.0e-12_dp, "get_infos determinant encoding")

   print '(a)', 'test_controls: PASS'

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

end program test_controls
