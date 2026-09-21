program test_core
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_quiet_nan, ieee_value
   use funData_api
   implicit none

   type(fun_data) :: f1
   type(fun_data) :: f2
   type(fun_data) :: mean1
   type(fun_data) :: interp
   type(fun_data) :: tprod
   type(fun_data) :: flipped
   type(fun_data) :: reg
   type(fun_data) :: f2d
   type(fun_data) :: f3d
   type(irreg_fun_data) :: i1
   type(irreg_fun_data) :: back
   type(irreg_fun_data) :: mixed_flipped
   type(irreg_fun_data) :: mixed_new
   type(multi_fun_data) :: m1
   type(real_vector) :: ia(3)
   type(real_vector) :: ix(3)
   type(real_vector) :: grids2(2)
   type(real_vector) :: grids3(3)
   real(dp), allocatable :: values(:)
   real(dp), allocatable :: values2(:)
   real(dp), allocatable :: weights(:)
   real(dp) :: x1(4, 5)
   real(dp) :: x2(4, 6)
   real(dp) :: qnan
   integer :: i
   integer :: j
   logical :: ok

   qnan = ieee_value(0.0_dp, ieee_quiet_nan)
   do j = 1, 5
      do i = 1, 4
         x1(i, j) = real(i + (j - 1) * 4, dp)
      end do
   end do
   do j = 1, 6
      do i = 1, 4
         x2(i, j) = real(i + (j - 1) * 4, dp)
      end do
   end do
   call create_fun_data_1d([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp], x1, f1, ok)
   call assert_true(ok, "create f1")
   call create_fun_data_1d([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp], x2, f2, ok)
   call assert_true(ok, "create f2")

   call int_weights(f1%argvals(1)%values, "trapezoidal", weights, ok)
   call assert_true(ok, "weights")
   call assert_close_vec(weights, [0.5_dp, 1.0_dp, 1.0_dp, 1.0_dp, 0.5_dp], 1.0e-12_dp, "weights values")

   call integrate_fun_data(f1, values, ok = ok)
   call assert_true(ok, "integrate f1")
   call assert_close_vec(values, [36.0_dp, 40.0_dp, 44.0_dp, 48.0_dp], 1.0e-12_dp, "integrate f1 values")

   grids2(1)%values = [0.0_dp, 0.5_dp, 1.0_dp]
   grids2(2)%values = [0.0_dp, 1.0_dp, 2.0_dp]
   call create_fun_data(grids2, [1, 3, 3], [(2.0_dp, i = 1, 9)], f2d, ok)
   call assert_true(ok, "create 2D constant")
   call integrate_fun_data(f2d, values2, ok = ok)
   call assert_true(ok, "integrate 2D")
   call assert_close(values2(1), 4.0_dp, 1.0e-12_dp, "2D constant integral")

   grids3(1)%values = [0.0_dp, 0.5_dp, 1.0_dp]
   grids3(2)%values = [0.0_dp, 1.0_dp, 2.0_dp]
   grids3(3)%values = [0.0_dp, 1.5_dp, 3.0_dp]
   call create_fun_data(grids3, [1, 3, 3, 3], [(1.0_dp, i = 1, 27)], f3d, ok)
   call assert_true(ok, "create 3D constant")
   call integrate_fun_data(f3d, values2, ok = ok)
   call assert_true(ok, "integrate 3D")
   call assert_close(values2(1), 6.0_dp, 1.0e-12_dp, "3D constant integral")

   call norm_fun_data(f1, values, ok = ok)
   call assert_true(ok, "norm f1")
   call assert_close(values(1), 420.0_dp, 1.0e-12_dp, "squared norm f1 first")
   call scalar_product_fun_data(f1, f1, values2, ok = ok)
   call assert_true(ok, "scalar product f1")
   call assert_close_vec(values2, values, 1.0e-12_dp, "scalar product equals squared norm")

   call mean_fun_data(f1, mean1, ok = ok)
   call assert_true(ok, "mean f1")
   call assert_close_vec(mean1%x, [2.5_dp, 6.5_dp, 10.5_dp, 14.5_dp, 18.5_dp], 1.0e-12_dp, "mean f1 values")

   ia(1)%values = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
   ix(1)%values = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
   ia(2)%values = [2.0_dp, 3.0_dp, 4.0_dp]
   ix(2)%values = [2.0_dp, 3.0_dp, 4.0_dp]
   ia(3)%values = [3.0_dp, 4.0_dp, 5.0_dp]
   ix(3)%values = [-3.0_dp, -2.0_dp, -1.0_dp]
   call create_irreg_fun_data(ia, ix, i1, ok)
   call assert_true(ok, "create irregular")
   call integrate_irreg_fun_data(i1, values, ok = ok)
   call assert_true(ok, "integrate irregular")
   call assert_close_vec(values, [12.0_dp, 6.0_dp, -4.0_dp], 1.0e-12_dp, "irregular integrals")
   call integrate_irreg_fun_data(i1, values, full_domain = .true., ok = ok)
   call assert_true(ok, "integrate irregular full domain")
   call assert_close_vec(values, [12.0_dp, 12.0_dp, -12.0_dp], 1.0e-12_dp, "irregular full-domain integrals")
   call norm_irreg_fun_data(i1, values, ok = ok)
   call assert_true(ok, "norm irregular")
   call assert_close_vec(values, [42.0_dp, 19.0_dp, 9.0_dp], 1.0e-12_dp, "irregular squared norms")
   call norm_irreg_fun_data(i1, values, full_domain = .true., ok = ok)
   call assert_true(ok, "norm irregular full")
   call assert_close_vec(values, [42.0_dp, 42.0_dp, 43.0_dp], 1.0e-12_dp, "irregular full squared norms")

   call as_fun_data(i1, reg, ok)
   call assert_true(ok, "irregular to regular")
   call assert_true(all(reg%dims == [3, 5]), "irregular union grid dimensions")
   call as_irreg_fun_data(reg, back, ok)
   call assert_true(ok, "regular to irregular")
   call assert_true(size(back%curves) == 3, "regular to irregular observation count")
   call assert_close_vec(back%curves(3)%x, [-3.0_dp, -2.0_dp, -1.0_dp], 1.0e-12_dp, "round-trip irregular")

   allocate(mixed_new%curves(1))
   mixed_new%curves(1) = i1%curves(1)
   mixed_new%curves(1)%x = -mixed_new%curves(1)%x
   call extract_fun_data(f1, [1], result = reg, ok = ok)
   call assert_true(ok, "regular reference extraction")
   call flip_fun_irreg_data(reg, mixed_new, mixed_flipped, ok)
   call assert_true(ok, "mixed regular-irregular flip")
   call assert_close_vec(mixed_flipped%curves(1)%x, i1%curves(1)%x, 1.0e-12_dp, "mixed flipped orientation")

   interp = f1
   interp%x(1 + (2 - 1) * 4) = qnan
   interp%x(1 + (3 - 1) * 4) = qnan
   call approx_na(interp, reg, ok)
   call assert_true(ok, "approx NA")
   call assert_close(fun_value(reg, 1, [2]), 5.0_dp, 1.0e-12_dp, "interpolated first gap")
   call assert_close(fun_value(reg, 1, [3]), 9.0_dp, 1.0e-12_dp, "interpolated second gap")

   call tensor_product2(f1, mean1, tprod, ok)
   call assert_true(ok, "tensor product")
   call assert_true(all(tprod%dims == [4, 5, 5]), "tensor dimensions")
   call assert_close(fun_value(tprod, 1, [2, 3]), 5.0_dp * 10.5_dp, 1.0e-12_dp, "tensor value")

   flipped = f1
   flipped%x = -flipped%x
   call flip_fun_data(f1, flipped, reg, ok)
   call assert_true(ok, "flip")
   call assert_close_vec(reg%x, f1%x, 1.0e-12_dp, "flipped orientation")

   allocate(m1%components(2))
   m1%components(1) = f1
   m1%components(2) = f2
   call integrate_multi_fun_data(m1, values, ok = ok)
   call assert_true(ok, "integrate multi")
   call integrate_fun_data(f1, values2, ok = ok)
   call integrate_fun_data(f2, weights, ok = ok)
   call assert_close_vec(values, values2 + weights, 1.0e-12_dp, "multi integral sum")

   print '(a)', "All funData core tests passed"

contains

   subroutine assert_true(condition, message)
      logical, intent(in) :: condition !! Condition that must evaluate true for the test to continue.
      character(len=*), intent(in) :: message !! Human-readable label printed when the assertion fails.
      if (.not. condition) then
         print '(a)', "FAIL: " // trim(message)
         error stop 1
      end if
   end subroutine assert_true

   subroutine assert_close(actual, expected, tol, message)
      real(dp), intent(in) :: actual !! Computed scalar result being checked.
      real(dp), intent(in) :: expected !! Reference scalar value.
      real(dp), intent(in) :: tol !! Maximum accepted absolute difference.
      character(len=*), intent(in) :: message !! Human-readable assertion label.
      if (abs(actual - expected) > tol) then
         print '(a,2es24.15)', "FAIL: " // trim(message) // " actual/expected ", actual, expected
         error stop 1
      end if
   end subroutine assert_close

   subroutine assert_close_vec(actual, expected, tol, message)
      real(dp), intent(in) :: actual(:) !! Computed vector result being checked.
      real(dp), intent(in) :: expected(:) !! Reference vector with the same shape as actual.
      real(dp), intent(in) :: tol !! Maximum accepted absolute elementwise difference.
      character(len=*), intent(in) :: message !! Human-readable assertion label.
      if (size(actual) /= size(expected)) then
         print '(a)', "FAIL size: " // trim(message)
         error stop 1
      end if
      if (size(actual) > 0) then
         if (maxval(abs(actual - expected)) > tol) then
            print '(a)', "FAIL values: " // trim(message)
            error stop 1
         end if
      end if
   end subroutine assert_close_vec

end program test_core
