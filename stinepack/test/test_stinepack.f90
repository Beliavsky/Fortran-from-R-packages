program test_stinepack
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_quiet_nan, ieee_value
   use stinepack_api, only : dp, matrix_interp_result, na_stinterp, parabola_slopes, &
                             stineman_slopes, stinterp, stinterp_result, vector_interp_result
   implicit none

   call test_parabola_slopes()
   call test_stineman_slopes()
   call test_stinterp_known_slopes()
   call test_stinterp_methods()
   call test_stinterp_bounds_and_validation()
   call test_na_stinterp_vector()
   call test_na_stinterp_matrix()

   print '(a)', 'All stinepack tests passed.'

contains

   subroutine test_parabola_slopes()
      real(dp), allocatable :: slopes(:)
      real(dp), parameter :: x(3) = [0.0_dp, 1.0_dp, 2.0_dp]
      real(dp), parameter :: y(3) = [0.0_dp, 1.0_dp, 4.0_dp]

      slopes = parabola_slopes(x, y)
      call assert_close_vector(slopes, [0.0_dp, 2.0_dp, 4.0_dp], 1.0e-14_dp, 'parabola slopes')
   end subroutine test_parabola_slopes

   subroutine test_stineman_slopes()
      real(dp), allocatable :: slopes(:)
      real(dp), parameter :: x(4) = [0.0_dp, 1.0_dp, 2.0_dp, 4.0_dp]
      real(dp), parameter :: y(4) = [0.0_dp, 2.0_dp, 1.0_dp, 3.0_dp]

      slopes = stineman_slopes(x, y, scale=.false.)
      call assert_close_vector(slopes, &
                               [4.142857142857143_dp, -0.14285714285714285_dp, &
                                -0.3333333333333333_dp, 2.3333333333333335_dp], &
                               2.0e-14_dp, 'unscaled Stineman slopes')

      slopes = stineman_slopes(x, y, scale=.true.)
      call assert_close_vector(slopes, &
                               [4.234693877551020_dp, -0.23469387755102034_dp, &
                                -0.33333333333333337_dp, 2.333333333333334_dp], &
                               2.0e-14_dp, 'scaled Stineman slopes')
   end subroutine test_stineman_slopes

   subroutine test_stinterp_known_slopes()
      type(stinterp_result) :: fit
      real(dp), parameter :: expected(5) = [0.0_dp, 0.0625_dp, 0.25_dp, 0.5625_dp, 1.0_dp]
      real(dp), parameter :: x(2) = [0.0_dp, 1.0_dp]
      real(dp), parameter :: xout(5) = [0.0_dp, 0.25_dp, 0.5_dp, 0.75_dp, 1.0_dp]
      real(dp), parameter :: y(2) = [0.0_dp, 1.0_dp]
      real(dp), parameter :: yp(2) = [0.0_dp, 2.0_dp]

      fit = stinterp(x, y, xout, yp=yp)
      call assert_status(fit%status, 0, 'stinterp known slopes status')
      call assert_close_vector(fit%x, xout, 0.0_dp, 'stinterp returns xout')
      call assert_close_vector(fit%y, expected, 2.0e-14_dp, 'stinterp known slopes')
   end subroutine test_stinterp_known_slopes

   subroutine test_stinterp_methods()
      type(stinterp_result) :: fit
      real(dp), parameter :: x(4) = [0.0_dp, 1.0_dp, 2.0_dp, 4.0_dp]
      real(dp), parameter :: xout(9) = [0.0_dp, 0.25_dp, 0.5_dp, 0.75_dp, 1.0_dp, &
                                        1.5_dp, 2.0_dp, 3.0_dp, 4.0_dp]
      real(dp), parameter :: y(4) = [0.0_dp, 2.0_dp, 1.0_dp, 3.0_dp]

      fit = stinterp(x, y, xout)
      call assert_status(fit%status, 0, 'scaledstineman status')
      call assert_close_vector(fit%y, &
                               [0.0_dp, 0.9190051020408163_dp, 1.558673469387755_dp, &
                                1.9190051020408163_dp, 2.0_dp, 1.5_dp, 1.0_dp, &
                                1.333333333333333_dp, 3.0_dp], &
                               3.0e-14_dp, 'scaledstineman interpolation')

      fit = stinterp(x, y, xout, method='st')
      call assert_status(fit%status, 0, 'stineman prefix status')
      call assert_close_vector(fit%y, &
                               [0.0_dp, 0.9017857142857143_dp, 1.5357142857142856_dp, &
                                1.9017857142857142_dp, 2.0_dp, 1.5_dp, 1.0_dp, &
                                1.3333333333333335_dp, 3.0_dp], &
                               3.0e-14_dp, 'stineman interpolation')

      fit = stinterp(x, y, xout, method='pa')
      call assert_status(fit%status, 0, 'parabola prefix status')
      call assert_close_vector(fit%y, &
                               [0.0_dp, 0.78125_dp, 1.375_dp, 1.78125_dp, 2.0_dp, &
                                1.5_dp, 1.0_dp, 1.3333333333333335_dp, 3.0_dp], &
                               3.0e-14_dp, 'parabola interpolation')
   end subroutine test_stinterp_methods

   subroutine test_stinterp_bounds_and_validation()
      type(stinterp_result) :: fit
      real(dp) :: tol
      real(dp), parameter :: x(2) = [0.0_dp, 1.0_dp]
      real(dp), parameter :: y(2) = [0.0_dp, 1.0_dp]

      tol = 5.0_dp * epsilon(1.0_dp)
      fit = stinterp(x, y, [-2.0_dp * tol, -0.5_dp * tol, 1.0_dp + 0.5_dp * tol, &
                            1.0_dp + 2.0_dp * tol])
      call assert_true(ieee_is_nan(fit%y(1)), 'far lower extrapolation is NaN')
      call assert_close_scalar(fit%y(2), -0.5_dp * tol, 1.0e-30_dp, 'tiny lower extrapolation')
      call assert_close_scalar(fit%y(3), 1.0_dp + 0.5_dp * tol, 1.0e-30_dp, 'tiny upper extrapolation')
      call assert_true(ieee_is_nan(fit%y(4)), 'far upper extrapolation is NaN')

      fit = stinterp([0.0_dp, 1.0_dp, 1.0_dp], [0.0_dp, 1.0_dp, 2.0_dp], [0.5_dp])
      call assert_true(fit%status /= 0, 'non-increasing x rejected')

      fit = stinterp([0.0_dp, 1.0_dp], y, [0.5_dp], method='s')
      call assert_true(fit%status /= 0, 'ambiguous method prefix rejected')

      fit = stinterp([0.0_dp, 1.0_dp], y, [0.5_dp], yp=[1.0_dp, 1.0_dp], method='parabola')
      call assert_true(fit%status /= 0, 'method with explicit slopes rejected')
   end subroutine test_stinterp_bounds_and_validation

   subroutine test_na_stinterp_vector()
      type(vector_interp_result) :: filled
      real(dp) :: nan
      real(dp), allocatable :: object(:)

      nan = ieee_value(0.0_dp, ieee_quiet_nan)
      object = [2.0_dp, nan, 1.0_dp, 4.0_dp, 5.0_dp, 2.0_dp]
      filled = na_stinterp(object)
      call assert_status(filled%status, 0, 'vector gap fill status')
      call assert_close_vector(filled%values, &
                               [2.0_dp, 1.227232924693520_dp, 1.0_dp, 4.0_dp, 5.0_dp, 2.0_dp], &
                               4.0e-14_dp, 'vector gap fill')

      object = [nan, 2.0_dp, nan, 1.0_dp, 4.0_dp, 5.0_dp, 2.0_dp, nan]
      filled = na_stinterp(object, na_rm=.false.)
      call assert_status(filled%status, 0, 'vector edge NaN preserve status')
      call assert_true(ieee_is_nan(filled%values(1)), 'leading NaN retained')
      call assert_true(ieee_is_nan(filled%values(8)), 'trailing NaN retained')
      call assert_close_scalar(filled%values(3), 1.227232924693520_dp, 4.0e-14_dp, 'interior NaN filled')

      filled = na_stinterp(object)
      call assert_status(filled%status, 0, 'vector edge NaN removal status')
      call assert_close_vector(filled%values, &
                               [2.0_dp, 1.227232924693520_dp, 1.0_dp, 4.0_dp, 5.0_dp, 2.0_dp], &
                               4.0e-14_dp, 'leading/trailing NaNs removed')
   end subroutine test_na_stinterp_vector

   subroutine test_na_stinterp_matrix()
      type(matrix_interp_result) :: filled
      real(dp) :: nan
      real(dp), allocatable :: object(:, :)

      nan = ieee_value(0.0_dp, ieee_quiet_nan)
      allocate (object(8, 2))
      object(:, 1) = [nan, 2.0_dp, nan, 1.0_dp, 4.0_dp, 5.0_dp, 2.0_dp, nan]
      object(:, 2) = [1.0_dp, nan, 3.0_dp, 4.0_dp, nan, 6.0_dp, 7.0_dp, 8.0_dp]

      filled = na_stinterp(object, na_rm=.false.)
      call assert_status(filled%status, 0, 'matrix gap fill status')
      call assert_true(all(shape(filled%values) == [8, 2]), 'matrix shape retained')
      call assert_close_scalar(filled%values(3, 1), 1.227232924693520_dp, 4.0e-14_dp, 'matrix first column fill')
      call assert_close_scalar(filled%values(2, 2), 2.0_dp, 4.0e-14_dp, 'matrix second column first fill')
      call assert_close_scalar(filled%values(5, 2), 5.0_dp, 4.0e-14_dp, 'matrix second column second fill')

      filled = na_stinterp(object)
      call assert_status(filled%status, 0, 'matrix row omission status')
      call assert_true(all(shape(filled%values) == [6, 2]), 'rows with unresolved edge NaNs omitted')
      call assert_close_vector(filled%values(:, 1), &
                               [2.0_dp, 1.227232924693520_dp, 1.0_dp, 4.0_dp, 5.0_dp, 2.0_dp], &
                               4.0e-14_dp, 'matrix retained first column')
      call assert_close_vector(filled%values(:, 2), &
                               [2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp, 7.0_dp], &
                               4.0e-14_dp, 'matrix retained second column')
   end subroutine test_na_stinterp_matrix

   subroutine assert_status(actual, expected, label)
      integer, intent(in) :: actual !! Actual integer status code.
      integer, intent(in) :: expected !! Expected integer status code.
      character(len=*), intent(in) :: label !! Description printed if the assertion fails.

      if (actual /= expected) then
         print '(a,2(1x,i0))', trim(label) // ':', actual, expected
         error stop 1
      end if
   end subroutine assert_status

   subroutine assert_true(condition, label)
      logical, intent(in) :: condition !! Condition that must evaluate to true.
      character(len=*), intent(in) :: label !! Description printed if the assertion fails.

      if (.not. condition) then
         print '(a)', trim(label)
         error stop 1
      end if
   end subroutine assert_true

   subroutine assert_close_scalar(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual !! Actual scalar value.
      real(dp), intent(in) :: expected !! Expected scalar reference value.
      real(dp), intent(in) :: tolerance !! Maximum allowed absolute error.
      character(len=*), intent(in) :: label !! Description printed if the assertion fails.

      if (abs(actual - expected) > tolerance) then
         print '(a,3(1x,es24.16))', trim(label) // ':', actual, expected, abs(actual - expected)
         error stop 1
      end if
   end subroutine assert_close_scalar

   subroutine assert_close_vector(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual(:) !! Actual vector values.
      real(dp), intent(in) :: expected(:) !! Expected vector reference values.
      real(dp), intent(in) :: tolerance !! Maximum allowed elementwise absolute error.
      character(len=*), intent(in) :: label !! Description printed if the assertion fails.

      if (size(actual) /= size(expected)) then
         print '(a,2(1x,i0))', trim(label) // ' size:', size(actual), size(expected)
         error stop 1
      end if
      if (size(actual) > 0) then
         if (maxval(abs(actual - expected)) > tolerance) then
            print '(a,1x,es24.16)', trim(label) // ' max error:', maxval(abs(actual - expected))
            error stop 1
         end if
      end if
   end subroutine assert_close_vector

end program test_stinepack
