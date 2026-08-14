program test_splines
   use splines, only : dp, b_spline_t, poly_spline_t, spline_design, &
      linear_interp, fit_interpolating_spline, fit_periodic_spline, &
      to_polynomial_spline, inverse_monotone_spline, bs_basis, &
      natural_spline_basis, type7_quantile
   implicit none

   integer :: failures
   failures = 0
   call test_design(failures)
   call test_quantile(failures)
   call test_bs_reference(failures)
   call test_linear(failures)
   call test_interpolation(failures)
   call test_periodic(failures)
   call test_regression_bases(failures)
   call test_inverse(failures)

   if (failures /= 0) then
      write(*, '(a,i0)') 'FAILED tests: ', failures
      error stop 1
   end if
   print '(a)', 'All splines tests passed.'

contains

   subroutine check(condition, message, failures)
      logical, intent(in) :: condition
      character(*), intent(in) :: message
      integer, intent(inout) :: failures
      if (.not. condition) then
         failures = failures + 1
         write(*, '(a)') 'FAIL: ' // message
      end if
   end subroutine check

   subroutine check_close(actual, expected, tolerance, message, failures)
      real(dp), intent(in) :: actual, expected, tolerance
      character(*), intent(in) :: message
      integer, intent(inout) :: failures
      call check(abs(actual - expected) <= tolerance, message, failures)
      if (abs(actual - expected) > tolerance) then
         write(*, '(a,es24.16,a,es24.16)') '  actual=', actual, ' expected=', expected
      end if
   end subroutine check_close

   subroutine test_design(failures)
      integer, intent(inout) :: failures
      real(dp), parameter :: knots(10) = [0.0_dp,0.0_dp,0.0_dp,0.0_dp,1.0_dp, &
         2.0_dp,3.0_dp,3.0_dp,3.0_dp,3.0_dp]
      real(dp), parameter :: x(7) = [0.0_dp,0.25_dp,0.5_dp,1.0_dp,1.5_dp,2.5_dp,3.0_dp]
      real(dp), allocatable :: d(:, :)
      integer, allocatable :: derivs(:)
      integer :: status

      call spline_design(knots, x, 4, d, status=status)
      call check(status == 0, 'spline_design status', failures)
      call check_close(d(2,1), 0.421875_dp, 1.0e-13_dp, 'basis reference 1', failures)
      call check_close(d(2,2), 0.49609375_dp, 1.0e-13_dp, 'basis reference 2', failures)
      call check_close(d(4,3), 7.0_dp/12.0_dp, 1.0e-13_dp, 'basis reference 3', failures)
      call check_close(d(7,6), 1.0_dp, 1.0e-13_dp, 'right endpoint basis', failures)
      call check(maxval(abs(sum(d, dim=2) - 1.0_dp)) < 1.0e-13_dp, 'partition of unity', failures)

      allocate(derivs(size(x))); derivs = 1
      call spline_design(knots, x, 4, d, derivs, status)
      call check(status == 0, 'derivative design status', failures)
      call check_close(d(1,1), -3.0_dp, 1.0e-12_dp, 'first derivative left', failures)
      call check_close(d(5,3), -0.5625_dp, 1.0e-12_dp, 'first derivative middle', failures)
      call check_close(d(7,6), 3.0_dp, 1.0e-12_dp, 'first derivative right', failures)
      derivs = 2
      call spline_design(knots, x, 4, d, derivs, status)
      call check_close(d(1,2), -9.0_dp, 1.0e-11_dp, 'second derivative reference', failures)
      call check_close(d(5,5), 0.75_dp, 1.0e-11_dp, 'second derivative middle', failures)
   end subroutine test_design

   subroutine test_quantile(failures)
      integer, intent(inout) :: failures
      real(dp) :: q
      q = type7_quantile([1.0_dp,4.0_dp,2.0_dp,3.0_dp], 0.25_dp)
      call check_close(q, 1.75_dp, 1.0e-14_dp, 'R type-7 quantile', failures)
   end subroutine test_quantile


   subroutine test_bs_reference(failures)
      integer, intent(inout) :: failures
      real(dp), parameter :: x(4) = [-0.5_dp,0.25_dp,1.5_dp,3.5_dp]
      real(dp), allocatable :: b(:, :)
      integer :: status
      call bs_basis(x, b, degree=3, knots=[1.0_dp,2.0_dp], intercept=.true., &
         boundary_knots=[0.0_dp,3.0_dp], status=status)
      call check(status == 0, 'explicit B-spline basis status', failures)
      call check_close(b(1,1), 3.375_dp, 1.0e-13_dp, 'B-spline left extrapolation', failures)
      call check_close(b(1,2), -2.84375_dp, 1.0e-13_dp, 'B-spline left reference', failures)
      call check_close(b(2,3), 0.079427083333333333_dp, 1.0e-13_dp, &
         'B-spline interior reference', failures)
      call check_close(b(3,4), 0.46875_dp, 1.0e-13_dp, 'B-spline center reference', failures)
      call check_close(b(4,6), 3.375_dp, 1.0e-13_dp, 'B-spline right extrapolation', failures)
   end subroutine test_bs_reference

   subroutine test_linear(failures)
      integer, intent(inout) :: failures
      real(dp), allocatable :: y0(:)
      integer :: status
      call linear_interp([2.0_dp,0.0_dp,1.0_dp], [4.0_dp,0.0_dp,1.0_dp], &
         [0.5_dp,1.5_dp,2.0_dp], y0, status)
      call check(status == 0, 'linear interpolation status', failures)
      call check(maxval(abs(y0 - [0.5_dp,2.5_dp,4.0_dp])) < 1.0e-14_dp, &
         'linear interpolation values', failures)
   end subroutine test_linear

   subroutine test_interpolation(failures)
      integer, intent(inout) :: failures
      real(dp), parameter :: x(6) = [0.0_dp,0.5_dp,1.0_dp,1.5_dp,2.0_dp,3.0_dp]
      real(dp), parameter :: y(6) = [0.0_dp,0.25_dp,1.0_dp,2.25_dp,4.0_dp,9.0_dp]
      type(b_spline_t) :: bsp
      type(poly_spline_t) :: poly
      real(dp), allocatable :: fitted(:), pfit(:)
      integer :: status

      call fit_interpolating_spline(x, y, bsp, status)
      call check(status == 0, 'natural interpolation fit', failures)
      fitted = bsp%evaluate(x, status=status)
      call check(status == 0, 'natural interpolation evaluation', failures)
      call check(maxval(abs(fitted-y)) < 2.0e-12_dp, 'interpolation reproduces data', failures)
      call check_close(bsp%evaluate(0.1_dp), 0.0298130841121495_dp, 2.0e-13_dp, &
         'natural interpolation reference', failures)
      call check_close(bsp%evaluate(x(1), 2), 0.0_dp, 2.0e-11_dp, &
         'left natural constraint', failures)
      call check_close(bsp%evaluate(x(6), 2), 0.0_dp, 2.0e-11_dp, &
         'right natural constraint', failures)
      call check_close(bsp%evaluate(-1.0_dp), bsp%evaluate(0.0_dp) - bsp%evaluate(0.0_dp,1), &
         2.0e-12_dp, 'natural linear extrapolation', failures)

      call to_polynomial_spline(bsp, poly, status)
      call check(status == 0, 'B-spline to polynomial conversion', failures)
      pfit = poly%evaluate([0.1_dp,0.75_dp,1.75_dp,2.75_dp])
      fitted = bsp%evaluate([0.1_dp,0.75_dp,1.75_dp,2.75_dp])
      call check(maxval(abs(pfit-fitted)) < 2.0e-11_dp, 'polynomial representation agrees', failures)
   end subroutine test_interpolation

   subroutine test_periodic(failures)
      integer, intent(inout) :: failures
      real(dp), parameter :: pi = acos(-1.0_dp)
      real(dp) :: x(8), y(8), v1, v2
      type(b_spline_t) :: bsp
      integer :: i, status
      do i = 1, 8
         x(i) = real(i-1, dp) * 2.0_dp*pi/8.0_dp
         y(i) = sin(x(i))
      end do
      call fit_periodic_spline(x, y, bsp, period=2.0_dp*pi, status=status)
      call check(status == 0, 'periodic interpolation fit', failures)
      call check(maxval(abs(bsp%evaluate(x)-y)) < 2.0e-12_dp, &
         'periodic interpolation reproduces data', failures)
      v1 = bsp%evaluate(0.37_dp)
      v2 = bsp%evaluate(0.37_dp + 4.0_dp*pi)
      call check_close(v1, v2, 2.0e-13_dp, 'periodic wrapping', failures)
      call check(abs(v1-sin(0.37_dp)) < 3.0e-3_dp, 'periodic sine approximation', failures)
   end subroutine test_periodic

   subroutine test_regression_bases(failures)
      integer, intent(inout) :: failures
      real(dp), parameter :: x(9) = [-2.0_dp,-1.0_dp,0.0_dp,0.5_dp,1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp]
      real(dp), allocatable :: b(:, :), n(:, :), ik(:)
      integer :: status

      call bs_basis(x, b, degree=3, df=5, boundary_knots=[0.0_dp,3.0_dp], &
         interior_knots=ik, status=status)
      call check(status == 0, 'bs basis status', failures)
      call check(size(b,1)==9 .and. size(b,2)==5, 'bs basis dimensions', failures)
      call check(size(ik)==2, 'bs automatic knot count', failures)

      call natural_spline_basis(x, n, df=4, boundary_knots=[0.0_dp,3.0_dp], &
         interior_knots=ik, status=status)
      call check(status == 0, 'natural basis status', failures)
      call check(size(n,1)==9 .and. size(n,2)==4, 'natural basis dimensions', failures)
      ! Outside each boundary every natural-basis column is exactly linear.
      call check(maxval(abs(n(1,:) - (2.0_dp*n(2,:) - n(3,:)))) < 2.0e-11_dp, &
         'natural basis left linear extrapolation', failures)
      call check(maxval(abs(n(9,:) - (2.0_dp*n(8,:) - n(7,:)))) < 2.0e-11_dp, &
         'natural basis right linear extrapolation', failures)
   end subroutine test_regression_bases

   subroutine test_inverse(failures)
      integer, intent(inout) :: failures
      type(poly_spline_t) :: p, inv
      integer :: status
      p%knots = [0.0_dp,1.0_dp,2.0_dp]
      allocate(p%coefficients(3,4)); p%coefficients = 0.0_dp
      p%coefficients(:,1) = p%knots
      p%coefficients(:,2) = 1.0_dp
      p%natural = .true.
      call inverse_monotone_spline(p, inv, status)
      call check(status == 0, 'inverse spline construction', failures)
      call check_close(inv%evaluate(0.25_dp), 0.25_dp, 2.0e-13_dp, &
         'inverse identity lower', failures)
      call check_close(inv%evaluate(1.75_dp), 1.75_dp, 2.0e-13_dp, &
         'inverse identity upper', failures)
   end subroutine test_inverse

end program test_splines
