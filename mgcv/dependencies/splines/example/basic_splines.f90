program basic_splines
   use splines, only : dp, b_spline_t, poly_spline_t, &
      fit_interpolating_spline, to_polynomial_spline, bs_basis
   implicit none

   real(dp), parameter :: x(6) = [0.0_dp, 0.5_dp, 1.0_dp, 1.5_dp, 2.0_dp, 3.0_dp]
   real(dp), parameter :: y(6) = [0.0_dp, 0.25_dp, 1.0_dp, 2.25_dp, 4.0_dp, 9.0_dp]
   real(dp), parameter :: x_new(5) = [0.0_dp, 0.75_dp, 1.25_dp, 2.25_dp, 3.0_dp]
   real(dp), allocatable :: fitted(:), basis(:, :)
   type(b_spline_t) :: spline
   type(poly_spline_t) :: polynomial
   integer :: i, status

   call fit_interpolating_spline(x, y, spline, status)
   if (status /= 0) error stop 'fit_interpolating_spline failed'

   fitted = spline%evaluate(x_new, status=status)
   if (status /= 0) error stop 'spline evaluation failed'

   print '(a)', 'Natural cubic interpolation:'
   do i = 1, size(x_new)
      print '(f7.3,2x,f12.6)', x_new(i), fitted(i)
   end do

   call to_polynomial_spline(spline, polynomial, status)
   if (status /= 0) error stop 'to_polynomial_spline failed'
   print '(a,f12.6)', 'Polynomial representation at x=1.25: ', &
      polynomial%evaluate(1.25_dp)

   call bs_basis(x_new, basis, degree=3, knots=[1.0_dp, 2.0_dp], &
      boundary_knots=[0.0_dp, 3.0_dp], intercept=.true., status=status)
   if (status /= 0) error stop 'bs_basis failed'
   print '(a,i0,a,i0)', 'B-spline design matrix dimensions: ', &
      size(basis, 1), ' x ', size(basis, 2)
end program basic_splines
