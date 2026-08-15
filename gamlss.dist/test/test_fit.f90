program test_fit
   use gamlss_dist
   implicit none
   integer, parameter :: n = 8
   real(dp) :: y(n), x_mu(n,2), x_sigma(n,1), x(n), e(n)
   type(gamlss_fit_result_t) :: result

   x = [-3.5_dp, -2.5_dp, -1.5_dp, -0.5_dp, &
         0.5_dp,  1.5_dp,  2.5_dp,  3.5_dp]
   e = [1.0_dp, -1.0_dp, -1.0_dp, 1.0_dp, &
        1.0_dp, -1.0_dp, -1.0_dp, 1.0_dp]
   x_mu(:,1) = 1.0_dp
   x_mu(:,2) = x
   x_sigma(:,1) = 1.0_dp
   y = 1.2_dp - 0.7_dp*x + e

   call fit_gamlss(y, x_mu, GAMLSS_NO, result, x_sigma=x_sigma, &
      max_iter=500, tol=1.0e-9_dp)
   if (.not. result%converged) error stop 'normal fit did not converge'
   call assert_close(result%beta_mu(1), 1.2_dp, 2.0e-5_dp, 'normal intercept')
   call assert_close(result%beta_mu(2), -0.7_dp, 2.0e-5_dp, 'normal slope')
   call assert_close(exp(result%beta_sigma(1)), 1.0_dp, 2.0e-5_dp, 'normal sigma')
   if (size(result%covariance, 1) /= 3) error stop 'missing covariance'

   print '(a)', 'test_fit: PASS'
contains
   subroutine assert_close(actual, expected, tolerance, name)
      real(dp), intent(in) :: actual, expected, tolerance
      character(*), intent(in) :: name
      if (abs(actual - expected) > tolerance) then
         print '(a)', 'FAIL: '//trim(name)
         print '(a,es24.16)', ' actual   = ', actual
         print '(a,es24.16)', ' expected = ', expected
         error stop 1
      end if
   end subroutine assert_close
end program test_fit
