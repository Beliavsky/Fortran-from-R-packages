program test_links_thresholds
   use ordinal, only : dp, link_logit, link_probit, link_cloglog, link_loglog, link_cauchit, &
      link_cdf, link_pdf, link_pdf_gradient, pgumbel, dgumbel, qgumbel, &
      threshold_flexible, threshold_symmetric, threshold_symmetric2, threshold_equidistant, &
      threshold_jacobian, thresholds_from_alpha, threshold_start
   implicit none
   real(dp), allocatable :: tjac(:, :), alpha(:), theta(:)
   real(dp) :: p, q
   integer :: status

   call assert_close(link_cdf(0.0_dp, link_logit, 1.0_dp, .true.), 0.5_dp, 1.0e-14_dp, 'logit cdf')
   call assert_close(link_pdf(0.0_dp, link_logit, 1.0_dp), 0.25_dp, 1.0e-14_dp, 'logit pdf')
   call assert_close(link_pdf_gradient(0.0_dp, link_logit, 1.0_dp), 0.0_dp, 1.0e-14_dp, 'logit gradient')
   call assert_close(link_cdf(0.0_dp, link_probit, 1.0_dp, .true.), 0.5_dp, 1.0e-14_dp, 'probit cdf')
   call assert_close(link_cdf(0.0_dp, link_cauchit, 1.0_dp, .true.), 0.5_dp, 1.0e-14_dp, 'cauchit cdf')
   call assert_close(link_cdf(0.0_dp, link_cloglog, 1.0_dp, .true.), 1.0_dp - exp(-1.0_dp), 1.0e-14_dp, 'cloglog cdf')
   call assert_close(link_cdf(0.0_dp, link_loglog, 1.0_dp, .true.), exp(-1.0_dp), 1.0e-14_dp, 'loglog cdf')

   p = pgumbel(0.7_dp, 0.2_dp, 1.3_dp, .true., .true.)
   q = qgumbel(p, 0.2_dp, 1.3_dp, .true., .true.)
   call assert_close(q, 0.7_dp, 2.0e-13_dp, 'gumbel quantile inverse')
   if (dgumbel(0.0_dp, 0.0_dp, 1.0_dp, .false., .true.) <= 0.0_dp) error stop 'gumbel density'

   call threshold_jacobian(5, threshold_flexible, tjac, status)
   if (status /= 0 .or. size(tjac, 1) /= 4 .or. size(tjac, 2) /= 4) error stop 'flexible threshold jacobian'
   call threshold_start(5, threshold_equidistant, alpha, status)
   if (status /= 0 .or. size(alpha) /= 2) error stop 'equidistant start'
   call thresholds_from_alpha(alpha, 5, threshold_equidistant, theta, status)
   if (status /= 0 .or. any(theta(2:) <= theta(:size(theta) - 1))) error stop 'equidistant thresholds'
   call threshold_start(5, threshold_symmetric, alpha, status)
   if (status /= 0) error stop 'symmetric start'
   call thresholds_from_alpha(alpha, 5, threshold_symmetric, theta, status)
   if (status /= 0 .or. any(theta(2:) <= theta(:size(theta) - 1))) error stop 'symmetric thresholds'
   call threshold_start(5, threshold_symmetric2, alpha, status)
   if (status /= 0) error stop 'symmetric2 start'
   call thresholds_from_alpha(alpha, 5, threshold_symmetric2, theta, status)
   if (status /= 0 .or. any(theta(2:) <= theta(:size(theta) - 1))) error stop 'symmetric2 thresholds'

   print *, 'test_links_thresholds: PASS'
contains
   subroutine assert_close(actual, expected, tol, label)
      real(dp), intent(in) :: actual !! Computed scalar result to validate.
      real(dp), intent(in) :: expected !! Deterministic reference value.
      real(dp), intent(in) :: tol !! Maximum permitted absolute error.
      character(len=*), intent(in) :: label !! Human-readable assertion name used on failure.
      if (abs(actual - expected) > tol) then
         print *, trim(label), actual, expected
         error stop 1
      end if
   end subroutine assert_close
end program test_links_thresholds
