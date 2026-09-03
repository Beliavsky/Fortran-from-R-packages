program test_logit
   use geepack
   implicit none
   real(dp) :: y(24)
   real(dp) :: x(24, 2)
   real(dp) :: xv(24)
   integer :: cs(8)
   real(dp) :: b0(2)
   real(dp), parameter :: bref(2) = [0.102387559999_dp, 0.333436059999_dp]
   type(gee_spec) :: spec
   type(gee_result) :: fit

   cs = 3
   xv = [-1.2_dp, -0.8_dp, -0.4_dp, -0.1_dp, 0.3_dp, 0.7_dp, 1.0_dp, 1.3_dp, 1.6_dp, &
      -1.0_dp, -0.5_dp, 0.2_dp, 0.6_dp, 1.1_dp, 1.5_dp, -1.3_dp, -0.7_dp, 0.1_dp, &
      0.5_dp, 0.9_dp, 1.4_dp, -0.9_dp, -0.2_dp, 0.8_dp]
   y = [0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
      0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, &
      0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp]
   x(:, 1) = 1.0_dp
   x(:, 2) = xv
   b0 = 0.0_dp
   spec%corstr = COR_INDEPENDENCE
   spec%scale_fixed = .true.
   spec%tolerance = 1.0e-10_dp
   allocate(spec%mean_links(1), spec%variance_codes(1), spec%scale_links(1))
   spec%mean_links = LINK_LOGIT
   spec%variance_codes = VAR_BINOMIAL
   spec%scale_links = LINK_IDENTITY
   call fit_geese(y, x, cs, b0, spec, fit)
   if (fit%error /= GEE_OK) error stop 1
   if (maxval(abs(fit%beta - bref)) > 2.0e-8_dp) error stop 2
   if (fit%iterations < 2) error stop 3
   print *, 'test_logit: PASS'
end program test_logit
