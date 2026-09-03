program test_scale
   use geepack
   implicit none
   real(dp) :: y(20)
   real(dp) :: x(20, 1)
   integer :: cs(10)
   real(dp) :: b0(1)
   real(dp) :: expected_scale
   type(gee_spec) :: spec
   type(gee_result) :: fit

   cs = 2
   x = 1.0_dp
   y = [0.6_dp, 1.3_dp, 0.8_dp, 1.2_dp, 0.5_dp, 1.4_dp, 0.7_dp, 1.1_dp, &
      0.9_dp, 1.5_dp, 0.4_dp, 1.0_dp, 0.8_dp, 1.6_dp, 0.6_dp, 1.2_dp, &
      0.7_dp, 1.3_dp, 0.5_dp, 1.1_dp]
   b0 = 1.0_dp
   spec%corstr = COR_INDEPENDENCE
   spec%scale_fixed = .false.
   spec%tolerance = 1.0e-10_dp
   allocate(spec%mean_links(1), spec%variance_codes(1), spec%scale_links(1))
   spec%mean_links = LINK_IDENTITY
   spec%variance_codes = VAR_GAUSSIAN
   spec%scale_links = LINK_IDENTITY
   call fit_geese(y, x, cs, b0, spec, fit)
   if (fit%error /= GEE_OK) error stop 1
   expected_scale = sum((y - sum(y) / real(size(y), dp)) ** 2) / real(size(y), dp)
   if (abs(fit%beta(1) - sum(y) / real(size(y), dp)) > 1.0e-10_dp) error stop 2
   if (abs(fit%gamma(1) - expected_scale) > 1.0e-9_dp) error stop 3
   if (fit%vgamma(1, 1) <= 0.0_dp) error stop 4
   print *, 'test_scale: PASS'
end program test_scale
