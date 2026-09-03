program test_gaussian
   use geepack
   implicit none
   real(dp) :: y(12)
   real(dp) :: x(12, 2)
   integer :: cs(4)
   real(dp) :: b0(2)
   real(dp) :: xtx(2, 2)
   real(dp) :: xty(2)
   real(dp) :: bref(2)
   type(gee_spec) :: spec
   type(gee_result) :: fit
   integer :: i

   cs = 3
   y = [1.0_dp, 2.0_dp, 2.5_dp, 1.5_dp, 2.2_dp, 3.0_dp, 2.0_dp, 3.1_dp, 3.4_dp, &
      2.4_dp, 3.0_dp, 4.1_dp]
   x(:, 1) = 1.0_dp
   do i = 1, 12
      x(i, 2) = real(i - 1, dp) / 5.0_dp
   end do
   b0 = [1.0_dp, 1.0_dp]
   xtx = matmul(transpose(x), x)
   xty = matmul(transpose(x), y)
   bref(1) = (xty(1) * xtx(2, 2) - xty(2) * xtx(1, 2)) / &
      (xtx(1, 1) * xtx(2, 2) - xtx(1, 2) * xtx(2, 1))
   bref(2) = (xtx(1, 1) * xty(2) - xtx(2, 1) * xty(1)) / &
      (xtx(1, 1) * xtx(2, 2) - xtx(1, 2) * xtx(2, 1))
   spec%corstr = COR_INDEPENDENCE
   spec%scale_fixed = .true.
   allocate(spec%mean_links(1), spec%variance_codes(1), spec%scale_links(1))
   spec%mean_links = LINK_IDENTITY
   spec%variance_codes = VAR_GAUSSIAN
   spec%scale_links = LINK_IDENTITY
   call fit_geese(y, x, cs, b0, spec, fit)
   if (fit%error /= GEE_OK) error stop 1
   if (maxval(abs(fit%beta - bref)) > 1.0e-10_dp) error stop 2
   if (any([(fit%vbeta(i, i) <= 0.0_dp, i=1,2)])) error stop 3
   print *, 'test_gaussian: PASS'
end program test_gaussian
