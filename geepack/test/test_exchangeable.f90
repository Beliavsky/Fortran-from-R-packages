program test_exchangeable
   use geepack
   implicit none
   real(dp) :: y(24)
   real(dp) :: x(24, 2)
   integer :: cs(8)
   real(dp) :: b0(2)
   type(gee_spec) :: spec
   type(gee_result) :: fit
   integer :: g
   integer :: j
   integer :: i

   cs = 3
   x(:, 1) = 1.0_dp
   i = 0
   do g = 1, 8
      do j = 1, 3
         i = i + 1
         x(i, 2) = real(j - 2, dp)
         y(i) = 1.0_dp + 0.6_dp * x(i, 2) + 0.18_dp * real(g - 4, dp) + &
            0.04_dp * real(mod(g + j, 3) - 1, dp)
      end do
   end do
   b0 = [1.0_dp, 0.5_dp]
   spec%corstr = COR_EXCHANGEABLE
   spec%scale_fixed = .true.
   spec%tolerance = 1.0e-9_dp
   spec%approximate_jackknife = .true.
   spec%one_step_jackknife = .true.
   spec%fully_iterated_jackknife = .true.
   allocate(spec%mean_links(1), spec%variance_codes(1), spec%scale_links(1))
   spec%mean_links = LINK_IDENTITY
   spec%variance_codes = VAR_GAUSSIAN
   spec%scale_links = LINK_IDENTITY
   call fit_geese(y, x, cs, b0, spec, fit)
   if (fit%error /= GEE_OK) error stop 1
   if (size(fit%alpha) /= 1) error stop 2
   if (fit%alpha(1) <= -0.5_dp .or. fit%alpha(1) >= 1.0_dp) error stop 3
   if (maxval(abs(fit%beta - [1.09_dp, 0.6_dp])) > 0.12_dp) error stop 4
   if (size(fit%influence, 2) /= 8) error stop 5
   if (any([(fit%vbeta_ajs(i, i) < 0.0_dp, i=1,2)])) error stop 6
   if (any([(fit%vbeta_j1s(i, i) < 0.0_dp, i=1,2)])) error stop 7
   if (any([(fit%vbeta_fij(i, i) < 0.0_dp, i=1,2)])) error stop 8
   print *, 'test_exchangeable: PASS'
end program test_exchangeable
