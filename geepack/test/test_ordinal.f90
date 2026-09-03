program test_ordinal
   use geepack
   implicit none
   integer :: y(30)
   real(dp) :: x(30, 1)
   integer :: cs(10)
   type(ordinal_spec) :: spec
   type(ordinal_result) :: fit
   type(ordinal_result) :: fit_exchangeable
   type(ordinal_result) :: fit_reverse
   real(dp) :: psi
   real(dp) :: mu1
   real(dp) :: mu2
   real(dp) :: h
   real(dp) :: fd
   real(dp) :: deriv(2)
   real(dp) :: fd1
   real(dp) :: fd2
   integer :: i

   psi = 1.7_dp
   mu1 = 0.35_dp
   mu2 = 0.62_dp
   h = 1.0e-6_dp
   fd = (odds_to_p11(psi + h, mu1, mu2) - odds_to_p11(psi - h, mu1, mu2)) / (2.0_dp * h)
   if (abs(fd - p11_odds_derivative(psi, mu1, mu2)) > 2.0e-8_dp) error stop 1
   call p11_mean_derivatives(psi, mu1, mu2, deriv)
   fd1 = (odds_to_p11(psi, mu1 + h, mu2) - odds_to_p11(psi, mu1 - h, mu2)) / (2.0_dp * h)
   fd2 = (odds_to_p11(psi, mu1, mu2 + h) - odds_to_p11(psi, mu1, mu2 - h)) / (2.0_dp * h)
   if (maxval(abs(deriv - [fd1, fd2])) > 2.0e-8_dp) error stop 2

   cs = 3
   y = [1, 2, 1, 2, 2, 3, 1, 1, 2, 2, 3, 3, 1, 2, 3, &
      2, 1, 3, 1, 2, 2, 2, 3, 1, 3, 2, 1, 3, 3, 2]
   do i = 1, 30
      x(i, 1) = real(mod(i - 1, 5) - 2, dp) / 2.0_dp
   end do
   spec%mean_link = LINK_LOGIT
   spec%corstr = COR_INDEPENDENCE
   spec%constant_intercepts = .true.
   spec%tolerance = 1.0e-7_dp
   spec%max_iterations = 50
   call fit_ordgee(y, x, cs, 3, spec, fit)
   if (fit%error /= GEE_OK) error stop 3
   if (size(fit%beta) /= 3 .or. size(fit%alpha) /= 0) error stop 4
   if (any(abs(fit%beta) > 20.0_dp)) error stop 5
   if (any([(fit%vbeta(i, i) <= 0.0_dp, i=1,3)])) error stop 6
   spec%corstr = COR_EXCHANGEABLE
   call fit_ordgee(y, x, cs, 3, spec, fit_exchangeable)
   if (fit_exchangeable%error /= GEE_OK .or. size(fit_exchangeable%alpha) /= 1) error stop 7
   if (abs(fit_exchangeable%alpha(1)) > 10.0_dp) error stop 8
   spec%corstr = COR_INDEPENDENCE
   spec%reverse_coding = .true.
   spec%constant_intercepts = .false.
   call fit_ordgee(y, x, cs, 3, spec, fit_reverse)
   if (fit_reverse%error /= GEE_OK .or. size(fit_reverse%beta) /= 7) error stop 9
   print *, 'test_ordinal: PASS'
end program test_ordinal
