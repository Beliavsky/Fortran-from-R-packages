program test_core
   use coda
   implicit none
   real(dp) :: x(10,2), y(6,1)
   real(dp), allocatable :: ac(:,:,:), cc(:,:), hpd(:,:), rr(:)
   type(mcmc_chain) :: ch, ch2
   integer :: i, fails

   fails = 0
   do i = 1, 10
      x(i,1) = real(i,dp)
      x(i,2) = 2.0_dp*real(i,dp)
   end do
   ch = make_mcmc(x, start=5, thin=2, var_names=[character(len=2) :: 'x1','x2'])
   if (ch%niter() /= 10 .or. ch%nvar() /= 2) fails = fails + 1
   if (ch%finish /= 23) fails = fails + 1

   ac = autocorr(ch, [0,1])
   if (abs(ac(1,1,1)-1.0_dp) > 1.0e-12_dp) fails = fails + 1
   if (abs(ac(1,1,2)-1.0_dp) > 1.0e-12_dp) fails = fails + 1

   cc = crosscorr(ch)
   if (maxval(abs(cc-1.0_dp)) > 1.0e-12_dp) fails = fails + 1

   hpd = hpd_interval(ch, 0.5_dp)
   if (abs(hpd(1,1)-1.0_dp) > 1.0e-12_dp .or. abs(hpd(1,2)-6.0_dp) > 1.0e-12_dp) fails = fails + 1

   y(:,1) = [1.0_dp,1.0_dp,2.0_dp,2.0_dp,2.0_dp,3.0_dp]
   ch2 = make_mcmc(y)
   rr = rejection_rate(ch2)
   if (abs(rr(1)-0.6_dp) > 1.0e-12_dp) fails = fails + 1

   if (abs(cramer_cdf(0.1_dp)-0.415126561593196_dp) > 2.0e-7_dp) fails = fails + 1
   if (abs(cramer_cdf(1.0_dp)-0.997539547819866_dp) > 2.0e-7_dp) fails = fails + 1

   if (fails == 0) then
      print *, 'test_core: PASS'
   else
      print *, 'test_core: FAIL', fails
      error stop 1
   end if
end program test_core
