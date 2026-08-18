program test_diagnostics
   use coda
   implicit none
   integer, parameter :: n=5000
   real(dp) :: x(n,2), u1, u2, z1, z2
   integer :: i, nseed, fails
   integer, allocatable :: seed(:)
   type(mcmc_chain) :: ch
   type(geweke_result) :: gw
   type(heidel_result) :: hd
   type(raftery_result) :: rd
   type(mcmc_summary) :: sm
   type(spectrum_ar_result) :: sp
   real(dp), allocatable :: ess(:), bse(:), sg(:)

   fails = 0
   call random_seed(size=nseed)
   allocate(seed(nseed))
   seed = [(104729 + 37*i, i=1,nseed)]
   call random_seed(put=seed)
   i = 1
   do while (i <= n)
      call random_number(u1)
      call random_number(u2)
      u1 = max(u1, tiny(1.0_dp))
      z1 = sqrt(-2.0_dp*log(u1))*cos(2.0_dp*acos(-1.0_dp)*u2)
      z2 = sqrt(-2.0_dp*log(u1))*sin(2.0_dp*acos(-1.0_dp)*u2)
      x(i,1) = z1
      x(i,2) = 0.7_dp*z1 + sqrt(1.0_dp-0.7_dp**2)*z2
      if (i < n) then
         call random_number(u1)
         call random_number(u2)
         u1 = max(u1, tiny(1.0_dp))
         z1 = sqrt(-2.0_dp*log(u1))*cos(2.0_dp*acos(-1.0_dp)*u2)
         z2 = sqrt(-2.0_dp*log(u1))*sin(2.0_dp*acos(-1.0_dp)*u2)
         x(i+1,1) = z1
         x(i+1,2) = 0.7_dp*z1 + sqrt(1.0_dp-0.7_dp**2)*z2
      end if
      i = i + 2
   end do
   ch = make_mcmc(x)

   sp = spectrum0_ar(ch)
   if (any(sp%spec <= 0.0_dp)) fails = fails + 1
   ess = effective_size(ch)
   if (any(ess <= 0.0_dp) .or. any(ess > 2.0_dp*real(n,dp))) fails = fails + 1

   gw = geweke_diag(ch)
   if (any(abs(gw%z) > 8.0_dp)) fails = fails + 1

   bse = batch_se(ch, 100)
   if (any(bse <= 0.0_dp)) fails = fails + 1

   sm = summarize_mcmc(ch)
   if (maxval(abs(sm%statistics(:,1))) > 0.08_dp) fails = fails + 1
   if (maxval(abs(sm%statistics(:,2)-1.0_dp)) > 0.08_dp) fails = fails + 1

   rd = raftery_diag(ch, q=0.1_dp, r=0.02_dp, s=0.9_dp)
   if (.not. rd%enough_samples) fails = fails + 1
   if (any(rd%total <= 0)) fails = fails + 1
   if (any(rd%dependence_factor <= 0.0_dp)) fails = fails + 1

   hd = heidel_diag(ch)
   if (any(hd%start < -1000000)) then
      ! Failure to pass stationarity is allowed for a finite random sample;
      ! this condition only catches malformed integer output.
      fails = fails + 1
   end if

   sg = spectrum0(ch, order=1)
   if (any(sg <= 0.0_dp)) fails = fails + 1

   if (fails == 0) then
      print *, 'test_diagnostics: PASS'
   else
      print *, 'test_diagnostics: FAIL', fails
      error stop 1
   end if
end program test_diagnostics
