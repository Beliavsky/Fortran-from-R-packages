program basic_diagnostics
   use coda
   implicit none
   integer, parameter :: n = 1000
   real(dp) :: x(n,2)
   type(mcmc_chain) :: chain
   type(geweke_result) :: gw
   type(mcmc_summary) :: sm
   real(dp), allocatable :: ess(:), hpd(:,:)
   integer :: i

   do i = 1, n
      x(i,1) = sin(0.017_dp*real(i,dp)) + 0.2_dp*cos(0.113_dp*real(i,dp))
      x(i,2) = cos(0.023_dp*real(i,dp)) + 0.1_dp*sin(0.071_dp*real(i,dp))
   end do

   chain = make_mcmc(x, var_names=[character(len=5) :: 'alpha','beta '])
   sm = summarize_mcmc(chain)
   ess = effective_size(chain)
   hpd = hpd_interval(chain)
   gw = geweke_diag(chain)

   print '(a)', 'Variable      Mean          SD           ESS       Geweke z'
   do i = 1, chain%nvar()
      print '(a8,4(1x,f12.5))', chain%var_names(i), sm%statistics(i,1), &
         sm%statistics(i,2), ess(i), gw%z(i)
      print '(a,2f12.5)', '  95% HPD: ', hpd(i,1), hpd(i,2)
   end do
end program basic_diagnostics
