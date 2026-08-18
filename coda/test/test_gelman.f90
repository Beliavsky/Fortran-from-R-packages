program test_gelman
   use coda
   implicit none
   integer, parameter :: n=100, m=4, p=2
   type(mcmc_chain) :: chains(m)
   type(mcmc_list) :: lst
   type(gelman_result) :: gr
   real(dp) :: x(n,p)
   integer :: i, k, fails
   real(dp), parameter :: ref(2,2) = reshape([ &
      0.9950963163687483_dp, 0.9967687139623238_dp, &
      0.9951865745386762_dp, 1.0005196418711937_dp], [2,2])

   fails = 0
   do k = 1, m
      do i = 1, n
         x(i,1) = sin(0.11_dp*real(i,dp)+0.2_dp*real(k,dp)) + 0.02_dp*real(k,dp)
         x(i,2) = cos(0.07_dp*real(i,dp)-0.15_dp*real(k,dp)) + 0.03_dp*real(k,dp)
      end do
      chains(k) = make_mcmc(x)
   end do
   lst = make_mcmc_list(chains)
   gr = gelman_diag(lst, autoburnin=.false.)
   if (maxval(abs(gr%psrf-ref)) > 2.0e-12_dp) fails = fails + 1
   if (.not. gr%has_mpsrf) fails = fails + 1
   if (abs(gr%mpsrf-0.9972575791017942_dp) > 2.0e-12_dp) fails = fails + 1

   if (fails == 0) then
      print *, 'test_gelman: PASS'
   else
      print *, 'test_gelman: FAIL', fails
      error stop 1
   end if
end program test_gelman
