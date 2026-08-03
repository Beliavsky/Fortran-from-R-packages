! SPDX-License-Identifier: GPL-3.0-only
program test_tail_estimators
   use fitheavytail
   use test_support, only: check, make_data
   implicit none
   real(dp) :: x(240,3), varx(3), r2(240), nu
   integer :: i, status

   call make_data(x)
   call check(abs(log_bessel_k(0.01_dp,5.0_dp) - &
      28.9764872325347_dp) < 1.0e-9_dp,'log Bessel K')
   call check(abs(bessel_k_ratio(10.0_dp,70.0_dp) - &
      14.0720817849628_dp) < 1.0e-9_dp,'Bessel K ratio')

   varx=[1.5_dp,1.5_dp,1.5_dp]
   nu=nu_opp_estimator(varx,3.0_dp,nu_min=2.5_dp, &
      nu_max=100.0_dp,status=status)
   call check(status==ht_success,'OPP status')
   call check(abs(nu-6.0_dp)<1.0e-12_dp,'OPP value')

   do i=1,size(x,1)
   r2(i)=sum((x(i,:)-sum(x,dim=1)/240.0_dp)**2)
   end do
   nu=nu_pop_estimator(r2,7.0_dp,3,'POP-approx-2', &
      nu_min=2.5_dp,nu_max=100.0_dp,status=status)
   call check(status==ht_success,'POP status')
   call check(nu>=2.5_dp.and.nu<=100.0_dp,'POP bounds')
   call check(nu_from_average_marginal_kurtosis(x)>=2.5_dp, &
      'marginal kurtosis estimator')
   call check(nu_from_cross_cumulants(x)>=2.5_dp, &
      'cross-cumulant estimator')
   call check(nu_hill_estimator(x)>=2.5_dp,'Hill estimator')
   call check(nu_pareto_tail_index(x)>=2.5_dp, &
      'Pareto estimator')
   write(*,'(a)') 'test_tail_estimators: PASS'
end program test_tail_estimators
