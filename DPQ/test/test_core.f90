program test_core
   use r_compat, only: dp, r_lgamma, r_gamma, dpois, dbinom, dnbinom
   use dpq_core
   use dpq_gamma_discrete
   implicit none

   call check_close(log1mexp(-log(0.2_dp)), log(0.8_dp), 2.0e-14_dp, 'log1mexp')
   call check_close(logspace_add(log(0.2_dp), log(0.3_dp)), log(0.5_dp), 2.0e-14_dp, 'logspace_add')
   call check_close(logspace_sub(log(0.8_dp), log(0.3_dp)), log(0.5_dp), 2.0e-14_dp, 'logspace_sub')
   call check_close(lgamma1p(0.2_dp), r_lgamma(1.2_dp), 1.0e-9_dp, 'lgamma1p')
   call check_close(rexpm1(1.0e-6_dp), exp(1.0e-6_dp)-1.0_dp, 2.0e-15_dp, 'rexpm1')
   call check_close(rlog1(1.0e-6_dp), 1.0e-6_dp-log(1.0_dp+1.0e-6_dp), 2.0e-15_dp, 'rlog1')
   call check_close(dpois_raw(7.0_dp, 4.2_dp), dpois(7.0_dp,4.2_dp), 2.0e-14_dp, 'dpois_raw')
   call check_close(dgamma_r(2.3_dp,3.2_dp,0.7_dp), (2.3_dp**2.2_dp)*exp(-2.3_dp/0.7_dp)/(r_gamma(3.2_dp)*0.7_dp**3.2_dp), &
      2.0e-13_dp, 'dgamma_r')
   call check_close(dbinom_raw(4.0_dp,10.0_dp,0.3_dp,0.7_dp), dbinom(4.0_dp,10,0.3_dp), 2.0e-14_dp, 'dbinom_raw')
   call check_close(dnbinom_r(4.0_dp,7.0_dp,0.4_dp), dnbinom(4.0_dp,7,0.4_dp), 2.0e-13_dp, 'dnbinom_r')
   call check_close(bd0(10.0_dp,10.0_dp), 0.0_dp, 1.0e-15_dp, 'bd0 equality')
   call check_close(stirlerr(10.0_dp), r_lgamma(11.0_dp)-(10.5_dp)*log(10.0_dp)+10.0_dp &
      -0.5_dp*log(2.0_dp*acos(-1.0_dp)), 2.0e-14_dp, 'stirlerr')
   print '(a)', 'test_core: PASS'
contains
   subroutine check_close(got,want,tol,label)
      real(dp),intent(in)::got,want,tol
      character(len=*),intent(in)::label
      if (.not.(abs(got-want) <= tol*max(1.0_dp,abs(want)))) then
         write(*,'(a,2es24.15)') trim(label)//' failed: ',got,want
         error stop 1
      end if
   end subroutine
end program
