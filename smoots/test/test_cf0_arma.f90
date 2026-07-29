! SPDX-License-Identifier: GPL-3.0-only
program test_cf0_arma
   use smoots
   implicit none
   integer,parameter::n=2000
   real(dp)::x(n),cf0
   real(dp),allocatable::series(:),mat(:,:),psi(:)
   integer::i,status,l0,lg,p,q
   type(arma_model)::fit
   do i=1,100
      x(i)=sin(0.17_dp*real(i,dp))+0.2_dp*cos(0.71_dp*real(i,dp))
   end do
   call lag_window_variance(x(1:100),cf0,l0,lg,status)
   call assert_close(cf0,2.9741584433486246_dp,2.0e-12_dp,'cf0')
   if(l0/=6.or.lg/=5)error stop 'lag-window orders'
   call seed_rng(123456789_8)
   call simulate_arma([0.65_dp],[-0.35_dp],2.0_dp,n,500,series,status)
   call fit_arma(series,1,1,.true.,fit,status)
   if(status/=sm_ok.and.status/=sm_iteration_limit)error stop 'ARMA fit status'
   if(abs(fit%ar(1)-0.65_dp)>0.15_dp)error stop 'AR recovery'
   if(abs(fit%ma(1)+0.35_dp)>0.20_dp)error stop 'MA recovery'
   call ma_infinity([0.5_dp],[0.2_dp],4,psi,status)
   call assert_close(psi(1),1.0_dp,1.0e-14_dp,'psi0')
   call assert_close(psi(2),0.7_dp,1.0e-14_dp,'psi1')
   call assert_close(psi(3),0.35_dp,1.0e-14_dp,'psi2')
   call information_criterion_matrix(series(1:250),2,2,.true.,.true.,mat,status)
   call optimal_order(mat,p,q,status=status)
   if(p<0.or.q<0)error stop 'order selection'
   print '(a)','test_cf0_arma: PASS'
contains
   subroutine assert_close(a,b,t,msg)
      real(dp),intent(in)::a,b,t
      character(len=*),intent(in)::msg
      if(abs(a-b)>t)then;print *,trim(msg),a,b;error stop 1;end if
   end subroutine assert_close
end program test_cf0_arma
