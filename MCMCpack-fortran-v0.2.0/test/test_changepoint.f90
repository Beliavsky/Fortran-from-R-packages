program test_changepoint
   use mcmcpack
   implicit none
   integer,parameter::n=12,ns=2,k=2
   integer::yb(n),yp(n),st(n),i
   real(dp)::phi(ns),p(ns,ns),a0(ns,ns),x(n,k),bs(ns,k),b0(k),bp(k,k),yr(n),sig(ns)
   type(change_result)::r
   call set_seed(97531)
   yb=[0,0,0,1,0,0,1,1,1,1,0,1]
   phi=[0.2_dp,0.8_dp];p=0.0_dp;p(1,1)=0.9_dp;p(1,2)=0.1_dp;p(2,2)=1.0_dp
   a0=0.0_dp;a0(1,1)=5.0_dp;a0(1,2)=1.0_dp;a0(2,2)=1.0_dp
   r=mcmc_binary_change(yb,phi,p,a0,1.0_dp,1.0_dp,10,20,2)
   if(r%status/=0.or.size(r%draws,1)/=10.or.size(r%states,2)/=n) error stop 'binary change'
   if(any(r%states(:,n)/=ns)) error stop 'binary final state'

   do i=1,n;x(i,1)=1.0_dp;x(i,2)=real(i-6,dp)/5.0_dp;end do
   bs=0.0_dp;bs(1,1)=-0.5_dp;bs(2,1)=0.5_dp;b0=0.0_dp;bp=0.0_dp;bp(1,1)=0.2_dp;bp(2,2)=0.2_dp
   r=mcmc_probit_change(yb,x,bs,p,b0,bp,a0,10,20,2)
   if(r%status/=0.or.size(r%draws,2)/=ns*k+ns*ns) error stop 'probit change'

   st(1:6)=1;st(7:n)=2;sig=1.0_dp
   do i=1,n
      if(i<=6)then;yr(i)=0.5_dp+0.3_dp*x(i,2);else;yr(i)=1.5_dp-0.2_dp*x(i,2);end if
   end do
   bs=0.0_dp
   r=mcmc_regress_change(yr,x,bs,sig,p,st,b0,bp,2.0_dp,1.0_dp,a0,10,20,2)
   if(r%status/=0.or.size(r%draws,2)/=ns*k+ns+ns*ns) error stop 'regress change'
   if(any(r%prob_state<0.0_dp).or.any(r%prob_state>1.0_dp)) error stop 'change probabilities'

   yp=[0,1,0,1,1,2,3,2,4,3,5,4];bs=0.0_dp
   r=mcmc_poisson_change(yp,x,bs,p,b0,bp,a0,10,20,2,0.15_dp)
   if(r%status/=0.or.size(r%draws,2)/=ns*k+ns*ns) error stop 'poisson change'
   print '(a)','test_changepoint: PASS'
end program
