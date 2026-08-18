program test_multinomial
   use mcmcpack
   implicit none
   integer,parameter::n=8,c=3,k=2
   real(dp)::y(n,c),x(n*c,k),b(k),bm(k),b0(k),bprec(k,k),v(k,k),tune(k,k)
   type(mcmc_result)::r
   integer::i,j,row
   call set_seed(777)
   y=0.0_dp
   do i=1,n;y(i,1+mod(i-1,c))=1.0_dp;end do
   do i=1,n;do j=1,c
      row=(i-1)*c+j;x(row,1)=merge(0.0_dp,1.0_dp,j==1);x(row,2)=real(j-1,dp)*real(i-4,dp)/5.0_dp
   end do;end do
   b=0.0_dp;bm=0.0_dp;b0=0.0_dp;bprec=0.0_dp;bprec(1,1)=0.1_dp;bprec(2,2)=0.1_dp
   v=0.0_dp;v(1,1)=0.2_dp;v(2,2)=0.2_dp;tune=0.0_dp;tune(1,1)=0.5_dp;tune(2,2)=0.5_dp
   r=mcmc_mnl(y,x,b,bm,b0,bprec,v,tune,10,20,2,.false.,6.0_dp)
   if(r%status/=0.or.size(r%draws,1)/=10.or.size(r%draws,2)/=k) error stop 'mnl'
   if(r%accept_rate<0.0_dp.or.r%accept_rate>1.0_dp) error stop 'mnl accept'
   print '(a)','test_multinomial: PASS'
end program
