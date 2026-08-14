program test_cnm
   use nspmix
   implicit none
   integer,parameter::n=13
   real(dp)::v(n),w(n),lam,p,pmf
   integer::i,j
   type(nsp_data)::x
   type(nspmix_result)::r
   type(disc_dist)::init
   v=[(real(i-1,dp),i=1,n)]; w=0.0_dp
   do i=1,n
      do j=1,2
         if(j==1) then; lam=1.5_dp; p=0.65_dp; else; lam=6.0_dp; p=0.35_dp; end if
         pmf=exp(v(i)*log(lam)-lam-log_gamma(v(i)+1.0_dp)); w(i)=w(i)+10000.0_dp*p*pmf
      end do
   end do
   call make_nppois_data(v,x,w)
   call make_disc([0.8_dp,3.0_dp,8.0_dp],d=init)
   call cnm(x,r,init_mix=init,maxit=30,tol=1.0e-8_dp,ngrid=120,kmax=12)
   if(abs(sum(r%mix%pr)-1.0_dp)>1.0e-11_dp) error stop "cnm masses"
   if(minval(r%mix%pt)<0.0_dp) error stop "cnm support"
   if(.not.(r%ll<0.0_dp)) error stop "cnm ll"
   if(size(r%mix%pt)>12) error stop "cnm kmax"
   print *, "test_cnm: PASS", r%mix%pt, r%mix%pr
end program
