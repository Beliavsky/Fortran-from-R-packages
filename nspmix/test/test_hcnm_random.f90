program test_hcnm_random
   use nspmix
   implicit none
   integer, parameter :: nr=80,n=7
   real(dp) :: d(n,2),w(n),p0(2),pstar,llstar
   type(hcnm_result) :: r
   integer :: seed,i,j
   seed=24681357
   do j=1,nr
      do i=1,n
         d(i,1)=0.05_dp+0.9_dp*urand(seed)
         d(i,2)=0.05_dp+0.9_dp*urand(seed)
         w(i)=0.5_dp+4.0_dp*urand(seed)
      end do
      p0=[0.5_dp,0.5_dp]
      call hcnm(d,p0,w,r,maxit=500,tol=1.0e-12_dp)
      call scalar_opt(d,w,pstar,llstar)
      if(llstar-r%ll>2.0e-8_dp) error stop "random hcnm mismatch"
   end do
   print *, "test_hcnm_random: PASS", nr
contains
   real(dp) function urand(s)
      integer,intent(inout)::s
      integer(kind=8)::z
      z=mod(1103515245_8*int(s,8)+12345_8,2147483647_8); s=int(z); urand=real(s,dp)/2147483647.0_dp
   end function
   real(dp) function f(p,d,w)
      real(dp),intent(in)::p,d(:,:),w(:)
      f=sum(w*log(max(p*d(:,1)+(1.0_dp-p)*d(:,2),1.0e-300_dp)))
   end function
   subroutine scalar_opt(d,w,pbest,fbest)
      real(dp),intent(in)::d(:,:),w(:); real(dp),intent(out)::pbest,fbest
      real(dp)::a,b,c,e,fc,fe,gr
      integer::it
      gr=(sqrt(5.0_dp)-1.0_dp)/2.0_dp; a=0.0_dp; b=1.0_dp
      c=b-gr*(b-a); e=a+gr*(b-a); fc=f(c,d,w); fe=f(e,d,w)
      do it=1,120
         if(fc>fe) then; b=e; e=c; fe=fc; c=b-gr*(b-a); fc=f(c,d,w)
         else; a=c; c=e; fc=fe; e=a+gr*(b-a); fe=f(e,d,w); end if
      end do
      pbest=0.5_dp*(a+b); fbest=max(f(pbest,d,w),f(0.0_dp,d,w),f(1.0_dp,d,w))
   end subroutine
end program
