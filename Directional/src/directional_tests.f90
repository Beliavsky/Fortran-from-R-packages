module directional_tests
   use directional_kinds, only : dp, pi
   implicit none
   private
   type, public :: test_result
      real(dp) :: statistic=0.0_dp, p_value=1.0_dp
      integer :: df=0
   end type
   public :: rayleigh_test, kuiper_test
contains
   function rayleigh_test(x,modified,mc_reps) result(res)
      real(dp),intent(in)::x(:,:);logical,intent(in),optional::modified;integer,intent(in),optional::mc_reps
      type(test_result)::res
      logical::modif;integer::b,i,j,n,p,exceed;real(dp)::m(size(x,2)),z(size(x,2)),nr,tb
      modif=.true.;if(present(modified))modif=modified;b=1;if(present(mc_reps))b=mc_reps
      n=size(x,1);p=size(x,2);m=sum(x,dim=1);res%statistic=p*sum(m*m)/n
      if(modif)res%statistic=(1-1/(2.0_dp*n))*res%statistic+res%statistic**2/(2.0_dp*n*(p+2.0_dp));res%df=p
      if(b<=1)then;res%p_value=gamma_q(0.5_dp*p,0.5_dp*res%statistic);return;end if
      exceed=0
      do i=1,b
         m=0
         do j=1,n;call random_normal_vec(z);nr=sqrt(sum(z*z));m=m+z/nr;end do
         tb=p*sum(m*m)/n;if(tb>res%statistic)exceed=exceed+1
      end do
      res%p_value=real(exceed+1,dp)/real(b+1,dp);res%df=0
   end function

   function kuiper_test(u,rads,mc_reps) result(res)
      real(dp),intent(in)::u(:);logical,intent(in),optional::rads;integer,intent(in),optional::mc_reps
      type(test_result)::res
      real(dp)::x(size(u)),sim(size(u)),f,vn,a1,a2,pv;integer::n,i,j,b,exceed
      logical::rr;rr=.false.;if(present(rads))rr=rads;x=u;if(.not.rr)x=x*pi/180;x=modulo(x,2*pi)/(2*pi);call sort_real(x)
      n=size(x);f=sqrt(real(n,dp));vn=f*(maxval(x-[(real(i-1,dp)/n,i=1,n)])+maxval([(real(i,dp)/n,i=1,n)]-x));res%statistic=vn
      b=1;if(present(mc_reps))b=mc_reps
      if(b<=1)then
         pv=0
         do j=1,50;a1=4.0_dp*j*j*vn*vn;a2=exp(-2.0_dp*j*j*vn*vn);pv=pv+2*(a1-1)*a2-8*vn/(3*f)*j*j*(a1-3)*a2;end do
         res%p_value=max(0.0_dp,min(1.0_dp,pv));return
      end if
      exceed=0
      do j=1,b
         call random_number(sim);call sort_real(sim)
         pv=f*(maxval(sim-[(real(i-1,dp)/n,i=1,n)])+maxval([(real(i,dp)/n,i=1,n)]-sim));if(pv>vn)exceed=exceed+1
      end do
      res%p_value=real(exceed+1,dp)/real(b+1,dp)
   end function

   subroutine random_normal_vec(z)
      real(dp),intent(out)::z(:);real(dp)::u1,u2;integer::i
      i=1;do while(i<=size(z));call random_number(u1);call random_number(u2);u1=max(u1,tiny(1.0_dp));z(i)=sqrt(-2*log(u1))*cos(2*pi*u2);if(i+1<=size(z))z(i+1)=sqrt(-2*log(u1))*sin(2*pi*u2);i=i+2;end do
   end subroutine
   subroutine sort_real(x)
      real(dp),intent(inout)::x(:);integer::i,j;real(dp)::t
      do i=2,size(x);t=x(i);j=i-1;do while(j>=1)
         if(x(j)<=t) exit
         x(j+1)=x(j);j=j-1
      end do
      x(j+1)=t
   end do
   end subroutine
   pure real(dp) function gamma_q(a,x) result(q)
      real(dp),intent(in)::a,x;real(dp)::ap,del,s,b,c,d,h,an;integer::i
      if(x<=0)then;q=1;return;end if
      if(x<a+1)then
         ap=a;del=1/a;s=del;do i=1,500;ap=ap+1;del=del*x/ap;s=s+del;if(abs(del)<abs(s)*1e-14_dp)exit;end do
         q=1-s*exp(-x+a*log(x)-log_gamma(a))
      else
         b=x+1-a;c=1/tiny(1.0_dp);d=1/b;h=d
         do i=1,500;an=-real(i,dp)*(real(i,dp)-a);b=b+2;d=an*d+b;if(abs(d)<tiny(1.0_dp))d=tiny(1.0_dp);c=b+an/c;if(abs(c)<tiny(1.0_dp))c=tiny(1.0_dp);d=1/d;del=d*c;h=h*del;if(abs(del-1)<1e-14_dp)exit;end do
         q=exp(-x+a*log(x)-log_gamma(a))*h
      end if
      q=max(0.0_dp,min(1.0_dp,q))
   end function
end module directional_tests
