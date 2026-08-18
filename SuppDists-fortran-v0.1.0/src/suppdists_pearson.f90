module suppdists_pearson
   use suppdists_kinds, only : dp, sqrt2
   use suppdists_special, only : normal_pdf, adaptive_integral, randn
   use suppdists_stats, only : dist_stats
   implicit none
   private
   public :: dpearson, ppearson, qpearson, rpearson, spearson
contains
   pure real(dp) function hyper2f1_half(c,x) result(sumv)
      real(dp),intent(in)::c,x
      real(dp)::term,dj,v,old
      integer::j
      sumv=0.0_dp;term=1.0_dp
      do j=1,500
         old=sumv;sumv=sumv+term
         if(abs(sumv-old)<=epsilon(sumv)*max(1.0_dp,abs(sumv)))exit
         dj=real(j,dp);v=real(2*j-1,dp)
         term=term*(0.25_dp*v*v/(c+dj-1.0_dp))*(x/dj)
         if(abs(term)<1e-16_dp*abs(sumv))exit
      end do
   end function hyper2f1_half

   pure real(dp) function dpearson(r,n,rho) result(f)
      real(dp),intent(in)::r,rho;integer,intent(in)::n
      real(dp)::dn,scale,lf,c,x
      if(n<3 .or. abs(r)>=1.0_dp .or. abs(rho)>1.0_dp)then;f=0.0_dp;return;end if
      dn=real(n,dp);scale=(dn-2.0_dp)/(sqrt2*(dn-1.0_dp))
      lf=0.5_dp*(dn-1.0_dp)*log(max(1e-300_dp,1.0_dp-rho*rho)) + &
         0.5_dp*(dn-4.0_dp)*log(max(1e-300_dp,1.0_dp-r*r)) + &
         (1.5_dp-dn)*log(max(1e-300_dp,1.0_dp-rho*r)) + &
         log_gamma(dn)-log_gamma(dn-0.5_dp)-0.5_dp*log(acos(-1.0_dp))
      c=dn-0.5_dp;x=0.5_dp*(1.0_dp+rho*r)
      f=scale*exp(lf)*hyper2f1_half(c,x)
   end function dpearson

   real(dp) function ppearson(r,n,rho) result(p)
      real(dp),intent(in)::r,rho;integer,intent(in)::n
      if(r<=-1.0_dp)then;p=0.0_dp;return;end if
      if(r>=1.0_dp)then;p=1.0_dp;return;end if
      p=adaptive_integral(fun,-1.0_dp,r,2e-9_dp)
      p=max(0.0_dp,min(1.0_dp,p))
   contains
      function fun(x) result(y)
         real(dp),intent(in)::x;real(dp)::y
         y=dpearson(x,n,rho)
      end function fun
   end function ppearson

   real(dp) function qpearson(p,n,rho) result(r)
      real(dp),intent(in)::p,rho;integer,intent(in)::n
      real(dp)::lo,hi,mid;integer::i
      if(p<=0.0_dp)then;r=-1.0_dp;return;end if
      if(p>=1.0_dp)then;r=1.0_dp;return;end if
      lo=-1.0_dp;hi=1.0_dp
      do i=1,70
         mid=0.5_dp*(lo+hi)
         if(ppearson(mid,n,rho)<p)then;lo=mid;else;hi=mid;end if
      end do
      r=0.5_dp*(lo+hi)
   end function qpearson

   real(dp) function rpearson(n,rho) result(r)
      integer,intent(in)::n;real(dp),intent(in)::rho
      real(dp),allocatable::x(:),y(:);real(dp)::mx,my,sxx,syy,sxy
      integer::i
      allocate(x(n),y(n))
      do i=1,n;x(i)=randn();y(i)=rho*x(i)+sqrt(max(0.0_dp,1.0_dp-rho*rho))*randn();end do
      mx=sum(x)/real(n,dp);my=sum(y)/real(n,dp)
      sxx=sum((x-mx)**2);syy=sum((y-my)**2);sxy=sum((x-mx)*(y-my))
      r=sxy/sqrt(sxx*syy)
   end function rpearson

   function spearson(n,rho) result(s)
      integer,intent(in)::n;real(dp),intent(in)::rho;type(dist_stats)::s
      real(dp)::m,m2,r2,r4,om,om2,om3,om4
      m=1.0_dp/real(n+6,dp);m2=m*m;r2=rho*rho;r4=r2*r2;om=1-r2
      om2=om*om;om3=om2*om;om4=om2*om2
      s%mean=rho-0.5_dp*m*rho*om*(1.0_dp+2.25_dp*m*(3.0_dp+r2)+ &
         0.375_dp*m2*(121.0_dp+70.0_dp*r2+25.0_dp*r4))
      s%median=qpearson(0.5_dp,n,rho)
      s%mode=rho
      s%third_central=-m2*rho*om3*(6.0_dp+m*(69.0_dp+88.0_dp*r2)+ &
         0.75_dp*m2*(797.0_dp+1691.0_dp*r2+1560.0_dp*r4))
      s%fourth_central=3.0_dp*m2*om4*(1.0_dp+m*(12.0_dp+35.0_dp*r2)+ &
         0.25_dp*m2*(436.0_dp+2028.0_dp*r2+3025.0_dp*r4))
      s%variance=m*om2*(1.0_dp+0.5_dp*m*(14.0_dp+11.0_dp*r2)+ &
         0.5_dp*m2*(98.0_dp+130.0_dp*r2+75.0_dp*r4))
   end function spearson
end module suppdists_pearson
