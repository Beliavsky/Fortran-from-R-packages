module suppdists_inverse_gaussian
   use suppdists_kinds, only : dp, sqrt2pi
   use suppdists_special, only : normal_cdf, randn
   use suppdists_stats, only : dist_stats
   implicit none
   private
   public :: dinvgauss, pinvgauss, qinvgauss, rinvgauss, sinvgauss
contains
   pure real(dp) function dinvgauss(x,nu,lambda) result(f)
      real(dp), intent(in) :: x,nu,lambda
      real(dp) :: delta,ratio
      if(x<=0.0_dp .or. nu<=0.0_dp .or. lambda<=0.0_dp)then
         f=0.0_dp; return
      end if
      delta=x/nu-1.0_dp; ratio=lambda/x
      f=sqrt(ratio/(2.0_dp*acos(-1.0_dp)*x*x))*exp(-0.5_dp*ratio*delta*delta)
   end function dinvgauss

   pure real(dp) function pinvgauss(x,nu,lambda) result(p)
      real(dp), intent(in) :: x,nu,lambda
      real(dp) :: a,b,z1,z2,c,lp2
      if(x<=0.0_dp)then;p=0.0_dp;return;end if
      if(nu<=0.0_dp .or. lambda<=0.0_dp)then;p=0.0_dp;return;end if
      a=sqrt(lambda/x); b=x/nu
      z1=a*(b-1.0_dp); z2=-a*(b+1.0_dp)
      p=normal_cdf(z1)
      c=2.0_dp*lambda/nu
      if(z2 > -38.0_dp .and. c < 700.0_dp) then
         lp2=log(max(normal_cdf(z2),tiny(1.0_dp)))+c
         if(lp2 < 0.0_dp) p=p+exp(lp2)
      end if
      p=max(0.0_dp,min(1.0_dp,p))
   end function pinvgauss

   pure real(dp) function qinvgauss(p,nu,lambda) result(x)
      real(dp), intent(in) :: p,nu,lambda
      real(dp) :: lo,hi,mid
      integer :: i
      if(p<=0.0_dp)then;x=0.0_dp;return;end if
      if(p>=1.0_dp)then;x=huge(1.0_dp);return;end if
      lo=0.0_dp; hi=max(nu,1.0_dp)
      do while(pinvgauss(hi,nu,lambda)<p .and. hi<huge(1.0_dp)/4.0_dp)
         hi=2.0_dp*hi
      end do
      do i=1,120
         mid=0.5_dp*(lo+hi)
         if(pinvgauss(mid,nu,lambda)<p)then;lo=mid;else;hi=mid;end if
      end do
      x=0.5_dp*(lo+hi)
   end function qinvgauss

   real(dp) function rinvgauss(nu,lambda) result(x)
      real(dp), intent(in) :: nu,lambda
      real(dp) :: y,w,u
      y=randn()**2
      w=nu + (nu*nu*y)/(2.0_dp*lambda) - &
        (nu/(2.0_dp*lambda))*sqrt(4.0_dp*nu*lambda*y+nu*nu*y*y)
      call random_number(u)
      if(u<=nu/(nu+w))then;x=w;else;x=nu*nu/w;end if
   end function rinvgauss

   pure function sinvgauss(nu,lambda) result(s)
      real(dp), intent(in) :: nu,lambda
      type(dist_stats) :: s
      real(dp) :: v
      s%mean=nu
      v=nu**3/lambda; s%variance=v
      s%median=qinvgauss_pure(0.5_dp,nu,lambda)
      s%mode=nu*(sqrt(1.0_dp+9.0_dp*nu*nu/(4.0_dp*lambda*lambda))- &
         3.0_dp*nu/(2.0_dp*lambda))
      s%third_central=3.0_dp*nu**5/lambda**2
      s%fourth_central=(3.0_dp+15.0_dp*nu/lambda)*v*v
   contains
      pure real(dp) function qinvgauss_pure(p,mu,lam) result(xx)
         real(dp), intent(in) :: p,mu,lam
         real(dp) :: lo,hi,mid
         integer :: j
         lo=0.0_dp; hi=max(mu,1.0_dp)
         do j=1,100
            if(pinvgauss(hi,mu,lam)>=p)exit
            hi=2.0_dp*hi
         end do
         do j=1,100
            mid=0.5_dp*(lo+hi)
            if(pinvgauss(mid,mu,lam)<p)then;lo=mid;else;hi=mid;end if
         end do
         xx=0.5_dp*(lo+hi)
      end function qinvgauss_pure
   end function sinvgauss
end module suppdists_inverse_gaussian
