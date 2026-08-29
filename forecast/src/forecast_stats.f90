module forecast_stats
   use forecast_kinds, only : dp, pi
   implicit none
   private
   public :: mean_value, variance_value, median_value, quantile_value, quantile_type8, acf_values, pacf_values
   public :: normal_cdf, normal_quantile, rmse_value, mae_value, autocorrelation
contains
   pure function mean_value(x) result(v)
      real(dp),intent(in) :: x(:)
      real(dp) :: v
      if(size(x)==0) then
      v=0.0_dp
      else
      v=sum(x)/real(size(x),dp)
      end if
   end function
   pure function variance_value(x,unbiased) result(v)
      real(dp),intent(in) :: x(:)
      logical,intent(in),optional :: unbiased
      real(dp) :: v,m,den
      logical :: u
      u=.true.
      if(present(unbiased)) u=unbiased
      if(size(x)<2) then
      v=0.0_dp
      return
      end if
      m=mean_value(x)
      den=real(size(x)-merge(1,0,u),dp)
      v=sum((x-m)**2)/den
   end function
   function median_value(x) result(v)
      real(dp),intent(in) :: x(:)
      real(dp) :: v
      v=quantile_value(x,0.5_dp)
   end function
   function quantile_value(x,p) result(v)
      real(dp),intent(in) :: x(:),p
      real(dp) :: v,h,g
      real(dp),allocatable :: z(:)
      integer :: n,j,i,k
      n=size(x)
      if(n==0) then
      v=0.0_dp
      return
      end if
      allocate(z(n))
      z=x
      do i=2,n
         v=z(i)
         j=i-1
         do while(j>=1)
            if(z(j)<=v) exit
            z(j+1)=z(j)
            j=j-1
         end do
         z(j+1)=v
      end do
      h=1.0_dp+(real(n-1,dp))*max(0.0_dp,min(1.0_dp,p))
      j=int(floor(h))
      g=h-real(j,dp)
      k=min(n,j+1)
      v=(1.0_dp-g)*z(max(1,j))+g*z(k)
   end function
   function quantile_type8(x,p) result(v)
      real(dp),intent(in)::x(:),p
      real(dp)::v,h,g,tmp,pp
      real(dp),allocatable::z(:)
      integer::n,j,i,k
      n=size(x)
      if(n==0)then
         v=0.0_dp
         return
      end if
      allocate(z(n))
      z=x
      do i=2,n
         tmp=z(i)
         j=i-1
         do while(j>=1)
            if(z(j)<=tmp)exit
            z(j+1)=z(j)
            j=j-1
         end do
         z(j+1)=tmp
      end do
      pp=max(0.0_dp,min(1.0_dp,p))
      if(pp<=0.0_dp)then
         v=z(1)
         return
      else if(pp>=1.0_dp)then
         v=z(n)
         return
      end if
      h=(real(n,dp)+1.0_dp/3.0_dp)*pp+1.0_dp/3.0_dp
      j=int(floor(h))
      g=h-real(j,dp)
      if(j<=0)then
         v=z(1)
      else if(j>=n)then
         v=z(n)
      else
         k=j+1
         v=(1.0_dp-g)*z(j)+g*z(k)
      end if
   end function quantile_type8

   function acf_values(x,maxlag,correlation) result(a)
      real(dp),intent(in) :: x(:)
      integer,intent(in) :: maxlag
      logical,intent(in),optional :: correlation
      real(dp),allocatable :: a(:)
      real(dp) :: m,c0
      integer :: k,n
      logical :: cor
      n=size(x)
      cor=.true.
      if(present(correlation)) cor=correlation
      allocate(a(0:maxlag))
      m=mean_value(x)
      c0=sum((x-m)**2)/real(n,dp)
      do k=0,maxlag
         if(k>=n) then
         a(k)=0.0_dp
         else
         a(k)=sum((x(1:n-k)-m)*(x(1+k:n)-m))/real(n,dp)
         end if
         if(cor .and. c0>0.0_dp) a(k)=a(k)/c0
      end do
   end function
   function pacf_values(x,maxlag) result(pacf)
      real(dp),intent(in) :: x(:)
      integer,intent(in) :: maxlag
      real(dp),allocatable :: pacf(:),a(:),phi(:,:),v(:)
      integer :: k,j
      a=acf_values(x,maxlag,.true.)
      allocate(pacf(0:maxlag),phi(maxlag,maxlag),v(maxlag))
      pacf=0.0_dp
      phi=0.0_dp
      pacf(0)=1.0_dp
      if(maxlag==0) return
      phi(1,1)=a(1)
      pacf(1)=phi(1,1)
      v(1)=1.0_dp-phi(1,1)**2
      do k=2,maxlag
         if(abs(v(k-1))<tiny(1.0_dp)) then
            phi(k,k)=0.0_dp
         else
            phi(k,k)=(a(k)-sum(phi(k-1,1:k-1)*a(k-1:1:-1)))/v(k-1)
         end if
         do j=1,k-1
         phi(k,j)=phi(k-1,j)-phi(k,k)*phi(k-1,k-j)
         end do
         v(k)=v(k-1)*(1.0_dp-phi(k,k)**2)
         pacf(k)=phi(k,k)
      end do
   end function
   pure function autocorrelation(x,lag) result(r)
      real(dp),intent(in) :: x(:)
      integer,intent(in) :: lag
      real(dp) :: r,m,d
      integer :: n
      n=size(x)
      if(lag<0 .or. lag>=n) then
      r=0.0_dp
      return
      end if
      m=mean_value(x)
      d=sum((x-m)**2)
      if(d<=0.0_dp) then
      r=0.0_dp
      else
      r=sum((x(1:n-lag)-m)*(x(1+lag:n)-m))/d
      end if
   end function
   pure function rmse_value(e) result(v)
      real(dp),intent(in) :: e(:)
      real(dp)::v
      if(size(e)==0) then
      v=0.0_dp
      else
      v=sqrt(sum(e*e)/real(size(e),dp))
      end if
   end function
   pure function mae_value(e) result(v)
      real(dp),intent(in) :: e(:)
      real(dp)::v
      if(size(e)==0) then
      v=0.0_dp
      else
      v=sum(abs(e))/real(size(e),dp)
      end if
   end function
   pure function normal_cdf(x) result(p)
      real(dp),intent(in)::x
      real(dp)::p
      p=0.5_dp*erfc(-x/sqrt(2.0_dp))
   end function
   pure function normal_quantile(p) result(x)
      real(dp),intent(in) :: p
      real(dp) :: x,q,r
      real(dp),parameter :: a1=-3.969683028665376e1_dp,a2=2.209460984245205e2_dp
      real(dp),parameter :: a3=-2.759285104469687e2_dp,a4=1.383577518672690e2_dp
      real(dp),parameter :: a5=-3.066479806614716e1_dp,a6=2.506628277459239_dp
      real(dp),parameter :: b1=-5.447609879822406e1_dp,b2=1.615858368580409e2_dp
      real(dp),parameter :: b3=-1.556989798598866e2_dp,b4=6.680131188771972e1_dp,b5=-1.328068155288572e1_dp
      real(dp),parameter :: c1=-7.784894002430293e-3_dp,c2=-3.223964580411365e-1_dp
      real(dp),parameter :: c3=-2.400758277161838_dp,c4=-2.549732539343734_dp
      real(dp),parameter :: c5=4.374664141464968_dp,c6=2.938163982698783_dp
      real(dp),parameter :: d1=7.784695709041462e-3_dp,d2=3.224671290700398e-1_dp
      real(dp),parameter :: d3=2.445134137142996_dp,d4=3.754408661907416_dp
      if(p<=0.0_dp) then
      x=-huge(1.0_dp)
      return
      end if
      if(p>=1.0_dp) then
      x=huge(1.0_dp)
      return
      end if
      if(p<0.02425_dp) then
         q=sqrt(-2.0_dp*log(p))
         x=(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
      else if(p>0.97575_dp) then
         q=sqrt(-2.0_dp*log(1.0_dp-p))
         x=-(((((c1*q+c2)*q+c3)*q+c4)*q+c5)*q+c6)/((((d1*q+d2)*q+d3)*q+d4)*q+1.0_dp)
      else
         q=p-0.5_dp
         r=q*q
         x=((((a1*r+a2)*r+a3)*r+a4)*r+a5)*r+a6
         x=x*q/(((((b1*r+b2)*r+b3)*r+b4)*r+b5)*r+1.0_dp)
      end if
   end function
end module forecast_stats
