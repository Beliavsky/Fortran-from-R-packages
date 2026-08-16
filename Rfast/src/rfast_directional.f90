module rfast_directional
   use rfast_special, only : dp, pi, f_cdf
   use rfast_arrays, only : mean_r, sort_real, median_r
   use rfast_distributions, only : vonmises_random
   implicit none
   private
   type, public :: circular_fit
      real(dp) :: mu = 0.0_dp
      real(dp) :: kappa = 0.0_dp
      real(dp) :: loglik = -huge(1.0_dp)
      integer :: iterations = 0
   end type circular_fit
   type, public :: vmf_fit
      real(dp), allocatable :: mu(:)
      real(dp) :: kappa = 0.0_dp
      real(dp) :: loglik = -huge(1.0_dp)
   end type vmf_fit
   public :: vm_mle, vmf_mle, rvonmises, watson_test, kuiper_test, circlin_cor
   public :: circular_mean, circular_resultant

contains

   pure real(dp) function bessi0(x) result(ans)
      real(dp), intent(in) :: x
      real(dp) :: ax,y
      ax=abs(x)
      if(ax<3.75_dp)then
         y=(x/3.75_dp)**2
         ans=1.0_dp+y*(3.5156229_dp+y*(3.0899424_dp+y*(1.2067492_dp+ &
             y*(0.2659732_dp+y*(0.0360768_dp+y*0.0045813_dp)))))
      else
         y=3.75_dp/ax
         ans=(exp(ax)/sqrt(ax))*(0.39894228_dp+y*(0.01328592_dp+y*(0.00225319_dp+ &
             y*(-0.00157565_dp+y*(0.00916281_dp+y*(-0.02057706_dp+y*(0.02635537_dp+ &
             y*(-0.01647633_dp+y*0.00392377_dp))))))))
      end if
   end function bessi0

   pure real(dp) function bessi1(x) result(ans)
      real(dp), intent(in) :: x
      real(dp) :: ax,y
      ax=abs(x)
      if(ax<3.75_dp)then
         y=(x/3.75_dp)**2
         ans=ax*(0.5_dp+y*(0.87890594_dp+y*(0.51498869_dp+y*(0.15084934_dp+ &
             y*(0.02658733_dp+y*(0.00301532_dp+y*0.00032411_dp))))))
      else
         y=3.75_dp/ax
         ans=0.39894228_dp+y*(-0.03988024_dp+y*(-0.00362018_dp+y*(0.00163801_dp+ &
             y*(-0.01031555_dp+y*(0.02282967_dp+y*(-0.02895312_dp+y*(0.01787654_dp- &
             y*0.00420059_dp)))))))
         ans=ans*exp(ax)/sqrt(ax)
      end if
      if(x<0)ans=-ans
   end function bessi1

   pure real(dp) function circular_mean(x) result(mu)
      real(dp), intent(in) :: x(:)
      mu=atan2(sum(sin(x)),sum(cos(x)))
      if(mu<0)mu=mu+2*pi
   end function circular_mean

   pure real(dp) function circular_resultant(x) result(r)
      real(dp), intent(in) :: x(:)
      real(dp)::c,s
      c=sum(cos(x))/size(x);s=sum(sin(x))/size(x);r=sqrt(c*c+s*s)
   end function circular_resultant

   function vm_mle(x,tol,maxiter) result(fit)
      real(dp), intent(in) :: x(:)
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: maxiter
      type(circular_fit) :: fit
      real(dp)::r,k,k2,a1,der,eps,con
      integer::it,mi,n
      n=size(x);fit%mu=circular_mean(x);r=circular_resultant(x);con=sum(cos(x-fit%mu))
      k=max(1e-10_dp,(1.28_dp-0.53_dp*r*r)*tan(0.5_dp*pi*min(r,0.999999_dp)));k2=k
      eps=1e-9_dp;if(present(tol))eps=tol;mi=100;if(present(maxiter))mi=maxiter
      do it=1,mi
         a1=bessi1(k)/bessi0(k);der=1.0_dp-a1*a1-a1/k;k2=max(0.0_dp,k-(a1-r)/der)
         if(abs(k2-k)<eps)exit;k=k2
      end do
      fit%kappa=k2;fit%iterations=it;fit%loglik=k2*con-real(n,dp)*log(2*pi*bessi0(k2))
   end function vm_mle

   function vmf_mle(x) result(fit)
      real(dp), intent(in) :: x(:,:)
      type(vmf_fit) :: fit
      real(dp)::m(size(x,2)),r,p
      integer::n
      n=size(x,1);p=real(size(x,2),dp);m=sum(x,dim=1);r=sqrt(sum(m*m))/real(n,dp)
      allocate(fit%mu(size(x,2)))
      if(r>0)then;fit%mu=m/(real(n,dp)*r);else;fit%mu=0;end if
      fit%kappa=max(0.0_dp,r*(p-r*r)/max(1e-12_dp,1-r*r))
      ! Log-likelihood normalization for arbitrary p requires I_(p/2-1); omitted from the approximation.
      fit%loglik=huge(1.0_dp)*(-1.0_dp)
   end function vmf_mle

   function rvonmises(n,mu,kappa) result(x)
      integer,intent(in)::n;real(dp),intent(in)::mu,kappa;real(dp)::x(n);integer::i
      do i=1,n;x(i)=vonmises_random(mu,kappa);end do
   end function rvonmises

   function watson_test(theta) result(out)
      real(dp),intent(in)::theta(:);real(dp)::out(2),u(size(theta)),avg,w,pv;integer::i,m,n
      n=size(theta);u=modulo(theta,2*pi)/(2*pi);call sort_real(u);avg=sum(u)/n;w=0.0_dp
      do i=1,n;w=w+(u(i)-real(i,dp)/n+0.5_dp/n-(avg-0.5_dp))**2;end do
      w=w+1.0_dp/(12.0_dp*n);pv=0.0_dp
      do m=1,20;pv=pv+2.0_dp*(-1.0_dp)**(m-1)*exp(-2.0_dp*m*m*pi*pi*w);end do
      out=[w,max(0.0_dp,min(1.0_dp,pv))]
   end function watson_test

   function kuiper_test(theta) result(out)
      real(dp),intent(in)::theta(:);real(dp)::out(2),u(size(theta)),dpv,dmv,v,pv,lam;integer::i,n,j
      n=size(theta);u=modulo(theta,2*pi)/(2*pi);call sort_real(u);dpv=0.0_dp;dmv=0.0_dp
      do i=1,n
         dpv=max(dpv,real(i,dp)/n-u(i));dmv=max(dmv,u(i)-real(i-1,dp)/n)
      end do
      v=dpv+dmv;lam=(sqrt(real(n,dp))+0.155_dp+0.24_dp/sqrt(real(n,dp)))*v;pv=0.0_dp
      do j=1,50
         pv=pv+2.0_dp*(4.0_dp*j*j*lam*lam-1.0_dp)*exp(-2.0_dp*j*j*lam*lam)
      end do
      out=[v,max(0.0_dp,min(1.0_dp,pv))]
   end function kuiper_test

   function circlin_cor(theta,x) result(out)
      real(dp),intent(in)::theta(:),x(:);real(dp)::out(2),c(size(theta)),s(size(theta)),rxc,rxs,rcs,r2,f
      c=cos(theta);s=sin(theta);rxc=cor_vec(c,x);rxs=cor_vec(s,x);rcs=cor_vec(c,s)
      r2=(rxc*rxc+rxs*rxs-2*rxc*rxs*rcs)/max(1e-15_dp,1-rcs*rcs);r2=max(0.0_dp,min(0.999999999_dp,r2))
      f=real(size(theta)-3,dp)*r2/(1-r2);out=[r2,1.0_dp-f_cdf(f,2.0_dp,real(size(theta)-3,dp))]
   end function circlin_cor

   pure real(dp) function cor_vec(x,y) result(r)
      real(dp),intent(in)::x(:),y(:);real(dp)::mx,my,sx,sy
      mx=sum(x)/size(x);my=sum(y)/size(y);sx=sqrt(sum((x-mx)**2));sy=sqrt(sum((y-my)**2))
      if(sx>0.and.sy>0)then;r=sum((x-mx)*(y-my))/(sx*sy);else;r=0;end if
   end function cor_vec

end module rfast_directional
