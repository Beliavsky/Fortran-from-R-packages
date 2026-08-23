module rfast_mle
   use rfast_special, only : dp, pi, digamma_r, trigamma_r, log_beta, log1p_r, expm1_r
   use rfast_arrays, only : mean_r, median_r, variance_r, nth_value
   implicit none
   private
   type, public :: mle_result
      real(dp), allocatable :: param(:)
      real(dp) :: loglik = -huge(1.0_dp)
      integer :: iters = 0
      integer :: status = 0
   end type mle_result
   public :: normal_mle, lognormal_mle, exponential_mle, exponential2_mle
   public :: poisson_mle, geometric_mle, ztp_mle, zip_mle, negbin_mle
   public :: gamma_mle, beta_mle, laplace_mle, cauchy_mle, logistic_mle
   public :: pareto_mle, rayleigh_mle, invgauss_mle, lindley_mle, maxboltz_mle
   public :: chisq_mle, wigner_mle, gumbel_mle

contains

   function normal_mle(x) result(res)
      real(dp),intent(in)::x(:);type(mle_result)::res;real(dp)::m,v;integer::n
      n=size(x);m=mean_r(x);v=sum((x-m)**2)/real(n,dp);allocate(res%param(2));res%param=[m,v]
      res%loglik=-0.5_dp*real(n,dp)*(log(2.0_dp*pi*v)+1.0_dp);res%iters=1
   end function normal_mle

   function lognormal_mle(x) result(res)
      real(dp),intent(in)::x(:);type(mle_result)::res;real(dp),allocatable::y(:);real(dp)::m,v;integer::n
      n=size(x);if(any(x<=0.0_dp))then;res%status=1;return;end if;y=log(x);m=mean_r(y);v=sum((y-m)**2)/real(n,dp)
      allocate(res%param(2));res%param=[m,v];res%loglik=-0.5_dp*real(n,dp)*(log(2*pi*v)+1.0_dp)-sum(y);res%iters=1
   end function lognormal_mle

   function exponential_mle(x) result(res)
      real(dp),intent(in)::x(:);type(mle_result)::res;real(dp)::scale;integer::n
      n=size(x);scale=mean_r(x);allocate(res%param(1));res%param=[scale];res%loglik=-real(n,dp)*(log(scale)+1.0_dp);res%iters=1
   end function exponential_mle

   function exponential2_mle(x) result(res)
      real(dp),intent(in)::x(:);type(mle_result)::res;real(dp)::a,b;integer::n
      n=size(x);a=minval(x);b=mean_r(x)-a;allocate(res%param(2));res%param=[a,b]
      res%loglik=-real(n,dp)*log(b)-sum((x-a)/b);res%iters=1
   end function exponential2_mle

   function poisson_mle(x) result(res)
      integer,intent(in)::x(:);type(mle_result)::res;real(dp)::lam,sx;integer::i,n
      n=size(x);sx=real(sum(x),dp);lam=sx/real(n,dp);allocate(res%param(1));res%param=[lam]
      res%loglik=-real(n,dp)*lam
      if(lam>0.0_dp)res%loglik=res%loglik+sx*log(lam)
      do i=1,n;res%loglik=res%loglik-log_gamma(real(x(i)+1,dp));end do;res%iters=1
   end function poisson_mle

   function geometric_mle(x,type) result(res)
      integer,intent(in)::x(:);integer,intent(in),optional::type;type(mle_result)::res;integer::t,n;real(dp)::p,sx
      t=1;if(present(type))t=type;n=size(x);sx=real(sum(x),dp)
      if(t==1)then;p=1.0_dp/(1.0_dp+sx/real(n,dp));res%loglik=real(n,dp)*log(p)+sx*log1p_r(-p)
      else;p=real(n,dp)/sx;res%loglik=real(n,dp)*log(p)+(sx-real(n,dp))*log1p_r(-p);end if
      allocate(res%param(1));res%param=[p];res%iters=1
   end function geometric_mle

   function ztp_mle(x,tol,maxiter) result(res)
      integer,intent(in)::x(:);real(dp),intent(in),optional::tol;integer,intent(in),optional::maxiter
      type(mle_result)::res;real(dp)::lam,lam2,ex,a1,a2,f,fp,eps,sx;integer::i,n,mi,j
      eps=1e-10_dp;if(present(tol))eps=tol;mi=100;if(present(maxiter))mi=maxiter
      n=size(x);sx=real(sum(x),dp);lam=sx/real(n,dp);lam2=lam
      do i=1,mi
         ex=exp(lam);a1=sx/lam;a2=real(n,dp)*ex/(ex-1.0_dp);f=a1-a2;fp=-a1/lam+a2/(ex-1.0_dp);lam2=lam-f/fp
         if(abs(lam2-lam)<eps)exit;lam=max(lam2,1e-12_dp)
      end do
      allocate(res%param(1));res%param=[lam2];res%iters=i
      res%loglik=sx*log(lam2)-real(n,dp)*log(expm1_r(lam2));do j=1,n;res%loglik=res%loglik-log_gamma(real(x(j)+1,dp));end do
   end function ztp_mle

   function zip_mle(x,tol,maxiter) result(res)
      integer,intent(in)::x(:);real(dp),intent(in),optional::tol;integer,intent(in),optional::maxiter
      type(mle_result)::res;integer::n,no,n1,i,j,mi;real(dp)::prop,sx,m,s,l1,l2,f,der,p,eps
      n=size(x);no=count(x==0);n1=n-no;prop=real(no,dp)/n;sx=real(sum(x),dp);m=sx/n
      s=(sum(real(x,dp)**2)-m*sx)/max(1.0_dp,real(n-1,dp));l1=max(1e-6_dp,s/m+m-1.0_dp)
      eps=1e-10_dp;if(present(tol))eps=tol;mi=100;if(present(maxiter))mi=maxiter;l2=l1
      do i=1,mi
         f=m-m*exp(-l1)-l1+prop*l1;der=m*exp(-l1)-1.0_dp+prop;l2=max(1e-12_dp,l1-f/der)
         if(abs(l2-l1)<eps)exit;l1=l2
      end do
      p=max(0.0_dp,min(1.0_dp,1.0_dp-m/l2));allocate(res%param(2));res%param=[l2,p];res%iters=i
      res%loglik=real(no,dp)*log(p+(1-p)*exp(-l2))+real(n1,dp)*log(max(tiny(1.0_dp),1-p))
      do j=1,n;if(x(j)>0)res%loglik=res%loglik-l2+real(x(j),dp)*log(l2)-log_gamma(real(x(j)+1,dp));end do
   end function zip_mle

   function negbin_mle(x,tol,maxiter) result(res)
      integer,intent(in)::x(:);real(dp),intent(in),optional::tol;integer,intent(in),optional::maxiter
      type(mle_result)::res;real(dp)::m,m2,p,r,lr,lr2,f,fp,eps,sx;integer::n,i,j,mi
      n=size(x);sx=real(sum(x),dp);m=sx/n;m2=sum(real(x,dp)**2)/n;p=1.0_dp-m/max(1e-15_dp,m2-m*m);r=abs(m/max(p,1e-8_dp)-m)
      r=max(r,1e-6_dp);lr=log(r);lr2=lr;eps=1e-9_dp;if(present(tol))eps=tol;mi=100;if(present(maxiter))mi=maxiter
      do i=1,mi
         r=exp(lr);f=0.0_dp;fp=0.0_dp
         do j=1,n;f=f+digamma_r(real(x(j),dp)+r);fp=fp+trigamma_r(real(x(j),dp)+r);end do
         f=r*(f-real(n,dp)*digamma_r(r)+real(n,dp)*log(r)-real(n,dp)*log(r+m))
         fp=f+r*r*(fp-real(n,dp)*trigamma_r(r))+real(n,dp)*r-real(n,dp)*r*r/(r+m)
         lr2=lr-f/fp;if(abs(lr2-lr)<eps)exit;lr=lr2
      end do
      r=exp(lr2);p=r/(r+m);allocate(res%param(3));res%param=[p,r,m];res%iters=i;res%loglik=0.0_dp
      do j=1,n
         res%loglik=res%loglik+log_gamma(real(x(j),dp)+r)-log_gamma(real(x(j)+1,dp))-log_gamma(r) &
                    +real(x(j),dp)*log1p_r(-p)+r*log(p)
      end do
   end function negbin_mle

   function gamma_mle(x,tol,maxiter) result(res)
      real(dp),intent(in)::x(:);real(dp),intent(in),optional::tol;integer,intent(in),optional::maxiter
      type(mle_result)::res;real(dp)::m,slx,s,a,a2,b,eps;integer::n,i,mi
      n=size(x);if(any(x<=0))then;res%status=1;return;end if;m=mean_r(x);slx=sum(log(x))/n;s=log(m)-slx
      if(s<=epsilon(1.0_dp))then;a=huge(1.0_dp)**0.25_dp;else;a=(3-s+sqrt((s-3)**2+24*s))/(12*s);end if
      eps=1e-10_dp;if(present(tol))eps=tol;mi=100;if(present(maxiter))mi=maxiter;a2=a
      do i=1,mi
         a2=a-(log(a)-digamma_r(a)-s)/(1.0_dp/a-trigamma_r(a));if(abs(a2-a)<eps)exit;a=max(a2,1e-10_dp)
      end do
      a=a2;b=a/m;allocate(res%param(2));res%param=[a,b];res%iters=i
      res%loglik=-b*real(n,dp)*m+(a-1)*real(n,dp)*slx+real(n,dp)*a*log(b)-real(n,dp)*log_gamma(a)
   end function gamma_mle

   function beta_mle(x,tol,maxiter) result(res)
      real(dp),intent(in)::x(:);real(dp),intent(in),optional::tol;integer,intent(in),optional::maxiter
      type(mle_result)::res;real(dp)::l1,l2,m,v,phi,a,b,da,db,dab,daa,dbb,det,na,nb,eps;integer::n,i,mi
      n=size(x);if(any(x<=0).or.any(x>=1))then;res%status=1;return;end if;l1=sum(log(x))/n;l2=sum(log1p_r(-x))/n
      m=mean_r(x);v=sum((x-m)**2)/n;phi=max(1e-3_dp,m*(1-m)/max(v,1e-15_dp)-1);a=max(1e-3_dp,m*phi);b=max(1e-3_dp,(1-m)*phi)
      eps=1e-10_dp;if(present(tol))eps=tol;mi=100;if(present(maxiter))mi=maxiter;na=a;nb=b
      do i=1,mi
         phi=a+b;da=l1-digamma_r(a)+digamma_r(phi);db=l2-digamma_r(b)+digamma_r(phi);dab=trigamma_r(phi)
         daa=-trigamma_r(a)+dab;dbb=-trigamma_r(b)+dab;det=daa*dbb-dab*dab
         na=a-(dbb*da-dab*db)/det;nb=b-(-dab*da+daa*db)/det
         if(na<=0.or.nb<=0)then;na=0.5_dp*a;nb=0.5_dp*b;end if
         if(abs(na-a)+abs(nb-b)<eps)exit;a=na;b=nb
      end do
      a=na;b=nb;allocate(res%param(2));res%param=[a,b];res%iters=i
      res%loglik=-real(n,dp)*log_beta(a,b)+(a-1)*sum(log(x))+(b-1)*sum(log1p_r(-x))
   end function beta_mle

   function laplace_mle(x) result(res)
      real(dp),intent(in)::x(:);type(mle_result)::res;real(dp)::m,b;integer::n
      n=size(x);m=median_r(x);b=sum(abs(x-m))/n;allocate(res%param(2));res%param=[m,b];res%loglik=-n*log(2*b)-n;res%iters=1
   end function laplace_mle

   function pareto_mle(x) result(res)
      real(dp),intent(in)::x(:);type(mle_result)::res;real(dp)::xm,a,sl;integer::n
      n=size(x);xm=minval(x);sl=sum(log(x));a=n/(sl-n*log(xm));allocate(res%param(2));res%param=[xm,a]
      res%loglik=n*log(a)+a*n*log(xm)-(a+1)*sl;res%iters=1
   end function pareto_mle

   function rayleigh_mle(x) result(res)
      real(dp),intent(in)::x(:);type(mle_result)::res;real(dp)::s;integer::n
      n=size(x);s=0.5_dp*sum(x*x)/n;allocate(res%param(1));res%param=[s];res%loglik=sum(log(x))-n*log(s)-n;res%iters=1
   end function rayleigh_mle

   function invgauss_mle(x) result(res)
      real(dp),intent(in)::x(:);type(mle_result)::res;real(dp)::m,lam,sx,sxi;integer::n
      n=size(x);sx=sum(x);sxi=sum(1.0_dp/x);m=sx/n;lam=1.0_dp/(sxi/n-1.0_dp/m);allocate(res%param(2));res%param=[m,lam]
      res%loglik=0.5_dp*n*log(lam/(2*pi))-1.5_dp*sum(log(x))-lam/(2*m*m)*(-sx+m*m*sxi);res%iters=1
   end function invgauss_mle

   function lindley_mle(x) result(res)
      real(dp),intent(in)::x(:);type(mle_result)::res;real(dp)::a,b,theta,sx;integer::n
      n=size(x);sx=sum(x);a=sx/n;b=a-1;theta=0.5_dp*(-b+sqrt(b*b+8*a))/a;allocate(res%param(1));res%param=[theta]
      res%loglik=2*n*log(theta)-n*log1p_r(theta)+sum(log1p_r(x))-theta*sx;res%iters=1
   end function lindley_mle

   function maxboltz_mle(x) result(res)
      real(dp),intent(in)::x(:);type(mle_result)::res;real(dp)::a;integer::n
      n=size(x);a=sqrt(sum(x*x)/(3*n));allocate(res%param(1));res%param=[a]
      res%loglik=0.5_dp*n*log(2/pi)+2*sum(log(x))-1.5_dp*n-3*n*log(a);res%iters=1
   end function maxboltz_mle

   function chisq_mle(x,tol,maxiter) result(res)
      real(dp),intent(in)::x(:);real(dp),intent(in),optional::tol;integer,intent(in),optional::maxiter
      type(mle_result)::res;real(dp)::v,v2,der,der2,sl,eps;integer::n,i,mi
      n=size(x);sl=0.5_dp*sum(log(x));v=mean_r(x);v2=v;eps=1e-10_dp;if(present(tol))eps=tol;mi=100;if(present(maxiter))mi=maxiter
      do i=1,mi
         der=-0.5_dp*n*log(2.0_dp)-0.5_dp*n*digamma_r(0.5_dp*v)+sl;der2=-0.25_dp*n*trigamma_r(0.5_dp*v)
         v2=max(1e-10_dp,v-der/der2);if(abs(v2-v)<eps)exit;v=v2
      end do
      allocate(res%param(1));res%param=[v2];res%iters=i
      res%loglik=-0.5_dp*n*v2*log(2.0_dp)-n*log_gamma(0.5_dp*v2)+(0.5_dp*v2-1)*sum(log(x))-0.5_dp*sum(x)
   end function chisq_mle

   function wigner_mle(x) result(res)
      real(dp),intent(in)::x(:);type(mle_result)::res;real(dp)::r2;integer::n,i
      n=size(x);r2=maxval(abs(x))**2;allocate(res%param(1));res%param=[r2];res%loglik=n*log(2.0_dp/(pi*r2))
      do i=1,n;if(r2>x(i)*x(i))res%loglik=res%loglik+0.5_dp*log(r2-x(i)*x(i));end do;res%iters=1
   end function wigner_mle

   function cauchy_mle(x,tol,maxiter) result(res)
      real(dp),intent(in)::x(:);real(dp),intent(in),optional::tol;integer,intent(in),optional::maxiter
      type(mle_result) :: res
      integer :: n, i, mi
      real(dp) :: m, s, logs, dm, ds, dmm, dss, dms, det, eps, mn, ln
      real(dp) :: y(size(x)), y2(size(x)), down(size(x)), down2(size(x))
      n=size(x);m=median_r(x);s=max(1e-8_dp,0.5_dp*(nth_value(x,max(1,3*n/4))-nth_value(x,max(1,n/4))));logs=log(s);mn=m;ln=logs
      eps=1e-9_dp;if(present(tol))eps=tol;mi=100;if(present(maxiter))mi=maxiter
      do i=1,mi
         s=exp(logs);y=x-m;y2=y*y;down=1.0_dp/(s*s+y2);down2=down*down
         dm=2*sum(y*down);ds=n-2*s*s*sum(down);dmm=2*sum((y2-s*s)*down2)
         dss=-2*s*s*(dmm+2*s*s*sum(down2));dms=-4*s*s*sum(y*down2);det=dmm*dss-dms*dms
         mn=m-(dss*dm-dms*ds)/det;ln=logs-(-dms*dm+dmm*ds)/det
         if(abs(mn-m)+abs(ln-logs)<eps)exit;m=mn;logs=ln
      end do
      m=mn;s=exp(ln);allocate(res%param(2));res%param=[m,s];res%iters=i
      res%loglik=n*log(s/pi)-sum(log((x-m)**2+s*s))
   end function cauchy_mle

   function logistic_mle(x,tol,maxiter) result(res)
      real(dp),intent(in)::x(:);real(dp),intent(in),optional::tol;integer,intent(in),optional::maxiter
      type(mle_result)::res;integer::n,i,mi;real(dp)::m,s,lgs,dm,ds,dmm,dss,dms,det,m2,l2,eps
      real(dp)::y(size(x)),h(size(x)),sech2(size(x))
      n=size(x);m=mean_r(x);s=sqrt(3.0_dp*variance_r(x))/pi;lgs=log(s);m2=m;l2=lgs;eps=1e-8_dp;if(present(tol))eps=tol
      mi=100;if(present(maxiter))mi=maxiter
      do i=1,mi
         s=exp(lgs);y=(x-m)/s;h=tanh(0.5_dp*y);sech2=1.0_dp/cosh(0.5_dp*y)**2
         dm=sum(h)/s;dmm=-0.5_dp*sum(sech2)/(s*s);ds=-n+sum(y*h);dss=-ds-n-sum((0.5_dp*y)**2*sech2)
         dms=-dm-sum(sech2*0.5_dp*y)/s;det=dmm*dss-dms*dms
         m2=m-(dss*dm-dms*ds)/det;l2=lgs-(-dms*dm+dmm*ds)/det
         if(abs(m2-m)+abs(l2-lgs)<eps)exit;m=m2;lgs=l2
      end do
      m=m2;s=exp(l2);y=(x-m)/s;allocate(res%param(2));res%param=[m,s];res%iters=i
      res%loglik=-n*log(4*s)+sum(log(1.0_dp/cosh(0.5_dp*y)**2))
   end function logistic_mle

   function gumbel_mle(x,tol,maxiter) result(res)
      real(dp),intent(in)::x(:);real(dp),intent(in),optional::tol;integer,intent(in),optional::maxiter
      type(mle_result)::res;integer::n,i,mi;real(dp)::m,s,s2,f,fp,sy,co,eps;real(dp)::y(size(x))
      n = size(x)
      m = mean_r(x)
      s = sqrt(6.0_dp*variance_r(x))/pi
      s2 = s
      eps = 1e-9_dp
      if (present(tol)) eps = tol
      mi = 100
      if (present(maxiter)) mi = maxiter
      do i=1,mi
         y=exp(-(x-m)/s);sy=sum(y);co=sum(x*y);f=s-m+co/sy
         fp=1.0_dp+(sum(x*x*y)*sy-co*co)/(s*s*sy*sy);s2=max(1e-12_dp,s-f/fp)
         if(abs(s2-s)<eps)exit;s=s2
      end do
      s=s2;y=exp(-x/s);m=-s*log(sum(y)/n);allocate(res%param(2));res%param=[m,s];res%iters=i
      y=(x-m)/s;res%loglik=-n*log(s)-sum(y)-sum(exp(-y))
   end function gumbel_mle

end module rfast_mle
