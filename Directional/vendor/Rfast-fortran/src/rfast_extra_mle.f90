module rfast_extra_mle
   use rfast_special, only : dp, pi, digamma_r, trigamma_r, log_beta, log1p_r
   use rfast_arrays, only : mean_r, variance_r
   use rfast_mle, only : mle_result, cauchy_mle, logistic_mle
   implicit none
   private
   public :: weibull_mle, halfnormal_mle, logcauchy_mle, loglogistic_mle
   public :: betaprime_mle, lomax_mle, binomial_mle, borel_mle, logseries_mle

contains

   function weibull_mle(x,tol,maxiter) result(res)
      real(dp),intent(in)::x(:)
      real(dp),intent(in),optional::tol
      integer,intent(in),optional::maxiter
      type(mle_result)::res
      real(dp)::lx(size(x)),k,k2,scale,sw,swl,swl2,g,gp,eps
      integer::n,it,mi
      if(size(x)==0.or.any(x<=0.0_dp))then;res%status=1;return;end if
      n=size(x);lx=log(x);k=max(0.1_dp,1.2_dp/max(1e-6_dp,sqrt(variance_r(lx))));k2=k
      eps=1e-9_dp;if(present(tol))eps=tol;mi=100;if(present(maxiter))mi=maxiter
      do it=1,mi
         sw=sum(exp(k*lx));swl=sum(exp(k*lx)*lx);swl2=sum(exp(k*lx)*lx*lx)
         g=1.0_dp/k+mean_r(lx)-swl/sw
         gp=-1.0_dp/(k*k)-(swl2/sw-(swl/sw)**2)
         k2=k-g/gp
         if(k2<=0.0_dp)k2=0.5_dp*k
         if(abs(k2-k)<eps*max(1.0_dp,k))exit
         k=k2
      end do
      k=k2;sw=sum(exp(k*lx));scale=(sw/real(n,dp))**(1.0_dp/k)
      allocate(res%param(2));res%param=[k,scale];res%iters=it
      res%loglik=real(n,dp)*log(k)-real(n,dp)*k*log(scale)+(k-1.0_dp)*sum(lx)-sum((x/scale)**k)
   end function weibull_mle

   function halfnormal_mle(x) result(res)
      real(dp),intent(in)::x(:);type(mle_result)::res;real(dp)::sigma;integer::n
      if(size(x)==0.or.any(x<0.0_dp))then;res%status=1;return;end if
      n=size(x);sigma=sqrt(sum(x*x)/real(n,dp));allocate(res%param(1));res%param=[sigma]
      res%loglik=real(n,dp)*(0.5_dp*log(2.0_dp/pi)-log(sigma)-0.5_dp);res%iters=1
   end function halfnormal_mle

   function logcauchy_mle(x,tol,maxiter) result(res)
      real(dp),intent(in)::x(:);real(dp),intent(in),optional::tol;integer,intent(in),optional::maxiter
      type(mle_result)::res;real(dp),allocatable::y(:);type(mle_result)::base
      if(size(x)==0.or.any(x<=0.0_dp))then;res%status=1;return;end if
      y=log(x);base=cauchy_mle(y,tol,maxiter);res=base
      if(res%status==0)res%loglik=res%loglik-sum(y)
   end function logcauchy_mle

   function loglogistic_mle(x,tol,maxiter) result(res)
      real(dp),intent(in)::x(:);real(dp),intent(in),optional::tol;integer,intent(in),optional::maxiter
      type(mle_result)::res;real(dp),allocatable::y(:);type(mle_result)::base;real(dp)::a,b;integer::n
      if(size(x)==0.or.any(x<=0.0_dp))then;res%status=1;return;end if
      n=size(x);y=log(x);base=logistic_mle(y,tol,maxiter)
      if(base%status/=0.or..not.allocated(base%param))then;res%status=1;return;end if
      a=exp(base%param(1));b=1.0_dp/base%param(2);allocate(res%param(2));res%param=[a,b];res%iters=base%iters
      res%loglik=real(n,dp)*log(b/a)+(b-1.0_dp)*sum(y)-real(n,dp)*(b-1.0_dp)*log(a) &
                   -2.0_dp*sum(log1p_r((x/a)**b))
   end function loglogistic_mle

   function betaprime_mle(x,tol,maxiter) result(res)
      real(dp),intent(in)::x(:);real(dp),intent(in),optional::tol;integer,intent(in),optional::maxiter
      type(mle_result)::res
      real(dp)::m,v,a,b,na,nb,sa,sb,sab,saa,sbb,det,slx,sl1,eps
      integer::n,it,mi
      if(size(x)==0.or.any(x<=0.0_dp))then;res%status=1;return;end if
      n=size(x);m=mean_r(x);v=sum((x-m)**2)/real(n,dp)
      b=max(2.1_dp,2.0_dp+m*(m+1.0_dp)/max(v,1e-8_dp));a=max(1e-6_dp,m*(b-1.0_dp))
      slx=sum(log(x));sl1=sum(log1p_r(x));na=a;nb=b;eps=1e-9_dp;if(present(tol))eps=tol;mi=100;if(present(maxiter))mi=maxiter
      do it=1,mi
         sa=real(n,dp)*(digamma_r(a+b)-digamma_r(a))+slx-sl1
         sb=real(n,dp)*(digamma_r(a+b)-digamma_r(b))-sl1
         sab=real(n,dp)*trigamma_r(a+b);saa=sab-real(n,dp)*trigamma_r(a);sbb=sab-real(n,dp)*trigamma_r(b)
         det=saa*sbb-sab*sab;if(abs(det)<=tiny(1.0_dp))then;res%status=2;return;end if
         na=a-(sbb*sa-sab*sb)/det;nb=b-(-sab*sa+saa*sb)/det
         if(na<=0.0_dp.or.nb<=0.0_dp)then;na=0.5_dp*a;nb=0.5_dp*b;end if
         if(abs(na-a)+abs(nb-b)<eps)exit;a=na;b=nb
      end do
      a=na;b=nb;allocate(res%param(2));res%param=[a,b];res%iters=it
      res%loglik=(a-1.0_dp)*slx-(a+b)*sl1-real(n,dp)*log_beta(a,b)
   end function betaprime_mle

   function lomax_mle(x,tol,maxiter) result(res)
      real(dp),intent(in)::x(:);real(dp),intent(in),optional::tol;integer,intent(in),optional::maxiter
      type(mle_result)::res
      real(dp)::m,v,a,lam,la,ll,nla,nll,sa,sl,saa,sal,sll,det,com,eps
      integer::n,it,mi
      if(size(x)==0.or.any(x<0.0_dp))then;res%status=1;return;end if
      n=size(x);m=mean_r(x);v=sum((x-m)**2)/real(n,dp)
      if(v>m*m.and.m>0.0_dp)then;a=max(1.01_dp,2.0_dp*v/(v-m*m));else;a=3.0_dp;end if
      lam=max(1e-8_dp,(a-1.0_dp)*max(m,1e-8_dp));la=log(a);ll=log(lam);nla=la;nll=ll
      eps=1e-9_dp;if(present(tol))eps=tol;mi=100;if(present(maxiter))mi=maxiter
      do it=1,mi
         a=exp(la);lam=exp(ll)
         saa=-a*sum(log1p_r(x/lam));sa=real(n,dp)+saa
         com=sum(x/(lam+x));sal=a*com;sl=-real(n,dp)+sal+com
         sll=-(a+1.0_dp)*lam*sum(x/(lam+x)**2);det=saa*sll-sal*sal
         if(abs(det)<=tiny(1.0_dp))then;res%status=2;return;end if
         nla=la-(sll*sa-sal*sl)/det;nll=ll-(-sal*sa+saa*sl)/det
         if(abs(nla-la)+abs(nll-ll)<eps)exit;la=nla;ll=nll
      end do
      a=exp(nla);lam=exp(nll);allocate(res%param(2));res%param=[a,lam];res%iters=it
      res%loglik=real(n,dp)*log(a/lam)-(a+1.0_dp)*sum(log1p_r(x/lam))
   end function lomax_mle

   function binomial_mle(x,ntrials) result(res)
      integer,intent(in)::x(:),ntrials;type(mle_result)::res
      real(dp)::p;integer::i,n
      if(size(x)==0.or.ntrials<=0.or.any(x<0).or.any(x>ntrials))then;res%status=1;return;end if
      n=size(x);p=real(sum(x),dp)/real(n*ntrials,dp);p=min(1.0_dp,max(0.0_dp,p))
      allocate(res%param(2));res%param=[real(ntrials,dp),p];res%iters=1;res%loglik=0.0_dp
      do i=1,n
         res%loglik=res%loglik+log_gamma(real(ntrials+1,dp))-log_gamma(real(x(i)+1,dp)) &
                      -log_gamma(real(ntrials-x(i)+1,dp))
         if(x(i)>0.and.p>0.0_dp)res%loglik=res%loglik+real(x(i),dp)*log(p)
         if(x(i)<ntrials.and.p<1.0_dp)res%loglik=res%loglik+real(ntrials-x(i),dp)*log1p_r(-p)
      end do
   end function binomial_mle

   function borel_mle(x) result(res)
      integer,intent(in)::x(:);type(mle_result)::res;real(dp)::m,sx;integer::i,n
      if(size(x)==0.or.any(x<1))then;res%status=1;return;end if
      n=size(x);sx=real(sum(x),dp);m=max(0.0_dp,min(1.0_dp-epsilon(1.0_dp),1.0_dp-real(n,dp)/sx))
      allocate(res%param(1));res%param=[m];res%iters=1;res%loglik=-sx+real(n,dp)
      do i=1,n
         if(x(i)>1.and.m>0.0_dp)res%loglik=res%loglik+real(x(i)-1,dp)*log(m*real(x(i),dp))
         res%loglik=res%loglik-log_gamma(real(x(i)+1,dp))
      end do
   end function borel_mle

   function logseries_mle(x,tol,maxiter) result(res)
      integer,intent(in)::x(:);real(dp),intent(in),optional::tol;integer,intent(in),optional::maxiter
      type(mle_result)::res;real(dp)::m,p,p1,a,a2,loga,com,der,der2,eps,sx;integer::n,it,mi,i
      if(size(x)==0.or.any(x<1))then;res%status=1;return;end if
      n=size(x);sx=real(sum(x),dp);m=sx/real(n,dp)
      if(m<=1.0_dp+epsilon(1.0_dp))then;p=tiny(1.0_dp);a=log(p);else;p=1.0_dp-1.0_dp/m;a=log(p/(1.0_dp-p));end if;a2=a
      eps=1e-9_dp;if(present(tol))eps=tol;mi=100;if(present(maxiter))mi=maxiter
      do it=1,mi
         p=1.0_dp/(1.0_dp+exp(-a));p1=1.0_dp-p;loga=log(p1);com=p/loga
         der=m*p1+com;der2=-m*p*p1+p*p1/loga+com*com
         a2=a-der/der2;if(abs(a2-a)<eps)exit;a=a2
      end do
      p=1.0_dp/(1.0_dp+exp(-a2));allocate(res%param(1));res%param=[p];res%iters=it
      res%loglik=sx*log(p)-real(n,dp)*log(-log1p_r(-p))
      do i=1,n;res%loglik=res%loglik-log(real(x(i),dp));end do
   end function logseries_mle

end module rfast_extra_mle
