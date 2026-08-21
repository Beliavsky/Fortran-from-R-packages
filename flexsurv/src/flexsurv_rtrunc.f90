! SPDX-License-Identifier: GPL-2.0-or-later
module flexsurv_rtrunc
  use flexsurv_kinds, only : dp
  use flexsurv_distributions, only : dist_npar, dist_logpdf, dist_cdf, &
    dist_gamma
  use flexsurv_fit, only : bfgs_minimize
  use flexsurv_math, only : integrate_gauss_legendre, normal_quantile, near_positive_definite, rng_normal
  use numderiv, only : hessian
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  private

  integer,parameter,public::rtrunc_joint=1,rtrunc_final=2

  type,public::survrtrunc_result
    real(dp),allocatable::time(:),surv(:),se_surv(:),se_log(:),lower(:),upper(:)
  end type survrtrunc_result

  type,public::flexsurvrtrunc_result
    integer::dist=0,method=rtrunc_joint
    real(dp),allocatable::parameters(:),theta_vec(:),covariance(:,:),gradient(:)
    logical,allocatable::fixed(:)
    real(dp)::theta=0.0_dp,loglik=-huge(1.0_dp),aic=huge(1.0_dp)
    integer::iterations=0,status=1
    logical::converged=.false.
  end type flexsurvrtrunc_result

  public::survrtrunc_fit,fit_flexsurvrtrunc,flexsurvrtrunc_loglik,rtrunc_parameter_draws

contains

  function survrtrunc_fit(t,rtrunc,tmax,eps,conf_level) result(res)
    real(dp),intent(in)::t(:),rtrunc(:),tmax
    real(dp),intent(in),optional::eps,conf_level
    type(survrtrunc_result)::res
    real(dp),allocatable::entry(:),trev(:),evtime(:),sl(:),d(:),risk(:),h(:),cdf(:)
    real(dp)::ee,cl,prod,z,varsum,lval
    integer::n,m,i,j,k
    n=size(t);ee=0.001_dp;if(present(eps))ee=eps
    cl=0.95_dp;if(present(conf_level))cl=conf_level
    allocate(entry(n),trev(n));entry=tmax-min(rtrunc,tmax)-ee;trev=tmax-t
    call unique_sorted(trev,evtime);m=size(evtime)
    allocate(sl(m),d(m),risk(m));prod=1.0_dp
    do j=1,m
      d(j)=real(count(abs(trev-evtime(j))<=16.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(evtime(j)))),dp)
      risk(j)=real(count(entry<evtime(j).and.trev>=evtime(j)),dp)
      if(risk(j)>0.0_dp)prod=prod*(1.0_dp-d(j)/risk(j))
      sl(j)=prod
    end do
    allocate(res%time(m),res%surv(m),res%se_surv(m),res%se_log(m),res%lower(m),res%upper(m))
    do j=1,m;res%time(j)=tmax-evtime(m-j+1);end do
    if(m==1)then
      res%surv=0.0_dp
    else
      do j=1,m-1;res%surv(j)=1.0_dp-sl(m-j);end do
      res%surv(m)=0.0_dp
    end if
    allocate(h(m));h=0.0_dp
    do j=1,m
      if(risk(j)>d(j).and.risk(j)>0.0_dp)h(j)=d(j)/(risk(j)*(risk(j)-d(j)))
    end do
    z=normal_quantile(0.5_dp*(1.0_dp-cl))
    res%se_surv=0.0_dp;res%se_log=0.0_dp;res%lower=0.0_dp;res%upper=1.0_dp
    do j=1,m-1
      cdf= [res%surv(j)] ! keeps notation close to upstream; cdf is 1-S after reversal
      varsum=0.0_dp
      do k=1,m-j-1;varsum=varsum+h(k);end do
      if(res%surv(j)>0.0_dp.and.res%surv(j)<1.0_dp)then
        res%se_surv(j)=res%surv(j)*sqrt(max(varsum,0.0_dp))
        if(abs(log(res%surv(j)))>tiny(1.0_dp)) &
          res%se_log(j)=sqrt(max(varsum,0.0_dp))/abs(log(res%surv(j)))
        lval=log(-log(max(res%surv(j),tiny(1.0_dp))))
        res%lower(j)=1.0_dp-exp(-exp(lval+z*res%se_log(j)))
        res%upper(j)=1.0_dp-exp(-exp(lval-z*res%se_log(j)))
      end if
    end do
    res%se_surv(m)=0.0_dp;res%se_log(m)=0.0_dp;res%lower(m)=0.0_dp;res%upper(m)=0.0_dp
  end function survrtrunc_fit

  function fit_flexsurvrtrunc(t,tinit,rtrunc,tmax,dist,par_init,theta_init, &
      method,fixed_theta,maxit,tol,fixed_params) result(res)
    real(dp),intent(in)::t(:),tinit(:),rtrunc(:),tmax,par_init(:),theta_init
    integer,intent(in)::dist
    integer,intent(in),optional::method,maxit
    logical,intent(in),optional::fixed_theta,fixed_params(:)
    real(dp),intent(in),optional::tol
    type(flexsurvrtrunc_result)::res
    real(dp),allocatable::zfull(:),free0(:),freehat(:),g(:),hh(:,:),hpd(:,:),inv(:,:)
    logical,allocatable::fix(:)
    integer,allocatable::freeidx(:)
    real(dp)::f,tt
    integer::nb,meth,mi,ist,nfree,i,j
    logical::ft
    nb=dist_npar(dist);meth=rtrunc_joint;if(present(method))meth=method
    ft=.true.;if(present(fixed_theta))ft=fixed_theta
    if(meth==rtrunc_final)ft=.true.
    allocate(zfull(nb+1),fix(nb+1));fix=.false.
    do i=1,nb;zfull(i)=transform_param(dist,i,par_init(i));end do
    zfull(nb+1)=theta_init
    if(present(fixed_params))then
      if(size(fixed_params)==nb)fix(1:nb)=fixed_params
    end if
    fix(nb+1)=ft;nfree=count(.not.fix);allocate(freeidx(nfree));j=0
    do i=1,nb+1
      if(.not.fix(i))then;j=j+1;freeidx(j)=i;end if
    end do
    mi=500;if(present(maxit))mi=maxit;tt=1.0e-7_dp;if(present(tol))tt=tol
    if(nfree>0)then
      allocate(free0(nfree));free0=zfull(freeidx)
      call bfgs_minimize(objective,free0,freehat,f,g,res%iterations,res%status,mi,tt)
      zfull(freeidx)=freehat
    else
      f=-flexsurvrtrunc_loglik(t,tinit,rtrunc,tmax,dist,par_init,theta_init,meth)
      allocate(g(0));res%status=0;res%iterations=0
    end if
    allocate(res%parameters(nb));do i=1,nb;res%parameters(i)=inverse_param(dist,i,zfull(i));end do
    res%theta=zfull(nb+1);res%theta_vec=zfull;res%fixed=fix
    res%dist=dist;res%method=meth;res%gradient=g;res%loglik=-f
    res%aic=-2.0_dp*res%loglik+2.0_dp*real(nfree,dp);res%converged=res%status==0
    allocate(res%covariance(nb+1,nb+1));res%covariance=0.0_dp
    if(nfree>0)then
      allocate(hh(nfree,nfree));call hessian(objective,freehat,hh)
      allocate(hpd(nfree,nfree));call near_positive_definite(hh,hpd,1.0e-10_dp)
      call invert_local(hpd,inv,ist)
      if(ist==0)then
        do i=1,nfree;do j=1,nfree
          res%covariance(freeidx(i),freeidx(j))=inv(i,j)
        end do;end do
      end if
    end if
  contains
    real(dp) function objective(zfree) result(v)
      real(dp),intent(in)::zfree(:)
      real(dp),allocatable::par(:),z(:)
      integer::jj
      z=zfull;z(freeidx)=zfree;allocate(par(nb))
      do jj=1,nb;par(jj)=inverse_param(dist,jj,z(jj));end do
      v=-flexsurvrtrunc_loglik(t,tinit,rtrunc,tmax,dist,par,z(nb+1),meth)
      if(.not.ieee_is_finite(v))v=1.0e100_dp
    end function objective
  end function fit_flexsurvrtrunc

  subroutine rtrunc_parameter_draws(res,nsim,draws,seed)
    type(flexsurvrtrunc_result),intent(in)::res
    integer,intent(in)::nsim
    real(dp),allocatable,intent(out)::draws(:,:)
    integer,intent(in),optional::seed
    real(dp),allocatable::pd(:,:),l(:,:),z(:),work(:)
    integer::n,i,k,st,nb
    if(present(seed))call set_seed_local(seed)
    n=size(res%theta_vec);nb=n-1;allocate(draws(n,nsim),pd(n,n),l(n,n),z(n),work(n))
    if(.not.allocated(res%covariance))then
      do k=1,nsim;draws(1:nb,k)=res%parameters;draws(n,k)=res%theta;end do
      return
    end if
    call near_positive_definite(res%covariance,pd,1.0e-12_dp);call chol_lower_local(pd,l,st)
    do k=1,nsim
      if(st==0)then
        do i=1,n;z(i)=rng_normal();end do;work=res%theta_vec+matmul(l,z)
        if(allocated(res%fixed))where(res%fixed)work=res%theta_vec
      else;work=res%theta_vec;end if
      do i=1,nb;draws(i,k)=inverse_param(res%dist,i,work(i));end do
      draws(n,k)=work(n)
    end do
  end subroutine rtrunc_parameter_draws

  subroutine chol_lower_local(a,l,status)
    real(dp),intent(in)::a(:,:);real(dp),intent(out)::l(size(a,1),size(a,2));integer,intent(out)::status
    real(dp)::ss;integer::i,j,k,n
    n=size(a,1);l=0.0_dp;status=0
    do i=1,n;do j=1,i
      ss=a(i,j);do k=1,j-1;ss=ss-l(i,k)*l(j,k);end do
      if(i==j)then;if(ss<=0.0_dp)then;status=1;return;end if;l(i,j)=sqrt(ss)
      else;l(i,j)=ss/l(j,j);end if
    end do;end do
  end subroutine chol_lower_local

  subroutine set_seed_local(seed)
    integer,intent(in)::seed;integer::n,i;integer,allocatable::put(:)
    call random_seed(size=n);allocate(put(n));do i=1,n;put(i)=mod(abs(seed)+65537*i,2147483646)+1;end do
    call random_seed(put=put)
  end subroutine set_seed_local

  real(dp) function flexsurvrtrunc_loglik(t,tinit,rtrunc,tmax,dist,par,theta,method) result(ll)
    real(dp),intent(in)::t(:),tinit(:),rtrunc(:),tmax,par(:),theta
    integer,intent(in)::dist,method
    real(dp)::integ,y
    integer::i
    ll=0.0_dp
    select case(method)
    case(rtrunc_joint)
      integ=integrate_gauss_legendre(joint_integrand,0.0_dp,tmax,96)
      if(integ<=0.0_dp)then;ll=-huge(1.0_dp);return;end if
      do i=1,size(t);ll=ll+theta*tinit(i)+dist_logpdf(dist,t(i),par);end do
      ll=ll-real(size(t),dp)*log(integ)
    case(rtrunc_final)
      if(dist==dist_gamma)then
        do i=1,size(t)
          y=tinit(i)+t(i)
          ll=ll+dist_logpdf(dist_gamma,t(i),[par(1),par(2)+theta]) &
            -log(max(dist_cdf(dist_gamma,y,[par(1),par(2)+theta]),tiny(1.0_dp)))
        end do
      else
        do i=1,size(t)
          y=tinit(i)+t(i)
          integ=integrate_gauss_legendre(final_integrand,0.0_dp,y,96)
          if(integ<=0.0_dp)then;ll=-huge(1.0_dp);return;end if
          ll=ll+dist_logpdf(dist,t(i),par)-theta*t(i)-log(integ)
        end do
      end if
    case default
      ll=-huge(1.0_dp)
    end select
  contains
    real(dp) function joint_integrand(a) result(v)
      real(dp),intent(in)::a
      v=exp(min(theta*a,700.0_dp))*dist_cdf(dist,tmax-a,par)
    end function joint_integrand
    real(dp) function final_integrand(a) result(v)
      real(dp),intent(in)::a
      v=exp(dist_logpdf(dist,a,par)-theta*a)
    end function final_integrand
  end function flexsurvrtrunc_loglik

  pure real(dp) function transform_param(dist,i,x) result(z)
    integer,intent(in)::dist,i;real(dp),intent(in)::x
    if(positive_param(dist,i))then;z=log(max(x,tiny(1.0_dp)));else;z=x;end if
  end function transform_param
  pure real(dp) function inverse_param(dist,i,z) result(x)
    integer,intent(in)::dist,i;real(dp),intent(in)::z
    if(positive_param(dist,i))then;x=exp(min(z,700.0_dp));else;x=z;end if
  end function inverse_param
  pure logical function positive_param(dist,i) result(p)
    use flexsurv_distributions, only : dist_exponential,dist_weibull,dist_weibull_ph,dist_gamma, &
      dist_lognormal,dist_gompertz,dist_loglogistic,dist_gengamma,dist_genf
    integer,intent(in)::dist,i
    select case(dist)
    case(dist_exponential);p=.true.
    case(dist_weibull,dist_weibull_ph,dist_gamma,dist_loglogistic);p=.true.
    case(dist_lognormal);p=i==2
    case(dist_gompertz);p=i==2
    case(dist_gengamma);p=i==2
    case(dist_genf);p=i==2.or.i==4
    case default;p=.false.
    end select
  end function positive_param

  subroutine unique_sorted(x,u)
    real(dp),intent(in)::x(:);real(dp),allocatable,intent(out)::u(:)
    real(dp),allocatable::y(:),tmp(:);real(dp)::key,tol
    integer::n,i,j,m
    n=size(x);allocate(y(n));y=x
    do i=2,n
      key=y(i);j=i-1
      do while(j>=1)
        if(y(j)<=key)exit
        y(j+1)=y(j);j=j-1
      end do
      y(j+1)=key
    end do
    allocate(tmp(n));m=0
    do i=1,n
      tol=16.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(y(i)))
      if(m==0)then;m=1;tmp(m)=y(i)
      else if(abs(y(i)-tmp(m))>tol)then;m=m+1;tmp(m)=y(i);end if
    end do
    allocate(u(m));u=tmp(1:m)
  end subroutine unique_sorted

  subroutine invert_local(a,ainv,status)
    real(dp),intent(in)::a(:,:);real(dp),allocatable,intent(out)::ainv(:,:);integer,intent(out)::status
    real(dp),allocatable::aug(:,:),tmp(:);real(dp)::piv,fac;integer::n,i,j,k,ip
    n=size(a,1);allocate(aug(n,2*n),tmp(2*n));aug=0.0_dp;aug(:,1:n)=a
    do i=1,n;aug(i,n+i)=1.0_dp;end do;status=0
    do i=1,n
      ip=i;do k=i+1,n;if(abs(aug(k,i))>abs(aug(ip,i)))ip=k;end do
      if(abs(aug(ip,i))<1.0e-12_dp)then;status=1;allocate(ainv(n,n));ainv=0.0_dp;return;end if
      if(ip/=i)then;tmp=aug(i,:);aug(i,:)=aug(ip,:);aug(ip,:)=tmp;end if
      piv=aug(i,i);aug(i,:)=aug(i,:)/piv
      do j=1,n;if(j/=i)then;fac=aug(j,i);aug(j,:)=aug(j,:)-fac*aug(i,:);end if;end do
    end do
    allocate(ainv(n,n));ainv=aug(:,n+1:)
  end subroutine invert_local

end module flexsurv_rtrunc
