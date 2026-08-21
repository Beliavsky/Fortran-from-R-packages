! SPDX-License-Identifier: GPL-2.0-or-later
module flexsurv_mixture_full
  use flexsurv_kinds, only : dp
  use flexsurv_fit, only : flexsurv_data, flexsurv_spec, flexsurv_result, &
    fit_flexsurvreg, flexsurv_loglik, initial_theta, parameter_count, &
    predict_survival, predict_density, parameter_row, bfgs_minimize
  use flexsurv_distributions, only : dist_mean, dist_random
  use flexsurv_math, only : logsumexp, near_positive_definite, rng_normal
  use numderiv, only : hessian, grad
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  private

  integer, parameter, public :: mix_method_em = 1
  integer, parameter, public :: mix_method_direct = 2
  integer, parameter, public :: mix_var_direct = 1
  integer, parameter, public :: mix_var_louis = 2

  type, public :: flexsurvmix_full_result
    integer :: k = 0
    integer :: nobs = 0
    integer :: nprob_cov = 0
    type(flexsurv_spec), allocatable :: specs(:)
    type(flexsurv_result), allocatable :: components(:)
    real(dp), allocatable :: alpha(:)
    real(dp), allocatable :: prob_beta(:,:)
    real(dp), allocatable :: prob_x(:,:)
    real(dp), allocatable :: posterior(:,:)
    real(dp), allocatable :: covariance(:,:)
    real(dp), allocatable :: covariance_direct(:,:)
    real(dp), allocatable :: covariance_louis(:,:)
    real(dp) :: loglik = -huge(1.0_dp)
    real(dp) :: aic = huge(1.0_dp)
    integer :: npar = 0
    integer :: iterations = 0
    logical :: converged = .false.
    integer :: status = 1
  end type flexsurvmix_full_result

  public :: fit_flexsurvmix_full, mix_full_loglik, mix_full_posterior
  public :: mix_probabilities, mix_probability_at, mix_full_survival
  public :: mix_full_density, mix_component_mean, mix_component_random
  public :: louis_information, louis_information_numeric, resample_flexsurvmix_full

contains

  subroutine resample_flexsurvmix_full(res,out,seed)
    type(flexsurvmix_full_result),intent(in)::res
    type(flexsurvmix_full_result),intent(out)::out
    integer,intent(in),optional::seed
    real(dp),allocatable::z(:),draw(:),pd(:,:),l(:,:),eps(:)
    integer,allocatable::lens(:)
    integer::j,n,st
    out=res
    if(present(seed))call set_seed_mix(seed)
    if(.not.allocated(res%covariance))return
    allocate(lens(res%k))
    do j=1,res%k
      if(.not.allocated(res%components(j)%theta))return
      lens(j)=size(res%components(j)%theta)
    end do
    z=pack_all(res%alpha,res%prob_beta,res%components,lens)
    n=size(z)
    if(size(res%covariance,1)/=n.or.size(res%covariance,2)/=n)return
    allocate(pd(n,n),l(n,n),eps(n),draw(n))
    call near_positive_definite(res%covariance,pd,1.0e-10_dp)
    call chol_mix(pd,l,st)
    if(st/=0)return
    do j=1,n;eps(j)=rng_normal();end do
    draw=z+matmul(l,eps)
    call unpack_all(draw,out%alpha,out%prob_beta,out%components,lens)
  end subroutine resample_flexsurvmix_full

  subroutine chol_mix(a,l,status)
    real(dp),intent(in)::a(:,:)
    real(dp),intent(out)::l(size(a,1),size(a,2))
    integer,intent(out)::status
    real(dp)::v
    integer::i,j,k,n
    n=size(a,1);l=0.0_dp;status=0
    do i=1,n
      do j=1,i
        v=a(i,j)
        do k=1,j-1;v=v-l(i,k)*l(j,k);end do
        if(i==j)then
          if(v<=0.0_dp)then;status=1;return;end if
          l(i,j)=sqrt(v)
        else
          l(i,j)=v/l(j,j)
        end if
      end do
    end do
  end subroutine chol_mix

  subroutine set_seed_mix(seed)
    integer,intent(in)::seed
    integer::n,i
    integer,allocatable::put(:)
    call random_seed(size=n);allocate(put(n))
    do i=1,n;put(i)=mod(abs(seed)+16127*i,2147483646)+1;end do
    call random_seed(put=put)
  end subroutine set_seed_mix

  function fit_flexsurvmix_full(data,specs,event,prob_x,allowed,component_data, &
      method,var_method,maxit,tol) result(res)
    type(flexsurv_data), intent(in) :: data
    type(flexsurv_spec), intent(in) :: specs(:)
    integer, intent(in), optional :: event(:)
    real(dp), intent(in), optional :: prob_x(:,:)
    logical, intent(in), optional :: allowed(:,:)
    type(flexsurv_data), intent(in), optional :: component_data(:)
    integer, intent(in), optional :: method,var_method,maxit
    real(dp), intent(in), optional :: tol
    type(flexsurvmix_full_result) :: res
    type(flexsurv_data), allocatable :: dcomp(:)
    real(dp), allocatable :: px(:,:),post(:,:),llmat(:,:),pp(:,:),z(:),zhat(:),g(:)
    real(dp), allocatable :: alpha(:),pb(:,:),oldtheta(:,:),par(:),h(:,:),hpd(:,:),cov(:,:)
    integer, allocatable :: ev(:),offset(:),lens(:)
    logical, allocatable :: allow(:,:)
    integer :: n,k,p,j,i,it,mi,meth,vm,st,ntmax
    real(dp) :: oldll,newll,ftol,f

    n=size(data%lower);k=size(specs);p=0
    if(present(prob_x))p=size(prob_x,2)
    res%k=k;res%nobs=n;res%nprob_cov=p
    if(k<1.or.n<1)return
    allocate(ev(n));ev=0
    if(present(event))then
      if(size(event)/=n)return
      ev=event
    end if
    allocate(allow(n,k));allow=.true.
    if(present(allowed))then
      if(size(allowed,1)/=n.or.size(allowed,2)/=k)return
      allow=allowed
    end if
    do i=1,n
      if(ev(i)>=1.and.ev(i)<=k)then
        allow(i,:)=.false.;allow(i,ev(i))=.true.
      end if
      if(.not.any(allow(i,:)))return
    end do
    allocate(px(n,p));if(p>0)px=prob_x
    allocate(dcomp(k));
    if(present(component_data))then
      if(size(component_data)/=k)return
      dcomp=component_data
    else
      do j=1,k;dcomp(j)=data;end do
    end if
    allocate(res%specs(k),res%components(k));res%specs=specs
    allocate(alpha(max(0,k-1)),pb(p,max(0,k-1)))
    alpha=0.0_dp;pb=0.0_dp
    allocate(post(n,k),llmat(n,k),pp(n,k))
    allocate(lens(k),offset(k));ntmax=0
    do j=1,k
      lens(j)=parameter_count(specs(j));ntmax=max(ntmax,lens(j))
    end do
    offset(1)=prob_parameter_count(k,p)+1
    do j=2,k;offset(j)=offset(j-1)+lens(j-1);end do
    res%npar=prob_parameter_count(k,p)+sum(lens)
    meth=mix_method_em;if(present(method))meth=method
    vm=mix_var_louis;if(present(var_method))vm=var_method
    mi=300;if(present(maxit))mi=maxit
    ftol=1.0e-8_dp;if(present(tol))ftol=tol

    ! Component starting fits.  These also improve direct-MLE starting values.
    do j=1,k
      res%components(j)=fit_flexsurvreg(dcomp(j),specs(j),control_maxit=min(mi,150),tol=1.0e-6_dp)
      if(.not.allocated(res%components(j)%theta))then
        res%components(j)%theta=initial_theta(specs(j))
      end if
    end do

    if(meth==mix_method_direct)then
      par=pack_all(alpha,pb,res%components,lens)
      call bfgs_minimize(objective_direct,par,zhat,f,g,res%iterations,st,mi,ftol)
      call unpack_all(zhat,alpha,pb,res%components,lens)
      res%status=st;res%converged=(st==0)
    else
      oldll=-huge(1.0_dp)
      do it=1,mi
        call fill_component_loglik(dcomp,specs,res%components,llmat)
        call probability_matrix(alpha,pb,px,pp)
        call posterior_from_parts(llmat,pp,allow,post)
        call fit_probability_mstep(post,px,alpha,pb,mi,ftol)
        do j=1,k
          dcomp(j)%weights=post(:,j)
          res%components(j)=fit_flexsurvreg(dcomp(j),specs(j),control_maxit=mi,tol=ftol)
        end do
        call fill_component_loglik(dcomp,specs,res%components,llmat,ignore_current_weights=.true.)
        call probability_matrix(alpha,pb,px,pp)
        newll=observed_loglik_from_parts(llmat,pp,allow)
        res%iterations=it
        if(it>1.and.abs(newll-oldll)<=ftol*(1.0_dp+abs(oldll)))then
          res%converged=.true.;res%status=0;exit
        end if
        oldll=newll
      end do
      if(.not.res%converged)res%status=2
    end if

    call fill_component_loglik(dcomp,specs,res%components,llmat,ignore_current_weights=.true.)
    call probability_matrix(alpha,pb,px,pp)
    call posterior_from_parts(llmat,pp,allow,post)
    res%loglik=observed_loglik_from_parts(llmat,pp,allow)
    res%aic=-2.0_dp*res%loglik+2.0_dp*real(res%npar,dp)
    res%alpha=alpha;res%prob_beta=pb;res%prob_x=px;res%posterior=post

    par=pack_all(alpha,pb,res%components,lens)
    allocate(h(size(par),size(par)))
    call hessian(objective_direct,par,h,status=st)
    if(st==0.and.all(ieee_is_finite(h)))then
      allocate(hpd(size(h,1),size(h,2)))
      call near_positive_definite(h,hpd,1.0e-9_dp)
      call invert_matrix_local(hpd,cov,st)
      if(st==0)res%covariance_direct=cov
    end if
    call louis_information(data,dcomp,specs,ev,allow,px,alpha,pb,res%components,post,cov,st)
    if(st==0)res%covariance_louis=cov
    if(vm==mix_var_louis.and.allocated(res%covariance_louis))then
      res%covariance=res%covariance_louis
    else if(allocated(res%covariance_direct))then
      res%covariance=res%covariance_direct
    end if

  contains
    real(dp) function objective_direct(x) result(v)
      real(dp),intent(in)::x(:)
      type(flexsurv_result),allocatable::cr(:)
      real(dp),allocatable::aa(:),bb(:,:),lm(:,:),pm(:,:)
      allocate(cr(k));cr=res%components
      allocate(aa(max(0,k-1)),bb(p,max(0,k-1)),lm(n,k),pm(n,k))
      call unpack_all(x,aa,bb,cr,lens)
      call fill_component_loglik(dcomp,specs,cr,lm,ignore_current_weights=.true.)
      call probability_matrix(aa,bb,px,pm)
      v=-observed_loglik_from_parts(lm,pm,allow)
      if(.not.ieee_is_finite(v))v=1.0e100_dp
    end function objective_direct
  end function fit_flexsurvmix_full

  real(dp) function mix_full_loglik(data,res,event,allowed,component_data) result(ll)
    type(flexsurv_data),intent(in)::data
    type(flexsurvmix_full_result),intent(in)::res
    integer,intent(in),optional::event(:)
    logical,intent(in),optional::allowed(:,:)
    type(flexsurv_data),intent(in),optional::component_data(:)
    type(flexsurv_data),allocatable::dc(:)
    real(dp),allocatable::lm(:,:),pm(:,:)
    logical,allocatable::al(:,:)
    integer::n,k,i
    n=size(data%lower);k=res%k;allocate(dc(k),lm(n,k),pm(n,k),al(n,k));al=.true.
    if(present(component_data))then;dc=component_data;else;do i=1,k;dc(i)=data;end do;end if
    if(present(allowed))al=allowed
    if(present(event))then
      do i=1,n
        if(event(i)>=1.and.event(i)<=k)then;al(i,:)=.false.;al(i,event(i))=.true.;end if
      end do
    end if
    call fill_component_loglik(dc,res%specs,res%components,lm,ignore_current_weights=.true.)
    call probability_matrix(res%alpha,res%prob_beta,res%prob_x,pm)
    ll=observed_loglik_from_parts(lm,pm,al)
  end function mix_full_loglik

  subroutine mix_full_posterior(data,res,post,event,allowed,component_data)
    type(flexsurv_data),intent(in)::data
    type(flexsurvmix_full_result),intent(in)::res
    real(dp),intent(out)::post(size(data%lower),res%k)
    integer,intent(in),optional::event(:)
    logical,intent(in),optional::allowed(:,:)
    type(flexsurv_data),intent(in),optional::component_data(:)
    type(flexsurv_data),allocatable::dc(:)
    real(dp),allocatable::lm(:,:),pm(:,:)
    logical,allocatable::al(:,:)
    integer::n,k,i
    n=size(data%lower);k=res%k;allocate(dc(k),lm(n,k),pm(n,k),al(n,k));al=.true.
    if(present(component_data))then;dc=component_data;else;do i=1,k;dc(i)=data;end do;end if
    if(present(allowed))al=allowed
    if(present(event))then
      do i=1,n
        if(event(i)>=1.and.event(i)<=k)then;al(i,:)=.false.;al(i,event(i))=.true.;end if
      end do
    end if
    call fill_component_loglik(dc,res%specs,res%components,lm,ignore_current_weights=.true.)
    call probability_matrix(res%alpha,res%prob_beta,res%prob_x,pm)
    call posterior_from_parts(lm,pm,al,post)
  end subroutine mix_full_posterior

  subroutine mix_probabilities(res,x,prob)
    type(flexsurvmix_full_result),intent(in)::res
    real(dp),intent(in)::x(:,:)
    real(dp),intent(out)::prob(size(x,1),res%k)
    call probability_matrix(res%alpha,res%prob_beta,x,prob)
  end subroutine mix_probabilities

  subroutine mix_probability_at(res,xrow,prob)
    type(flexsurvmix_full_result),intent(in)::res
    real(dp),intent(in)::xrow(:)
    real(dp),intent(out)::prob(res%k)
    real(dp)::xx(1,size(xrow)),pp(1,res%k)
    xx(1,:)=xrow;call probability_matrix(res%alpha,res%prob_beta,xx,pp);prob=pp(1,:)
  end subroutine mix_probability_at

  real(dp) function mix_full_survival(res,row,t,xprob) result(s)
    type(flexsurvmix_full_result),intent(in)::res
    integer,intent(in)::row
    real(dp),intent(in)::t
    real(dp),intent(in),optional::xprob(:)
    real(dp)::pr(res%k)
    integer::j
    call prediction_probs(res,row,xprob,pr);s=0.0_dp
    do j=1,res%k;s=s+pr(j)*predict_survival(res%specs(j),res%components(j)%theta,row,t);end do
  end function mix_full_survival

  real(dp) function mix_full_density(res,row,t,xprob) result(f)
    type(flexsurvmix_full_result),intent(in)::res
    integer,intent(in)::row
    real(dp),intent(in)::t
    real(dp),intent(in),optional::xprob(:)
    real(dp)::pr(res%k)
    integer::j
    call prediction_probs(res,row,xprob,pr);f=0.0_dp
    do j=1,res%k;f=f+pr(j)*predict_density(res%specs(j),res%components(j)%theta,row,t);end do
  end function mix_full_density

  real(dp) function mix_component_mean(res,component,row) result(v)
    type(flexsurvmix_full_result),intent(in)::res
    integer,intent(in)::component,row
    real(dp),allocatable::pa(:)
    pa=parameter_row(res%specs(component),res%components(component)%theta,row)
    v=dist_mean(res%specs(component)%dist,pa)
  end function mix_component_mean

  real(dp) function mix_component_random(res,component,row) result(v)
    type(flexsurvmix_full_result),intent(in)::res
    integer,intent(in)::component,row
    real(dp),allocatable::pa(:)
    pa=parameter_row(res%specs(component),res%components(component)%theta,row)
    v=dist_random(res%specs(component)%dist,pa)
  end function mix_component_random

  subroutine louis_information(data,dcomp,specs,event,allowed,px,alpha,pb,components,post,cov,status)
    type(flexsurv_data),intent(in)::data,dcomp(:)
    type(flexsurv_spec),intent(in)::specs(:)
    integer,intent(in)::event(:)
    logical,intent(in)::allowed(:,:)
    real(dp),intent(in)::px(:,:),alpha(:),pb(:,:),post(:,:)
    type(flexsurv_result),intent(in)::components(:)
    real(dp),allocatable,intent(out)::cov(:,:)
    integer,intent(out)::status
    integer::k,n,p,np,j,i,c,d,r,sr,npar,st,offj
    integer,allocatable::lens(:),off(:)
    real(dp),allocatable::pp(:,:),score(:,:),gbar(:),info(:,:),hfull(:,:), &
      gloc(:),hloc(:,:),ipd(:,:),feat(:),li(:)
    real(dp)::pc,pd,coef
    k=size(specs);n=size(data%lower);p=size(px,2);np=(k-1)*(p+1)
    allocate(lens(k),off(k));offj=np
    do j=1,k
      lens(j)=size(components(j)%theta);off(j)=offj;offj=offj+lens(j)
    end do
    npar=offj;allocate(pp(n,k),score(npar,k),gbar(npar),info(npar,npar), &
      hfull(npar,npar),feat(p+1),li(n))
    call probability_matrix(alpha,pb,px,pp);info=0.0_dp;status=0
    do i=1,n
      feat(1)=1.0_dp;if(p>0)feat(2:)=px(i,:);score=0.0_dp;gbar=0.0_dp
      do j=1,k
        if(post(i,j)<=0.0_dp)cycle
        call membership_score(j,pp(i,:),feat,score(:,j))
        allocate(gloc(lens(j)),hloc(lens(j),lens(j)))
        call grad(compfun,components(j)%theta,gloc,status=st)
        if(st/=0)then;status=st;return;end if
        call hessian(compfun,components(j)%theta,hloc,status=st)
        if(st/=0)then;status=st;return;end if
        score(off(j)+1:off(j)+lens(j),j)=gloc
        hfull=0.0_dp
        call membership_hessian(pp(i,:),feat,hfull(1:np,1:np))
        hfull(off(j)+1:off(j)+lens(j),off(j)+1:off(j)+lens(j))=hloc
        info=info-post(i,j)*hfull
        gbar=gbar+post(i,j)*score(:,j)
        deallocate(gloc,hloc)
      end do
      do j=1,k
        if(post(i,j)>0.0_dp) &
          info=info-post(i,j)*outer(score(:,j)-gbar,score(:,j)-gbar)
      end do
    end do
    allocate(ipd(npar,npar));call near_positive_definite(info,ipd,1.0e-9_dp)
    call invert_matrix_local(ipd,cov,status)
  contains
    real(dp) function compfun(z) result(v)
      real(dp),intent(in)::z(:)
      call component_individual_unweighted(dcomp(j),specs(j),z,li);v=li(i)
    end function compfun
    subroutine membership_score(comp,prob,feature,sc)
      integer,intent(in)::comp
      real(dp),intent(in)::prob(:),feature(:)
      real(dp),intent(inout)::sc(:)
      integer::cc,rr,idx
      do cc=2,k
        coef=merge(1.0_dp,0.0_dp,comp==cc)-prob(cc)
        idx=cc-1;sc(idx)=coef
        do rr=2,p+1
          idx=(k-1)+(cc-2)*p+(rr-1);sc(idx)=coef*feature(rr)
        end do
      end do
    end subroutine membership_score
    subroutine membership_hessian(prob,feature,hh)
      real(dp),intent(in)::prob(:),feature(:)
      real(dp),intent(out)::hh(:,:)
      integer::cc,dd,rr,ss,ic,id
      hh=0.0_dp
      do cc=2,k
        pc=prob(cc)
        do dd=2,k
          pd=prob(dd);coef=-pc*(merge(1.0_dp,0.0_dp,cc==dd)-pd)
          do rr=1,p+1
            if(rr==1)then;ic=cc-1;else;ic=(k-1)+(cc-2)*p+(rr-1);end if
            do ss=1,p+1
              if(ss==1)then;id=dd-1;else;id=(k-1)+(dd-2)*p+(ss-1);end if
              hh(ic,id)=coef*feature(rr)*feature(ss)
            end do
          end do
        end do
      end do
    end subroutine membership_hessian
  end subroutine louis_information

  subroutine louis_information_numeric(data,dcomp,specs,event,allowed,px,alpha,pb,components,post,cov,status)
    type(flexsurv_data),intent(in)::data,dcomp(:)
    type(flexsurv_spec),intent(in)::specs(:)
    integer,intent(in)::event(:)
    logical,intent(in)::allowed(:,:)
    real(dp),intent(in)::px(:,:),alpha(:),pb(:,:),post(:,:)
    type(flexsurv_result),intent(in)::components(:)
    real(dp),allocatable,intent(out)::cov(:,:)
    integer,intent(out)::status
    integer::k,n,p,j,i,q,npar,st
    integer,allocatable::lens(:)
    real(dp),allocatable::full(:),score(:,:,:),hs(:,:),info(:,:),gbar(:),hh(:,:),hneg(:,:),ipd(:,:)
    k=size(specs);n=size(data%lower);p=size(px,2);allocate(lens(k))
    do j=1,k;lens(j)=size(components(j)%theta);end do
    full=pack_all(alpha,pb,components,lens);npar=size(full)
    allocate(score(n,npar,k),hs(npar,npar),info(npar,npar),gbar(npar),hh(npar,npar))
    score=0.0_dp;info=0.0_dp
    do i=1,n
      gbar=0.0_dp
      do j=1,k
        if(post(i,j)<=0.0_dp)cycle
        call grad(cfun,full,score(i,:,j),status=st)
        if(st/=0)then;status=st;return;end if
        call hessian(cfun,full,hh,status=st)
        if(st/=0)then;status=st;return;end if
        info=info-post(i,j)*hh
        gbar=gbar+post(i,j)*score(i,:,j)
      end do
      do j=1,k
        if(post(i,j)<=0.0_dp)cycle
        info=info-post(i,j)*outer(score(i,:,j)-gbar,score(i,:,j)-gbar)
      end do
    end do
    allocate(ipd(npar,npar))
    call near_positive_definite(info,ipd,1.0e-9_dp)
    call invert_matrix_local(ipd,cov,status)
  contains
    real(dp) function cfun(x) result(v)
      real(dp),intent(in)::x(:)
      type(flexsurv_result),allocatable::cr(:)
      real(dp),allocatable::aa(:),bb(:,:),pp(:,:),li(:)
      integer::jj
      allocate(cr(k),aa(max(0,k-1)),bb(p,max(0,k-1)),pp(n,k),li(n));cr=components
      call unpack_all(x,aa,bb,cr,lens);call probability_matrix(aa,bb,px,pp)
      call component_individual_unweighted(dcomp(j),specs(j),cr(j)%theta,li)
      v=log(max(pp(i,j),tiny(1.0_dp)))+li(i)
    end function cfun
  end subroutine louis_information_numeric

  subroutine fit_probability_mstep(post,x,alpha,beta,maxit,tol)
    real(dp),intent(in)::post(:,:),x(:,:)
    real(dp),intent(inout)::alpha(:),beta(:,:)
    integer,intent(in)::maxit
    real(dp),intent(in)::tol
    real(dp),allocatable::z0(:),zh(:),g(:)
    real(dp)::f
    integer::st,it,k,p
    k=size(post,2);p=size(x,2)
    z0=pack_prob(alpha,beta)
    if(size(z0)==0)return
    call bfgs_minimize(pobj,z0,zh,f,g,it,st,maxit,tol)
    if(st==0.or.all(ieee_is_finite(zh)))call unpack_prob(zh,k,p,alpha,beta)
  contains
    real(dp) function pobj(z) result(v)
      real(dp),intent(in)::z(:)
      real(dp),allocatable::aa(:),bb(:,:),pm(:,:)
      integer::i,j
      allocate(aa(max(0,k-1)),bb(p,max(0,k-1)),pm(size(post,1),k))
      call unpack_prob(z,k,p,aa,bb);call probability_matrix(aa,bb,x,pm);v=0.0_dp
      do i=1,size(post,1);do j=1,k
        if(post(i,j)>0.0_dp)v=v-post(i,j)*log(max(pm(i,j),tiny(1.0_dp)))
      end do;end do
    end function pobj
  end subroutine fit_probability_mstep

  subroutine fill_component_loglik(dc,specs,comp,llmat,ignore_current_weights)
    type(flexsurv_data),intent(in)::dc(:)
    type(flexsurv_spec),intent(in)::specs(:)
    type(flexsurv_result),intent(in)::comp(:)
    real(dp),intent(out)::llmat(:,:)
    logical,intent(in),optional::ignore_current_weights
    integer::j
    do j=1,size(specs)
      call component_individual_unweighted(dc(j),specs(j),comp(j)%theta,llmat(:,j))
    end do
  end subroutine fill_component_loglik

  subroutine component_individual_unweighted(data,spec,theta,li)
    type(flexsurv_data),intent(in)::data
    type(flexsurv_spec),intent(in)::spec
    real(dp),intent(in)::theta(:)
    real(dp),intent(out)::li(:)
    type(flexsurv_data)::du
    real(dp)::dummy
    du=data;du%weights=1.0_dp
    dummy=flexsurv_loglik(du,spec,theta,individual=li)
  end subroutine component_individual_unweighted

  subroutine probability_matrix(alpha,beta,x,pm)
    real(dp),intent(in)::alpha(:),beta(:,:),x(:,:)
    real(dp),intent(out)::pm(size(x,1),size(alpha)+1)
    real(dp),allocatable::eta(:)
    real(dp)::mx,den
    integer::i,j,k,p
    k=size(alpha)+1;p=size(x,2);allocate(eta(k))
    do i=1,size(x,1)
      eta=0.0_dp
      do j=2,k
        eta(j)=alpha(j-1)
        if(p>0)eta(j)=eta(j)+dot_product(x(i,:),beta(:,j-1))
      end do
      mx=maxval(eta);eta=exp(eta-mx);den=sum(eta);pm(i,:)=eta/den
    end do
  end subroutine probability_matrix

  subroutine posterior_from_parts(llmat,pm,allowed,post)
    real(dp),intent(in)::llmat(:,:),pm(:,:)
    logical,intent(in)::allowed(:,:)
    real(dp),intent(out)::post(size(llmat,1),size(llmat,2))
    real(dp),allocatable::lw(:)
    real(dp)::ls
    integer::i,j,k
    k=size(llmat,2);allocate(lw(k))
    do i=1,size(llmat,1)
      lw=-huge(1.0_dp)
      do j=1,k
        if(allowed(i,j))lw(j)=log(max(pm(i,j),tiny(1.0_dp)))+llmat(i,j)
      end do
      ls=logsumexp(lw);post(i,:)=0.0_dp
      do j=1,k;if(allowed(i,j))post(i,j)=exp(lw(j)-ls);end do
    end do
  end subroutine posterior_from_parts

  real(dp) function observed_loglik_from_parts(llmat,pm,allowed) result(ll)
    real(dp),intent(in)::llmat(:,:),pm(:,:)
    logical,intent(in)::allowed(:,:)
    real(dp),allocatable::lw(:)
    integer::i,j,k
    k=size(llmat,2);allocate(lw(k));ll=0.0_dp
    do i=1,size(llmat,1)
      lw=-huge(1.0_dp)
      do j=1,k
        if(allowed(i,j))lw(j)=log(max(pm(i,j),tiny(1.0_dp)))+llmat(i,j)
      end do
      ll=ll+logsumexp(lw)
    end do
  end function observed_loglik_from_parts

  subroutine prediction_probs(res,row,xprob,pr)
    type(flexsurvmix_full_result),intent(in)::res
    integer,intent(in)::row
    real(dp),intent(in),optional::xprob(:)
    real(dp),intent(out)::pr(res%k)
    real(dp),allocatable::xx(:,:),pp(:,:)
    if(present(xprob))then
      allocate(xx(1,size(xprob)),pp(1,res%k));xx(1,:)=xprob
    else if(allocated(res%prob_x).and.row>=1.and.row<=size(res%prob_x,1))then
      allocate(xx(1,size(res%prob_x,2)),pp(1,res%k));xx(1,:)=res%prob_x(row,:)
    else
      allocate(xx(1,res%nprob_cov),pp(1,res%k));xx=0.0_dp
    end if
    call probability_matrix(res%alpha,res%prob_beta,xx,pp);pr=pp(1,:)
  end subroutine prediction_probs

  integer pure function prob_parameter_count(k,p) result(n)
    integer,intent(in)::k,p;n=max(0,k-1)*(p+1)
  end function prob_parameter_count

  function pack_prob(alpha,beta) result(z)
    real(dp),intent(in)::alpha(:),beta(:,:)
    real(dp),allocatable::z(:)
    integer::k1,p,j,pos
    k1=size(alpha);p=size(beta,1);allocate(z(k1*(p+1)));pos=0
    if(k1>0)then;z(1:k1)=alpha;pos=k1;end if
    do j=1,k1
      if(p>0)then;z(pos+1:pos+p)=beta(:,j);pos=pos+p;end if
    end do
  end function pack_prob

  subroutine unpack_prob(z,k,p,alpha,beta)
    real(dp),intent(in)::z(:)
    integer,intent(in)::k,p
    real(dp),intent(out)::alpha(max(0,k-1)),beta(p,max(0,k-1))
    integer::j,pos,k1
    k1=max(0,k-1);alpha=0.0_dp;beta=0.0_dp;pos=0
    if(k1>0)then;alpha=z(1:k1);pos=k1;end if
    do j=1,k1
      if(p>0)then;beta(:,j)=z(pos+1:pos+p);pos=pos+p;end if
    end do
  end subroutine unpack_prob

  function pack_all(alpha,beta,comp,lens) result(z)
    real(dp),intent(in)::alpha(:),beta(:,:)
    type(flexsurv_result),intent(in)::comp(:)
    integer,intent(in)::lens(:)
    real(dp),allocatable::z(:),zp(:)
    integer::j,pos
    zp=pack_prob(alpha,beta);allocate(z(size(zp)+sum(lens)));pos=0
    if(size(zp)>0)then;z(1:size(zp))=zp;pos=size(zp);end if
    do j=1,size(comp);z(pos+1:pos+lens(j))=comp(j)%theta;pos=pos+lens(j);end do
  end function pack_all

  subroutine unpack_all(z,alpha,beta,comp,lens)
    real(dp),intent(in)::z(:)
    real(dp),intent(out)::alpha(:),beta(:,:)
    type(flexsurv_result),intent(inout)::comp(:)
    integer,intent(in)::lens(:)
    integer::p,k1,j,pos,np
    k1=size(alpha);p=size(beta,1);np=k1*(p+1);call unpack_prob(z(1:np),k1+1,p,alpha,beta);pos=np
    do j=1,size(comp)
      if(.not.allocated(comp(j)%theta))allocate(comp(j)%theta(lens(j)))
      comp(j)%theta=z(pos+1:pos+lens(j));pos=pos+lens(j)
    end do
  end subroutine unpack_all

  pure function outer(a,b) result(c)
    real(dp),intent(in)::a(:),b(:)
    real(dp)::c(size(a),size(b))
    integer::i,j
    do i=1,size(a);do j=1,size(b);c(i,j)=a(i)*b(j);end do;end do
  end function outer

  subroutine invert_matrix_local(a,ainv,status)
    real(dp),intent(in)::a(:,:)
    real(dp),allocatable,intent(out)::ainv(:,:)
    integer,intent(out)::status
    real(dp),allocatable::aug(:,:),tmp(:)
    real(dp)::piv,fac
    integer::n,i,j,k,ip
    n=size(a,1);allocate(aug(n,2*n),tmp(2*n));aug=0.0_dp;aug(:,1:n)=a
    do i=1,n;aug(i,n+i)=1.0_dp;end do
    status=0
    do i=1,n
      ip=i;do k=i+1,n;if(abs(aug(k,i))>abs(aug(ip,i)))ip=k;end do
      if(abs(aug(ip,i))<1.0e-12_dp)then;status=1;allocate(ainv(n,n));ainv=0.0_dp;return;end if
      if(ip/=i)then;tmp=aug(i,:);aug(i,:)=aug(ip,:);aug(ip,:)=tmp;end if
      piv=aug(i,i);aug(i,:)=aug(i,:)/piv
      do j=1,n;if(j/=i)then;fac=aug(j,i);aug(j,:)=aug(j,:)-fac*aug(i,:);end if;end do
    end do
    allocate(ainv(n,n));ainv=aug(:,n+1:)
  end subroutine invert_matrix_local

end module flexsurv_mixture_full
