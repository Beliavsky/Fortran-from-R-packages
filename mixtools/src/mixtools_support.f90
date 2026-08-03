! SPDX-License-Identifier: GPL-2.0-or-later
module mixtools_support
  use mixtools_kinds, only : dp
  use mixtools_status
  use mixtools_types
  use mixtools_rng
  use mixtools_distributions, only : normalize_logweights, normal_logpdf
  use mixtools_distributions, only : rmvnorm, rnormmix
  use mixtools_linalg, only : inverse_spd
  use mixtools_parametric, only : normalmix_em, multmix_em, repnormmix_em
  use mixtools_regression, only : regmix_em
  use mixtools_utilities, only : wquantile
  implicit none
  private
  public :: make_constraints, makemultdata, normalmix_model_selection
  public :: regmix_model_selection, multmix_model_selection, repnormmix_model_selection
  public :: normalmix_boot_comp, normalmix_boot_se, regmix_mh, posterior_beta_intervals
  public :: likelihood_ratio_pvalue, mixture_information_criteria
contains
  subroutine make_constraints(category,alpha,result)
    integer,intent(in)::category(:)
    real(dp),intent(in),optional::alpha(:)
    type(constraint_result),intent(out)::result
    if(any(category<0))then;result%status=MIXTOOLS_INVALID_ARGUMENT;return;end if
    result%category=category;allocate(result%alpha(size(category)));result%alpha=1.0_dp
    if(present(alpha))then
      if(size(alpha)/=size(category))then;result%status=MIXTOOLS_DIMENSION_ERROR;return;end if
      result%alpha=alpha
    end if
    result%status=MIXTOOLS_SUCCESS
  end subroutine make_constraints

  subroutine makemultdata(x,cuts,y,status)
    real(dp),intent(in)::x(:,:)
    real(dp),intent(in)::cuts(:)
    real(dp),allocatable,intent(out)::y(:,:)
    integer,intent(out)::status
    integer::i,j,c
    if(size(cuts)<1)then;allocate(y(0,0));status=MIXTOOLS_INVALID_ARGUMENT;return;end if
    allocate(y(size(x,1),size(cuts)+1));y=0.0_dp
    do i=1,size(x,1);do j=1,size(x,2)
      c=1;do while(c<=size(cuts));if(x(i,j)<=cuts(c))exit;c=c+1;end do
      y(i,c)=y(i,c)+1.0_dp
    end do;end do
    status=MIXTOOLS_SUCCESS
  end subroutine makemultdata

  subroutine mixture_information_criteria(loglik,nparam,n,entropy,aic,bic,caic,icl)
    real(dp),intent(in)::loglik,entropy
    integer,intent(in)::nparam,n
    real(dp),intent(out)::aic,bic,caic,icl
    aic=-2.0_dp*loglik+2.0_dp*real(nparam,dp)
    bic=-2.0_dp*loglik+log(real(n,dp))*real(nparam,dp)
    caic=-2.0_dp*loglik+(log(real(n,dp))+1.0_dp)*real(nparam,dp)
    icl=bic+2.0_dp*entropy
  end subroutine mixture_information_criteria

  subroutine normalmix_model_selection(x,components,result,control)
    real(dp),intent(in)::x(:)
    integer,intent(in)::components(:)
    type(model_selection_result),intent(out)::result
    type(em_control),intent(in),optional::control
    type(mixture_result)::fit
    integer::i,k,p,n
    real(dp)::ent
    n=size(x);result%components=components
    allocate(result%aic(size(components)),result%bic(size(components)),result%caic(size(components)))
    allocate(result%icl(size(components)),result%loglik(size(components)))
    do i=1,size(components);k=components(i);call normalmix_em(x,k,fit,control)
      result%loglik(i)=fit%loglik;p=3*k-1;ent=-sum(fit%posterior*log(max(fit%posterior,tiny(1.0_dp))))
      call mixture_information_criteria(fit%loglik,p,n,ent,result%aic(i),result%bic(i),result%caic(i),result%icl(i))
    end do
    result%best_aic=components(minloc(result%aic,dim=1));result%best_bic=components(minloc(result%bic,dim=1));result%status=0
  end subroutine normalmix_model_selection

  subroutine regmix_model_selection(y,x,components,result,control)
    real(dp),intent(in)::y(:),x(:,:);integer,intent(in)::components(:)
    type(model_selection_result),intent(out)::result;type(em_control),intent(in),optional::control
    type(regression_mixture_result)::fit
    integer::i,k,p,n;real(dp)::ent
    n=size(y);result%components=components;allocate(result%aic(size(components)),result%bic(size(components)))
    allocate(result%caic(size(components)),result%icl(size(components)),result%loglik(size(components)))
    do i=1,size(components);k=components(i);call regmix_em(y,x,k,fit,control,.true.)
      result%loglik(i)=fit%loglik;p=k*(size(x,2)+2)+k-1;ent=-sum(fit%posterior*log(max(fit%posterior,tiny(1.0_dp))))
      call mixture_information_criteria(fit%loglik,p,n,ent,result%aic(i),result%bic(i),result%caic(i),result%icl(i))
    end do
    result%best_aic=components(minloc(result%aic,dim=1));result%best_bic=components(minloc(result%bic,dim=1));result%status=0
  end subroutine regmix_model_selection

  subroutine multmix_model_selection(y,components,result,control)
    real(dp),intent(in)::y(:,:);integer,intent(in)::components(:)
    type(model_selection_result),intent(out)::result;type(em_control),intent(in),optional::control
    type(multinomial_mixture_result)::fit
    integer::i,k,p,n;real(dp)::ent
    n=size(y,1);result%components=components;allocate(result%aic(size(components)),result%bic(size(components)))
    allocate(result%caic(size(components)),result%icl(size(components)),result%loglik(size(components)))
    do i=1,size(components);k=components(i);call multmix_em(y,k,fit,control)
      result%loglik(i)=fit%loglik;p=k*(size(y,2)-1)+k-1;ent=-sum(fit%posterior*log(max(fit%posterior,tiny(1.0_dp))))
      call mixture_information_criteria(fit%loglik,p,n,ent,result%aic(i),result%bic(i),result%caic(i),result%icl(i))
    end do
    result%best_aic=components(minloc(result%aic,dim=1));result%best_bic=components(minloc(result%bic,dim=1));result%status=0
  end subroutine multmix_model_selection

  subroutine repnormmix_model_selection(x,components,result,control)
    real(dp),intent(in)::x(:,:);integer,intent(in)::components(:)
    type(model_selection_result),intent(out)::result;type(em_control),intent(in),optional::control
    type(mixture_result)::fit
    integer::i,k,p,n;real(dp)::ent
    n=size(x,2);result%components=components;allocate(result%aic(size(components)),result%bic(size(components)))
    allocate(result%caic(size(components)),result%icl(size(components)),result%loglik(size(components)))
    do i=1,size(components);k=components(i);call repnormmix_em(x,k,fit,control)
      result%loglik(i)=fit%loglik;p=3*k-1;ent=-sum(fit%posterior*log(max(fit%posterior,tiny(1.0_dp))))
      call mixture_information_criteria(fit%loglik,p,n,ent,result%aic(i),result%bic(i),result%caic(i),result%icl(i))
    end do
    result%best_aic=components(minloc(result%aic,dim=1));result%best_bic=components(minloc(result%bic,dim=1));result%status=0
  end subroutine repnormmix_model_selection

  subroutine normalmix_boot_comp(x,k_null,k_alt,b,bootstrap,control,seed)
    real(dp),intent(in)::x(:);integer,intent(in)::k_null,k_alt,b
    type(bootstrap_result),intent(out)::bootstrap;type(em_control),intent(in),optional::control
    integer,intent(in),optional::seed
    type(mixture_result)::nullfit,altfit,nf,af
    type(rng_state)::rng
    real(dp),allocatable::xb(:)
    integer::i,s
    s=12345;if(present(seed))s=seed;call rng_seed(rng,s)
    call normalmix_em(x,k_null,nullfit,control);call normalmix_em(x,k_alt,altfit,control)
    bootstrap%observed=2.0_dp*(altfit%loglik-nullfit%loglik);allocate(bootstrap%statistic(b),xb(size(x)))
    do i=1,b;call rnormmix(rng,size(x),nullfit%lambda,nullfit%mu,nullfit%sigma,xb)
      call normalmix_em(xb,k_null,nf,control);call normalmix_em(xb,k_alt,af,control)
      bootstrap%statistic(i)=2.0_dp*(af%loglik-nf%loglik)
    end do
    bootstrap%successful=b;bootstrap%p_value=(1.0_dp+real(count(bootstrap%statistic>=bootstrap%observed),dp))/real(b+1,dp)
    bootstrap%status=0
  end subroutine normalmix_boot_comp

  subroutine normalmix_boot_se(x,k,b,parameter_se,control,seed)
    real(dp),intent(in)::x(:);integer,intent(in)::k,b
    real(dp),allocatable,intent(out)::parameter_se(:)
    type(em_control),intent(in),optional::control;integer,intent(in),optional::seed
    type(mixture_result)::fit,bfit
    type(rng_state)::rng
    real(dp),allocatable::xb(:),draws(:,:)
    integer::i,s,np
    s=12345;if(present(seed))s=seed;call rng_seed(rng,s);call normalmix_em(x,k,fit,control)
    np=3*k;allocate(draws(b,np),xb(size(x)))
    do i=1,b;call rnormmix(rng,size(x),fit%lambda,fit%mu,fit%sigma,xb);call normalmix_em(xb,k,bfit,control)
      draws(i,:)=[bfit%lambda,bfit%mu,bfit%sigma]
    end do
    allocate(parameter_se(np));parameter_se=sqrt(sum((draws-spread(sum(draws,dim=1)/real(b,dp),1,b))**2,dim=1)/real(max(1,b-1),dp))
  end subroutine normalmix_boot_se

  subroutine likelihood_ratio_pvalue(loglik_null,loglik_alt,df,statistic,p_value)
    real(dp),intent(in)::loglik_null,loglik_alt
    integer,intent(in)::df
    real(dp),intent(out)::statistic,p_value
    statistic=max(0.0_dp,2.0_dp*(loglik_alt-loglik_null))
    ! Wilson-Hilferty normal approximation to the chi-square upper tail.
    p_value=0.5_dp*erfc(((statistic/real(df,dp))**(1.0_dp/3.0_dp) &
      -(1.0_dp-2.0_dp/(9.0_dp*real(df,dp))))/sqrt(4.0_dp/(9.0_dp*real(df,dp))))
  end subroutine likelihood_ratio_pvalue

  subroutine regmix_mh(y,x,k,samples,result,seed,burnin,thin)
    real(dp),intent(in)::y(:),x(:,:);integer,intent(in)::k,samples
    type(mcmc_result),intent(out)::result;integer,intent(in),optional::seed,burnin,thin
    type(regression_mixture_result)::fit
    type(rng_state)::rng
    type(em_control)::ctl
    real(dp),allocatable::xa(:,:),lambda(:),beta(:,:),sigma(:),lw(:),prob(:),alpha(:),zreal(:),w(:)
    integer,allocatable::z(:),counts(:)
    integer::n,p,i,j,it,keep,nb,th,s,status
    real(dp)::ln,shape,scale,rss
    n=size(y);p=size(x,2)+1;nb=100;if(present(burnin))nb=burnin;th=1;if(present(thin))th=thin
    s=12345;if(present(seed))s=seed;call rng_seed(rng,s);ctl%max_iterations=500
    call regmix_em(y,x,k,fit,ctl,.true.);allocate(xa(n,p));xa(:,1)=1.0_dp;xa(:,2:)=x
    lambda=fit%lambda;beta=fit%beta;sigma=fit%sigma;allocate(lw(k),prob(k),alpha(k),z(n),counts(k),zreal(n),w(n))
    allocate(result%lambda_draws(samples,k),result%beta_draws(samples,p,k),result%sigma_draws(samples,k))
    allocate(result%allocation_draws(samples,n));keep=0
    do it=1,nb+samples*th
      counts=0
      do i=1,n;do j=1,k;lw(j)=log(max(lambda(j),tiny(1.0_dp)))+normal_logpdf(y(i),dot_product(xa(i,:),beta(:,j)),sigma(j));end do
        call normalize_logweights(lw,prob,ln);z(i)=draw_category(rng,prob);counts(z(i))=counts(z(i))+1;end do
      alpha=real(counts,dp)+1.0_dp;call random_dirichlet(rng,alpha,lambda)
      do j=1,k
        zreal=merge(1.0_dp,0.0_dp,z==j);w=zreal
        call posterior_beta_draw(rng,xa,y,w,sigma(j),beta(:,j),status)
        rss=sum(w*(y-matmul(xa,beta(:,j)))**2);shape=2.0_dp+0.5_dp*real(counts(j),dp);scale=2.0_dp+0.5_dp*rss
        sigma(j)=sqrt(scale/random_gamma(rng,shape))
      end do
      if(it>nb.and.mod(it-nb,th)==0)then;keep=keep+1;result%lambda_draws(keep,:)=lambda
        result%beta_draws(keep,:,:)=beta;result%sigma_draws(keep,:)=sigma;result%allocation_draws(keep,:)=z;end if
    end do
    result%status=0
  contains
    integer function draw_category(r,pv) result(c)
      type(rng_state),intent(inout)::r;real(dp),intent(in)::pv(:);real(dp)::u,ss;integer::jj
      u=random_uniform(r);ss=0.0_dp;c=size(pv);do jj=1,size(pv);ss=ss+pv(jj);if(u<=ss)then;c=jj;exit;end if;end do
    end function draw_category
    subroutine posterior_beta_draw(r,a,yy,ww,sig,coef,stat)
      type(rng_state),intent(inout)::r;real(dp),intent(in)::a(:,:),yy(:),ww(:),sig
      real(dp),intent(inout)::coef(:);integer,intent(out)::stat
      real(dp),allocatable::prec(:,:),cov(:,:),rhs(:),mean(:),one(:,:)
      integer::ii,jj,pp
      pp=size(a,2);allocate(prec(pp,pp),rhs(pp));prec=0.0_dp;rhs=0.0_dp
      do ii=1,size(yy);do jj=1,pp;rhs(jj)=rhs(jj)+ww(ii)*a(ii,jj)*yy(ii)/sig**2;end do
        prec=prec+ww(ii)*outer(a(ii,:))/sig**2;end do
      do ii=1,pp;prec(ii,ii)=prec(ii,ii)+1.0e-4_dp;end do
      call inverse_spd(prec,cov,stat)
      if(stat/=0)return
      mean=matmul(cov,rhs)
      call rmvnorm(r,1,mean,cov,one,stat)
      if(stat==0)coef=one(1,:)
    end subroutine posterior_beta_draw
    pure function outer(v) result(a)
      real(dp),intent(in)::v(:);real(dp)::a(size(v),size(v));integer::ii,jj
      do ii=1,size(v);do jj=1,size(v);a(ii,jj)=v(ii)*v(jj);end do;end do
    end function outer
  end subroutine regmix_mh

  subroutine posterior_beta_intervals(draws,probability,mean,lower,upper,status)
    real(dp),intent(in)::draws(:,:,:),probability
    real(dp),allocatable,intent(out)::mean(:,:),lower(:,:),upper(:,:)
    integer,intent(out)::status
    real(dp),allocatable::w(:)
    integer::i,j,n
    n=size(draws,1);allocate(mean(size(draws,2),size(draws,3)),lower(size(draws,2),size(draws,3)))
    allocate(upper(size(draws,2),size(draws,3)),w(n));w=1.0_dp
    do i=1,size(draws,2);do j=1,size(draws,3);mean(i,j)=sum(draws(:,i,j))/real(n,dp)
      lower(i,j)=wquantile(draws(:,i,j),w,0.5_dp*(1.0_dp-probability))
      upper(i,j)=wquantile(draws(:,i,j),w,1.0_dp-0.5_dp*(1.0_dp-probability));end do;end do
    status=0
  end subroutine posterior_beta_intervals
end module mixtools_support
