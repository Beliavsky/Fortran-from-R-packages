! SPDX-License-Identifier: MIT
module bayesianou_diagnostics
  use bayesianou_kinds, only : dp, status_ok, status_bad_input
  use bayesianou_types
  use bayesianou_math, only : split_rhat, effective_sample_size, percentile, sample_mean, sample_sd, &
                              sort_real, log_sum_exp, clamp
  use bayesianou_model, only : ou_mean_increment
  implicit none
  private
  public :: compute_fit_diagnostics, evaluate_oos, evaluate_oos_nested
  public :: psis_loo, compare_models_loo, kappa_stability_evidence
  public :: extract_convergence_evidence, count_divergences, validate_ou_fit

contains

  subroutine compute_fit_diagnostics(input,fit,horizons)
    type(ou_input),intent(in)::input
    type(ou_fit_result),intent(inout)::fit
    integer,intent(in),optional::horizons(:)
    integer,allocatable::hh(:)
    real(dp),allocatable::loglik(:,:)
    if(present(horizons))then;allocate(hh(size(horizons)));hh=horizons;else;allocate(hh(3));hh=[1,4,8];end if
    call compute_chain_diagnostics(fit)
    call build_loglik_matrix(input,fit,loglik)
    call psis_loo(loglik,fit%diagnostics%loo)
    if(fit%options%n_levels==1)then
      call evaluate_oos(input,fit,hh,fit%diagnostics%oos)
    else
      call evaluate_oos_nested(input,fit,hh,fit%diagnostics%oos)
    end if
  end subroutine compute_fit_diagnostics

  subroutine compute_chain_diagnostics(fit)
    type(ou_fit_result),intent(inout)::fit
    integer::s,p,nc,nper,i,offset
    real(dp),allocatable::chains(:,:)
    p=4*size(fit%summary%theta)+3
    allocate(fit%diagnostics%rhat(p),fit%diagnostics%ess(p));i=0
    nc=maxval(fit%draws%chain);nper=count(fit%draws%chain==1)
    allocate(chains(nper,nc))
    do s=1,size(fit%summary%theta)
      call diag_one(fit%draws%theta(:,s));call diag_one(fit%draws%kappa(:,s));call diag_one(fit%draws%a3(:,s));call diag_one(fit%draws%beta0(:,s))
    end do
    call diag_one(fit%draws%beta1);call diag_one(fit%draws%gamma);call diag_one(fit%draws%nu)
    fit%diagnostics%rhat_max=maxval(fit%diagnostics%rhat)
    fit%diagnostics%rhat_share=real(count(fit%diagnostics%rhat>1.01_dp),dp)/real(p,dp)
    if(fit%diagnostics%acceptance_rate<=0.0_dp)fit%diagnostics%acceptance_rate=1.0_dp
  contains
    subroutine diag_one(x)
      real(dp),intent(in)::x(:)
      integer::cc,jj
      i=i+1
      do cc=1,nc
        jj=0
        do offset=1,size(x)
          if(fit%draws%chain(offset)==cc)then;jj=jj+1;chains(jj,cc)=x(offset);end if
        end do
      end do
      fit%diagnostics%rhat(i)=split_rhat(chains)
      fit%diagnostics%ess(i)=0.0_dp
      do cc=1,nc;fit%diagnostics%ess(i)=fit%diagnostics%ess(i)+effective_sample_size(chains(:,cc));end do
    end subroutine diag_one
  end subroutine compute_chain_diagnostics

  subroutine build_loglik_matrix(input,fit,ll)
    type(ou_input),intent(in)::input
    type(ou_fit_result),intent(in)::fit
    real(dp),allocatable,intent(out)::ll(:,:)
    type(ou_summary)::s
    integer::d,t,j,k,nobs
    real(dp)::resid,sd
    nobs=(fit%t_lik-1)*size(input%y,2);allocate(ll(size(fit%draws%beta1),nobs))
    do d=1,size(ll,1)
      call draw_to_summary(fit,d,s);k=0
      do t=2,fit%t_lik;do j=1,size(input%y,2);k=k+1
        resid=fit%zy%mz(t,j)-fit%zy%mz(t-1,j)-ou_mean_increment(input,s,fit,t,j)
        sd=clamp(exp(0.5_dp*fit%h_median(t,j)),1e-8_dp,1e8_dp)
        if(fit%options%n_levels>=2.and..not.fit%options%level_spec%level1%student_t)then
          ll(d,k)=-0.5_dp*log(2.0_dp*acos(-1.0_dp))-log(sd)-0.5_dp*(resid/sd)**2
        else
          ll(d,k)=log_gamma(0.5_dp*(s%nu+1))-log_gamma(0.5_dp*s%nu)-0.5_dp*log(s%nu*acos(-1.0_dp))-log(sd) &
                  -0.5_dp*(s%nu+1)*log(1+(resid/sd)**2/s%nu)
        end if
      end do;end do
    end do
  end subroutine build_loglik_matrix

  subroutine draw_to_summary(fit,d,s)
    type(ou_fit_result),intent(in)::fit
    integer,intent(in)::d
    type(ou_summary),intent(out)::s
    integer::n
    n=size(fit%summary%theta)
    allocate(s%theta(n),s%kappa(n),s%a3(n),s%beta0(n),s%alpha(n),s%rho(n),s%sigma_eta(n))
    allocate(s%kappa_p(n),s%mu_const(n),s%sigma_p(n),s%a3_p(n))
    s%theta=fit%draws%theta(d,:);s%kappa=fit%draws%kappa(d,:);s%a3=fit%draws%a3(d,:);s%beta0=fit%draws%beta0(d,:)
    s%alpha=fit%draws%alpha(d,:);s%rho=fit%draws%rho(d,:);s%sigma_eta=fit%draws%sigma_eta(d,:)
    s%kappa_p=fit%draws%kappa_p(d,:);s%mu_const=fit%draws%mu_const(d,:);s%sigma_p=fit%draws%sigma_p(d,:);s%a3_p=fit%draws%a3_p(d,:)
    s%beta1=fit%draws%beta1(d);s%gamma=fit%draws%gamma(d);s%nu=fit%draws%nu(d);s%m1=fit%draws%m1(d);s%m_v=fit%draws%m_v(d)
    s%sigma_phi_meas=fit%draws%sigma_phi_meas(d)
  end subroutine draw_to_summary

  subroutine psis_loo(log_lik,result)
    real(dp),intent(in)::log_lik(:,:)
    type(loo_result),intent(out)::result
    integer::i,n,m,ntail
    real(dp),allocatable::lw(:),w(:),ws(:),elpd(:),lpd(:)
    real(dp)::mx,threshold,mean_exc,var_exc,k,sigma,cap,den,num
    n=size(log_lik,1);m=size(log_lik,2);allocate(result%pointwise(m),result%pareto_k(m),elpd(m),lpd(m),lw(n),w(n),ws(n))
    ntail=max(3,min(n/5,int(3.0_dp*sqrt(real(n,dp)))))
    do i=1,m
      lw=-log_lik(:,i);mx=maxval(lw);w=exp(lw-mx);ws=w;call sort_real(ws)
      threshold=ws(max(1,n-ntail));mean_exc=sample_mean(ws(n-ntail+1:n)-threshold);var_exc=sample_sd(ws(n-ntail+1:n)-threshold)**2
      if(var_exc>mean_exc**2.and.mean_exc>0)then
        k=0.5_dp*(1.0_dp-mean_exc**2/var_exc);sigma=0.5_dp*mean_exc*(1+mean_exc**2/var_exc)
      else;k=0.0_dp;sigma=max(mean_exc,1e-12_dp);end if
      result%pareto_k(i)=clamp(k,-0.5_dp,5.0_dp)
      cap=sample_mean(w)*real(n,dp)**0.75_dp;w=min(w,cap);den=sum(w)
      num=sum(w*exp(log_lik(:,i)-maxval(log_lik(:,i))))
      elpd(i)=maxval(log_lik(:,i))+log(max(num,1e-300_dp)/max(den,1e-300_dp))
      lpd(i)=log_sum_exp(log_lik(:,i))-log(real(n,dp));result%pointwise(i)=elpd(i)
    end do
    result%elpd_loo=sum(elpd);result%se_elpd=sqrt(real(m,dp))*sample_sd(elpd)
    result%p_loo=sum(lpd-elpd);result%looic=-2.0_dp*result%elpd_loo
  end subroutine psis_loo

  subroutine evaluate_oos(input,fit,horizons,metrics)
    type(ou_input),intent(in)::input
    type(ou_fit_result),intent(in)::fit
    integer,intent(in)::horizons(:)
    type(oos_metric),allocatable,intent(out)::metrics(:)
    integer::ih,h,t,step,s,last,nerr
    real(dp),allocatable::pred(:),err(:)
    real(dp)::dev,mu,comstd
    allocate(metrics(size(horizons)))
    do ih=1,size(horizons);h=horizons(ih);metrics(ih)%horizon=h;last=size(input%y,1)-h+1
      if(last<fit%t_train+1)cycle
      allocate(err((last-fit%t_train)*size(input%y,2)));nerr=0
      do t=fit%t_train+1,last
        allocate(pred(size(input%y,2)));pred=fit%zy%mz(t-1,:)
        do step=1,h
          do s=1,size(pred)
            dev=pred(s)-fit%summary%theta(s)
            mu=fit%summary%kappa(s)*(-dev+fit%summary%a3(s)*dev**3)+(fit%summary%beta0(s)+fit%summary%beta1*fit%ztmg(t-1+step))*fit%zx%mz(t-2+step,s)
            if(fit%options%com_in_mean)then;comstd=(input%com(t-2+step,s)-fit%com_wmean(s))/fit%com_wsd(s);mu=mu+fit%summary%gamma*comstd;end if
            pred(s)=pred(s)+mu
          end do
        end do
        do s=1,size(pred);nerr=nerr+1;err(nerr)=fit%zy%mz(t+h-1,s)-pred(s);end do
        deallocate(pred)
      end do
      metrics(ih)%n_obs=nerr;metrics(ih)%rmse=sqrt(sum(err(:nerr)**2)/real(nerr,dp));metrics(ih)%mae=sum(abs(err(:nerr)))/real(nerr,dp);deallocate(err)
    end do
  end subroutine evaluate_oos

  subroutine evaluate_oos_nested(input,fit,horizons,metrics)
    type(ou_input),intent(in)::input
    type(ou_fit_result),intent(in)::fit
    integer,intent(in)::horizons(:)
    type(oos_metric),allocatable,intent(out)::metrics(:)
    integer::ih,h,t,step,s,last,nerr
    real(dp),allocatable::pred(:),err(:)
    real(dp)::dev,mu,k,comstd
    allocate(metrics(size(horizons)))
    do ih=1,size(horizons);h=horizons(ih);metrics(ih)%horizon=h;last=size(input%y,1)-h+1
      if(last<fit%t_train+1)cycle
      allocate(err((last-fit%t_train)*size(input%y,2)));nerr=0
      do t=fit%t_train+1,last
        allocate(pred(size(input%y,2)));pred=fit%zy%mz(t-1,:)
        do step=1,h
          do s=1,size(pred)
            dev=pred(s)-fit%phi_median(t-2+step,s)
            k=fit%options%kappa_cap/(1+exp(-(log(fit%summary%kappa(s)/(fit%options%kappa_cap-fit%summary%kappa(s)))+fit%summary%beta1*fit%ztmg(t-1+step))))
            mu=k*(-dev+fit%summary%a3(s)*dev**3)
            if(fit%options%com_in_mean)then;comstd=(input%com(t-2+step,s)-fit%com_wmean(s))/fit%com_wsd(s);mu=mu+fit%summary%gamma*comstd;end if
            pred(s)=pred(s)+mu
          end do
        end do
        do s=1,size(pred);nerr=nerr+1;err(nerr)=fit%zy%mz(t+h-1,s)-pred(s);end do
        deallocate(pred)
      end do
      metrics(ih)%n_obs=nerr;metrics(ih)%rmse=sqrt(sum(err(:nerr)**2)/real(nerr,dp));metrics(ih)%mae=sum(abs(err(:nerr)))/real(nerr,dp);deallocate(err)
    end do
  end subroutine evaluate_oos_nested

  subroutine kappa_stability_evidence(fit,result)
    type(ou_fit_result),intent(in)::fit
    type(stability_result),intent(out)::result
    integer::s,d
    logical::allstable
    allocate(result%intervals(size(fit%summary%kappa),3))
    do s=1,size(fit%summary%kappa)
      result%intervals(s,:)=[percentile(fit%draws%kappa(:,s),0.025_dp),percentile(fit%draws%kappa(:,s),0.5_dp),percentile(fit%draws%kappa(:,s),0.975_dp)]
    end do
    result%stable=all(result%intervals(:,1)>0.and.result%intervals(:,3)<1)
    result%probability=0
    do d=1,size(fit%draws%kappa,1);allstable=all(fit%draws%kappa(d,:)>0.and.fit%draws%kappa(d,:)<1);if(allstable)result%probability=result%probability+1;end do
    result%probability=result%probability/real(size(fit%draws%kappa,1),dp)
  end subroutine kappa_stability_evidence

  subroutine extract_convergence_evidence(fit,result)
    type(ou_fit_result),intent(in)::fit
    type(stability_result),intent(out)::result
    call kappa_stability_evidence(fit,result)
  end subroutine extract_convergence_evidence

  integer function count_divergences(fit)
    type(ou_fit_result),intent(in)::fit
    count_divergences=fit%diagnostics%divergences
  end function count_divergences

  subroutine compare_models_loo(new_fit,base_fit,delta_elpd,se_delta,status)
    type(ou_fit_result),intent(in)::new_fit,base_fit
    real(dp),intent(out)::delta_elpd,se_delta
    integer,intent(out)::status
    real(dp),allocatable::d(:)
    if(size(new_fit%diagnostics%loo%pointwise)/=size(base_fit%diagnostics%loo%pointwise))then;status=status_bad_input;delta_elpd=0;se_delta=0;return;end if
    d=new_fit%diagnostics%loo%pointwise-base_fit%diagnostics%loo%pointwise
    delta_elpd=sum(d);se_delta=sqrt(real(size(d),dp))*sample_sd(d);status=status_ok
  end subroutine compare_models_loo

  logical function validate_ou_fit(fit)
    type(ou_fit_result),intent(in)::fit
    validate_ou_fit=fit%status==status_ok.and.allocated(fit%diagnostics%rhat).and.allocated(fit%diagnostics%loo%pointwise)
  end function validate_ou_fit

end module bayesianou_diagnostics
