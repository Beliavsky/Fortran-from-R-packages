! SPDX-License-Identifier: MIT
module bayesianou_mi
  use bayesianou_kinds, only : dp, status_ok, status_bad_input
  use bayesianou_types
  use bayesianou_math, only : sample_mean, sample_sd, t_quantile_approx
  use bayesianou_model, only : fit_ou_nested_core => fit_ou_nested
  use bayesianou_diagnostics, only : compute_fit_diagnostics
  implicit none
  private
  public :: fit_ou_nested_mi, rubin_combine, pack_pooled_draws

contains

  subroutine fit_ou_nested_mi(phi_draws,base_input,options,m,keep_fits,result)
    real(dp),intent(in)::phi_draws(:,:,:)
    type(ou_input),intent(in)::base_input
    type(ou_options),intent(in)::options
    integer,intent(in)::m
    logical,intent(in)::keep_fits
    type(ou_mi_result),intent(out)::result
    type(ou_input)::inp
    type(ou_fit_result)::fit
    integer::j,n_use,p
    real(dp),allocatable::qm(:,:),um(:,:),means(:),vars(:)
    if(m<2.or.size(phi_draws,3)<2)then;result%status=status_bad_input;return;end if
    n_use=min(m,size(phi_draws,3));p=pooled_parameter_count(options%n_levels,size(phi_draws,2))
    allocate(qm(n_use,p),um(n_use,p),means(p),vars(p))
    if(keep_fits)allocate(result%fits(n_use))
    do j=1,n_use
      inp=base_input
      if(allocated(inp%y))deallocate(inp%y)
      allocate(inp%y(size(phi_draws,1),size(phi_draws,2)));inp%y=phi_draws(:,:,1+int(real(j-1,dp)*real(size(phi_draws,3)-1,dp)/real(max(1,n_use-1),dp)))
      fit%options=options;fit%options%seed=options%seed+j
      call fit_ou_nested_core(inp,fit%options,fit)
      call compute_fit_diagnostics(inp,fit,[1])
      call pack_pooled_draws(fit,means,vars)
      qm(j,:)=means;um(j,:)=vars
      if(keep_fits)result%fits(j)=fit
    end do
    call rubin_combine(qm,um,result%pooled)
    result%status=status_ok
  end subroutine fit_ou_nested_mi

  integer function pooled_parameter_count(n_levels,s)
    integer,intent(in)::n_levels,s
    if(n_levels==1)then
      pooled_parameter_count=4*s+3
    else
      pooled_parameter_count=7*s+6
    end if
  end function pooled_parameter_count

  subroutine pack_pooled_draws(fit,means,vars)
    type(ou_fit_result),intent(in)::fit
    real(dp),intent(out)::means(:),vars(:)
    integer::s,i
    i=0
    do s=1,size(fit%summary%theta)
      call add(fit%draws%theta(:,s));call add(fit%draws%kappa(:,s));call add(fit%draws%a3(:,s));call add(fit%draws%beta0(:,s))
    end do
    call add(fit%draws%beta1);call add(fit%draws%gamma);call add(fit%draws%nu)
    if(fit%options%n_levels>=2)then
      do s=1,size(fit%summary%theta)
        call add(fit%draws%kappa_p(:,s));call add(fit%draws%mu_const(:,s));call add(fit%draws%sigma_p(:,s))
      end do
      call add(fit%draws%m1);call add(fit%draws%m_v);call add(fit%draws%sigma_phi_meas)
    end if
  contains
    subroutine add(x)
      real(dp),intent(in)::x(:)
      i=i+1;means(i)=sample_mean(x);vars(i)=sample_sd(x)**2
    end subroutine add
  end subroutine pack_pooled_draws

  subroutine rubin_combine(qm,um,result)
    real(dp),intent(in)::qm(:,:),um(:,:)
    type(rubin_result),intent(out)::result
    integer::p,j,m
    real(dp)::r,tcrit,eps
    m=size(qm,1);p=size(qm,2);eps=epsilon(1.0_dp)
    allocate(result%estimate(p),result%total_sd(p),result%within_var(p),result%between_var(p),result%df(p),result%fmi(p),result%lo(p),result%hi(p))
    do j=1,p
      result%estimate(j)=sample_mean(qm(:,j));result%within_var(j)=sample_mean(um(:,j));result%between_var(j)=sample_sd(qm(:,j))**2
      result%total_sd(j)=sqrt(max(0.0_dp,result%within_var(j)+(1.0_dp+1.0_dp/real(m,dp))*result%between_var(j)))
      r=(1.0_dp+1.0_dp/real(m,dp))*result%between_var(j)/max(result%within_var(j),eps)
      result%df(j)=real(m-1,dp)*(1.0_dp+1.0_dp/max(r,eps))**2
      result%fmi(j)=(r+2.0_dp/(result%df(j)+3.0_dp))/(r+1.0_dp)
      tcrit=t_quantile_approx(0.975_dp,result%df(j));result%lo(j)=result%estimate(j)-tcrit*result%total_sd(j);result%hi(j)=result%estimate(j)+tcrit*result%total_sd(j)
    end do
  end subroutine rubin_combine

end module bayesianou_mi
