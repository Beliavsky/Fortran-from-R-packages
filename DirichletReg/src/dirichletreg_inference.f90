! SPDX-License-Identifier: GPL-2.0-or-later
module dirichletreg_inference
  use dirichletreg_kinds, only : dp
  use dirichletreg_special, only : normal_cdf, normal_quantile, chi_square_sf
  use dirichletreg_types, only : dirichletreg_model
  implicit none
  private
  public :: standardized_residuals, raw_residuals, composite_residuals
  public :: wald_confint, coefficient_tests, likelihood_ratio_test

contains

  subroutine raw_residuals(y,model,resid,stat)
    real(dp),intent(in)::y(:,:)
    type(dirichletreg_model),intent(in)::model
    real(dp),intent(out)::resid(:,:)
    integer,intent(out),optional::stat
    if(present(stat)) stat=0
    if(any(shape(y)/=shape(model%mu)) .or. any(shape(resid)/=shape(y))) then
      if(present(stat)) stat=1; resid=0.0_dp; return
    end if
    resid=y-model%mu
  end subroutine raw_residuals


  subroutine standardized_residuals(y,model,resid,stat)
    real(dp),intent(in)::y(:,:)
    type(dirichletreg_model),intent(in)::model
    real(dp),intent(out)::resid(:,:)
    integer,intent(out),optional::stat
    integer :: j
    real(dp),allocatable::v(:,:)
    if(present(stat)) stat=0
    if(any(shape(y)/=shape(model%mu)) .or. any(shape(resid)/=shape(y)) .or. size(model%phi)/=size(y,1)) then
      if(present(stat)) stat=1; resid=0.0_dp; return
    end if
    allocate(v(size(y,1),size(y,2)))
    do j=1,size(y,2)
      v(:,j)=model%mu(:,j)*(1.0_dp-model%mu(:,j))/(1.0_dp+model%phi)
    end do
    resid=(y-model%mu)/sqrt(v)
  end subroutine standardized_residuals


  subroutine composite_residuals(y,model,resid,stat)
    real(dp),intent(in)::y(:,:)
    type(dirichletreg_model),intent(in)::model
    real(dp),intent(out)::resid(:)
    integer,intent(out),optional::stat
    real(dp),allocatable::r(:,:)
    integer::i,ierr
    if(size(resid)/=size(y,1)) then
      if(present(stat)) stat=1; resid=0.0_dp; return
    end if
    allocate(r(size(y,1),size(y,2)))
    call standardized_residuals(y,model,r,ierr)
    if(ierr/=0) then
      if(present(stat)) stat=ierr; resid=0.0_dp; return
    end if
    do i=1,size(y,1)
      resid(i)=sum(r(i,:)**2)
    end do
    if(present(stat)) stat=0
  end subroutine composite_residuals


  subroutine wald_confint(model,level,lower,upper,exp_transform,stat)
    type(dirichletreg_model),intent(in)::model
    real(dp),intent(in)::level
    real(dp),intent(out)::lower(:),upper(:)
    logical,intent(in),optional::exp_transform
    integer,intent(out),optional::stat
    real(dp)::z
    logical::doexp
    if(present(stat)) stat=0
    if(level<=0.0_dp .or. level>=1.0_dp .or. size(lower)/=model%npar .or. size(upper)/=model%npar) then
      if(present(stat)) stat=1; lower=0.0_dp; upper=0.0_dp; return
    end if
    z=normal_quantile(0.5_dp+0.5_dp*level)
    lower=model%coefficients-z*model%se
    upper=model%coefficients+z*model%se
    doexp=.false.; if(present(exp_transform)) doexp=exp_transform
    if(doexp) then
      lower=exp(lower); upper=exp(upper)
    end if
  end subroutine wald_confint


  subroutine coefficient_tests(model,zvalue,pvalue,stat)
    type(dirichletreg_model),intent(in)::model
    real(dp),intent(out)::zvalue(:),pvalue(:)
    integer,intent(out),optional::stat
    integer::i
    if(present(stat)) stat=0
    if(size(zvalue)/=model%npar .or. size(pvalue)/=model%npar) then
      if(present(stat)) stat=1; zvalue=0.0_dp; pvalue=0.0_dp; return
    end if
    zvalue=model%coefficients/model%se
    do i=1,model%npar
      pvalue(i)=2.0_dp*normal_cdf(-abs(zvalue(i)))
    end do
  end subroutine coefficient_tests


  subroutine likelihood_ratio_test(loglik_full,npar_full,loglik_reduced,npar_reduced,deviance,df,pvalue,stat)
    real(dp),intent(in)::loglik_full,loglik_reduced
    integer,intent(in)::npar_full,npar_reduced
    real(dp),intent(out)::deviance,df,pvalue
    integer,intent(out),optional::stat
    if(present(stat)) stat=0
    df=real(abs(npar_full-npar_reduced),dp)
    deviance=max(0.0_dp,2.0_dp*(max(loglik_full,loglik_reduced)-min(loglik_full,loglik_reduced)))
    if(df<=0.0_dp) then
      pvalue=1.0_dp
      if(present(stat)) stat=1
    else
      pvalue=chi_square_sf(deviance,df)
    end if
  end subroutine likelihood_ratio_test

end module dirichletreg_inference
