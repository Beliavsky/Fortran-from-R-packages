! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2018 Marius Hofert, Erik Hintz and Christiane Lemieux
module nvmix_compat
  use nvmix_kinds, only : dp,i8
  use nvmix_types
  use nvmix_core
  use nvmix_distributions
  use nvmix_fitting
  implicit none
  private
  public :: dnvmix,pnvmix,rnvmix,qnvmix,dgnvmix,pgnvmix,rgnvmix
  public :: dNorm,pNorm,rNorm,rNorm_sumconstr,fitNorm
  public :: dStudent,pStudent,rStudent,fitStudent
  public :: dgStudent,pgStudent,rgStudent
  public :: dStudentcopula,pStudentcopula,rStudentcopula
  public :: dgStudentcopula,pgStudentcopula,rgStudentcopula
  public :: fitStudentcopula,fitgStudentcopula
contains
  real(dp) function dnvmix(x,model,control,log_density) result(value)
    real(dp), intent(in) :: x(:)
    type(nvmix_model), intent(in) :: model
    type(integration_control), intent(in), optional :: control
    logical, intent(in), optional :: log_density
    logical :: lg
    lg=.false.; if(present(log_density))lg=log_density
    if(present(control))then; value=nvmix_logpdf(x,model,control); else; value=nvmix_logpdf(x,model); end if
    if(.not.lg)value=exp(value)
  end function
  function pnvmix(lower,upper,model,control) result(value)
    real(dp), intent(in) :: lower(:),upper(:)
    type(nvmix_model), intent(in) :: model
    type(integration_control), intent(in), optional :: control
    type(probability_result) :: value
    if(present(control))then
      value=nvmix_probability(lower,upper,model,control)
    else
      value=nvmix_probability(lower,upper,model)
    end if
  end function
  function rnvmix(n,model,seed) result(value)
    integer, intent(in) :: n
    type(nvmix_model), intent(in) :: model
    integer(i8), intent(in), optional :: seed
    type(sample_result) :: value
    if(present(seed))then; value=nvmix_random_sample(n,model,seed); else; value=nvmix_random_sample(n,model); end if
  end function
  real(dp) function qnvmix(p,model,control) result(value)
    real(dp), intent(in) :: p
    type(nvmix_model), intent(in) :: model
    type(integration_control), intent(in), optional :: control
    if(present(control))then; value=nvmix_quantile(p,model,control); else; value=nvmix_quantile(p,model); end if
  end function
  real(dp) function dgnvmix(x,model,control,log_density) result(value)
    real(dp), intent(in) :: x(:)
    type(nvmix_model), intent(in) :: model
    type(integration_control), intent(in), optional :: control
    logical, intent(in), optional :: log_density
    value=dnvmix(x,model,control,log_density)
  end function
  function pgnvmix(lower,upper,model,control) result(value)
    real(dp), intent(in) :: lower(:),upper(:)
    type(nvmix_model), intent(in) :: model
    type(integration_control), intent(in), optional :: control
    type(probability_result) :: value
    value=pnvmix(lower,upper,model,control)
  end function
  function rgnvmix(n,model,seed) result(value)
    integer, intent(in) :: n
    type(nvmix_model), intent(in) :: model
    integer(i8), intent(in), optional :: seed
    type(sample_result) :: value
    value=rnvmix(n,model,seed)
  end function

  real(dp) function dNorm(x,loc,scale,log_density) result(value)
    real(dp), intent(in) :: x(:),loc(:),scale(:,:)
    logical, intent(in), optional :: log_density
    value=dnorm_mv(x,loc,scale,log_density)
  end function
  function pNorm(lower,upper,loc,scale,control) result(value)
    real(dp), intent(in) :: lower(:),upper(:),loc(:),scale(:,:)
    type(integration_control), intent(in), optional :: control
    type(probability_result) :: value
    value=pnorm_mv(lower,upper,loc,scale,control)
  end function
  function rNorm(n,loc,scale,seed) result(value)
    integer, intent(in) :: n
    real(dp), intent(in) :: loc(:),scale(:,:)
    integer(i8), intent(in), optional :: seed
    type(sample_result) :: value
    value=rnorm_mv(n,loc,scale,seed)
  end function
  function rNorm_sumconstr(n,weights,s,seed) result(value)
    integer, intent(in) :: n
    real(dp), intent(in) :: weights(:),s(:)
    integer(i8), intent(in), optional :: seed
    type(sample_result) :: value
    value=rnorm_sum_constraint(n,weights,s,seed)
  end function
  function fitNorm(x) result(value)
    real(dp), intent(in) :: x(:,:)
    type(fit_result) :: value
    value=fit_norm(x)
  end function

  real(dp) function dStudent(x,df,loc,scale,log_density) result(value)
    real(dp), intent(in) :: x(:),df,loc(:),scale(:,:)
    logical, intent(in), optional :: log_density
    value=dstudent_mv(x,df,loc,scale,log_density)
  end function
  function pStudent(lower,upper,df,loc,scale,control) result(value)
    real(dp), intent(in) :: lower(:),upper(:),df,loc(:),scale(:,:)
    type(integration_control), intent(in), optional :: control
    type(probability_result) :: value
    value=pstudent_mv(lower,upper,df,loc,scale,control)
  end function
  function rStudent(n,df,loc,scale,seed) result(value)
    integer, intent(in) :: n
    real(dp), intent(in) :: df,loc(:),scale(:,:)
    integer(i8), intent(in), optional :: seed
    type(sample_result) :: value
    value=rstudent_mv(n,df,loc,scale,seed)
  end function
  function fitStudent(x,estimate_df,df_fixed) result(value)
    real(dp), intent(in) :: x(:,:)
    logical, intent(in), optional :: estimate_df
    real(dp), intent(in), optional :: df_fixed
    type(fit_result) :: value
    value=fit_student(x,estimate_df,df_fixed)
  end function

  real(dp) function dgStudent(x,groupings,df,loc,scale,control,log_density) result(value)
    real(dp), intent(in) :: x(:),df(:),loc(:),scale(:,:)
    integer, intent(in) :: groupings(:)
    type(integration_control), intent(in), optional :: control
    logical, intent(in), optional :: log_density
    value=dgrouped_student(x,groupings,df,loc,scale,control,log_density)
  end function
  function pgStudent(lower,upper,groupings,df,loc,scale,control) result(value)
    real(dp), intent(in) :: lower(:),upper(:),df(:),loc(:),scale(:,:)
    integer, intent(in) :: groupings(:)
    type(integration_control), intent(in), optional :: control
    type(probability_result) :: value
    value=pgrouped_student(lower,upper,groupings,df,loc,scale,control)
  end function
  function rgStudent(n,groupings,df,loc,scale,seed) result(value)
    integer, intent(in) :: n,groupings(:)
    real(dp), intent(in) :: df(:),loc(:),scale(:,:)
    integer(i8), intent(in), optional :: seed
    type(sample_result) :: value
    value=rgrouped_student(n,groupings,df,loc,scale,seed)
  end function

  real(dp) function dStudentcopula(u,df,scale,log_density) result(value)
    real(dp), intent(in) :: u(:),df,scale(:,:)
    logical, intent(in), optional :: log_density
    value=dstudent_copula(u,df,scale,log_density)
  end function
  function pStudentcopula(lower,upper,df,scale,control) result(value)
    real(dp), intent(in) :: lower(:),upper(:),df,scale(:,:)
    type(integration_control), intent(in), optional :: control
    type(probability_result) :: value
    value=pstudent_copula(lower,upper,df,scale,control)
  end function
  function rStudentcopula(n,df,scale,seed) result(value)
    integer, intent(in) :: n
    real(dp), intent(in) :: df,scale(:,:)
    integer(i8), intent(in), optional :: seed
    type(sample_result) :: value
    value=rstudent_copula(n,df,scale,seed)
  end function
  real(dp) function dgStudentcopula(u,groupings,df,scale,control,log_density) result(value)
    real(dp), intent(in) :: u(:),df(:),scale(:,:)
    integer, intent(in) :: groupings(:)
    type(integration_control), intent(in), optional :: control
    logical, intent(in), optional :: log_density
    value=dgrouped_student_copula(u,groupings,df,scale,control,log_density)
  end function
  function pgStudentcopula(lower,upper,groupings,df,scale,control) result(value)
    real(dp), intent(in) :: lower(:),upper(:),df(:),scale(:,:)
    integer, intent(in) :: groupings(:)
    type(integration_control), intent(in), optional :: control
    type(probability_result) :: value
    value=pgrouped_student_copula(lower,upper,groupings,df,scale,control)
  end function
  function rgStudentcopula(n,groupings,df,scale,seed) result(value)
    integer, intent(in) :: n,groupings(:)
    real(dp), intent(in) :: df(:),scale(:,:)
    integer(i8), intent(in), optional :: seed
    type(sample_result) :: value
    value=rgrouped_student_copula(n,groupings,df,scale,seed)
  end function
  function fitStudentcopula(u) result(value)
    real(dp), intent(in) :: u(:,:)
    type(fit_result) :: value
    value=fit_student_copula(u)
  end function
  function fitgStudentcopula(u,groupings) result(value)
    real(dp), intent(in) :: u(:,:)
    integer, intent(in) :: groupings(:)
    type(fit_result) :: value
    value=fit_grouped_student_copula(u,groupings)
  end function
end module nvmix_compat
