! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2018 Marius Hofert, Erik Hintz and Christiane Lemieux
module nvmix_distributions
  use nvmix_kinds, only : dp,i8,log_two_pi
  use nvmix_types
  use nvmix_core
  use nvmix_special, only : normal_pdf,normal_cdf,normal_quantile,student_pdf,student_cdf,student_quantile
  use nvmix_linalg, only : sample_mean_covariance,cholesky_lower
  use nvmix_random, only : seed_random,random_normal
  implicit none
  private
  public :: dnorm_mv,pnorm_mv,rnorm_mv,rnorm_sum_constraint,fit_norm
  public :: dstudent_mv,pstudent_mv,rstudent_mv
  public :: dgrouped_student,pgrouped_student,rgrouped_student
  public :: dnvmix_copula,pnvmix_copula,rnvmix_copula
  public :: dstudent_copula,pstudent_copula,rstudent_copula
  public :: dgrouped_student_copula,pgrouped_student_copula,rgrouped_student_copula
contains
  real(dp) function dnorm_mv(x,loc,scale,log_density) result(v)
    real(dp), intent(in) :: x(:),loc(:),scale(:,:)
    logical, intent(in), optional :: log_density
    type(nvmix_model) :: m
    logical :: lg
    m=make_nvmix_model(loc,scale,mix_constant,1.0_dp); v=nvmix_logpdf(x,m)
    lg=.false.; if(present(log_density))lg=log_density
    if(.not.lg)v=exp(v)
  end function
  function pnorm_mv(lower,upper,loc,scale,control) result(v)
    real(dp), intent(in) :: lower(:),upper(:),loc(:),scale(:,:)
    type(integration_control), intent(in), optional :: control
    type(probability_result) :: v
    type(nvmix_model) :: m
    m=make_nvmix_model(loc,scale,mix_constant,1.0_dp)
    if(present(control))then; v=nvmix_probability(lower,upper,m,control); else; v=nvmix_probability(lower,upper,m); end if
  end function
  function rnorm_mv(n,loc,scale,seed) result(v)
    integer, intent(in) :: n
    real(dp), intent(in) :: loc(:),scale(:,:)
    integer(i8), intent(in), optional :: seed
    type(sample_result) :: v
    type(nvmix_model) :: m
    m=make_nvmix_model(loc,scale,mix_constant,1.0_dp)
    if(present(seed))then; v=nvmix_random_sample(n,m,seed); else; v=nvmix_random_sample(n,m); end if
  end function
  function rnorm_sum_constraint(n,weights,s,seed) result(v)
    integer, intent(in) :: n
    real(dp), intent(in) :: weights(:),s(:)
    integer(i8), intent(in), optional :: seed
    type(sample_result) :: v
    real(dp), allocatable :: z(:)
    real(dp) :: normw,target,adjustment
    integer :: i,j,d
    d=size(weights); allocate(v%x(max(0,n),d),z(d))
    if(n<1 .or. d<2 .or. any(abs(weights)<=tiny(1.0_dp)) .or. (size(s)/=1 .and. size(s)/=n))then
      v%ok=.false.; v%message='invalid constraint dimensions'; return
    end if
    if(present(seed))call seed_random(seed)
    normw=dot_product(weights,weights)
    do i=1,n
      do j=1,d; z(j)=random_normal(); end do
      if(size(s)==1)then; target=s(1); else; target=s(i); end if
      adjustment=(target-dot_product(weights,z))/normw
      v%x(i,:)=z+adjustment*weights
    end do
  end function
  function fit_norm(x) result(fit)
    real(dp), intent(in) :: x(:,:)
    type(fit_result) :: fit
    logical :: ok
    integer :: n,d,i
    type(nvmix_model) :: m
    call sample_mean_covariance(x,fit%loc,fit%scale,ok)
    if(.not.ok)then; fit%ok=.false.; fit%message='at least two observations are required'; return; end if
    n=size(x,1); d=size(x,2); m=make_nvmix_model(fit%loc,fit%scale,mix_constant,1.0_dp)
    fit%log_likelihood=0.0_dp
    do i=1,n; fit%log_likelihood=fit%log_likelihood+nvmix_logpdf(x(i,:),m); end do
    fit%aic=-2.0_dp*fit%log_likelihood+2.0_dp*real(d+d*(d+1)/2,dp)
    fit%bic=-2.0_dp*fit%log_likelihood+log(real(n,dp))*real(d+d*(d+1)/2,dp)
    fit%converged=.true.
  end function

  real(dp) function dstudent_mv(x,df,loc,scale,log_density) result(v)
    real(dp), intent(in) :: x(:),df,loc(:),scale(:,:)
    logical, intent(in), optional :: log_density
    type(nvmix_model) :: m
    logical :: lg
    m=make_nvmix_model(loc,scale,mix_inverse_gamma,df); v=nvmix_logpdf(x,m)
    lg=.false.; if(present(log_density))lg=log_density
    if(.not.lg)v=exp(v)
  end function
  function pstudent_mv(lower,upper,df,loc,scale,control) result(v)
    real(dp), intent(in) :: lower(:),upper(:),df,loc(:),scale(:,:)
    type(integration_control), intent(in), optional :: control
    type(probability_result) :: v
    type(nvmix_model) :: m
    m=make_nvmix_model(loc,scale,mix_inverse_gamma,df)
    if(present(control))then; v=nvmix_probability(lower,upper,m,control); else; v=nvmix_probability(lower,upper,m); end if
  end function
  function rstudent_mv(n,df,loc,scale,seed) result(v)
    integer, intent(in) :: n
    real(dp), intent(in) :: df,loc(:),scale(:,:)
    integer(i8), intent(in), optional :: seed
    type(sample_result) :: v
    type(nvmix_model) :: m
    m=make_nvmix_model(loc,scale,mix_inverse_gamma,df)
    if(present(seed))then; v=nvmix_random_sample(n,m,seed); else; v=nvmix_random_sample(n,m); end if
  end function

  real(dp) function dgrouped_student(x,groupings,df,loc,scale,control,log_density) result(v)
    real(dp), intent(in) :: x(:),df(:),loc(:),scale(:,:)
    integer, intent(in) :: groupings(:)
    type(integration_control), intent(in), optional :: control
    logical, intent(in), optional :: log_density
    type(nvmix_model) :: m
    integer, allocatable :: families(:)
    logical :: lg
    allocate(families(size(df))); families=mix_inverse_gamma
    m=make_grouped_model(loc,scale,groupings,families,df)
    if(present(control))then; v=nvmix_logpdf(x,m,control); else; v=nvmix_logpdf(x,m); end if
    lg=.false.; if(present(log_density))lg=log_density
    if(.not.lg)v=exp(v)
  end function
  function pgrouped_student(lower,upper,groupings,df,loc,scale,control) result(v)
    real(dp), intent(in) :: lower(:),upper(:),df(:),loc(:),scale(:,:)
    integer, intent(in) :: groupings(:)
    type(integration_control), intent(in), optional :: control
    type(probability_result) :: v
    type(nvmix_model) :: m
    integer, allocatable :: families(:)
    allocate(families(size(df))); families=mix_inverse_gamma
    m=make_grouped_model(loc,scale,groupings,families,df)
    if(present(control))then; v=nvmix_probability(lower,upper,m,control); else; v=nvmix_probability(lower,upper,m); end if
  end function
  function rgrouped_student(n,groupings,df,loc,scale,seed) result(v)
    integer, intent(in) :: n,groupings(:)
    real(dp), intent(in) :: df(:),loc(:),scale(:,:)
    integer(i8), intent(in), optional :: seed
    type(sample_result) :: v
    type(nvmix_model) :: m
    integer, allocatable :: families(:)
    allocate(families(size(df))); families=mix_inverse_gamma
    m=make_grouped_model(loc,scale,groupings,families,df)
    if(present(seed))then; v=nvmix_random_sample(n,m,seed); else; v=nvmix_random_sample(n,m); end if
  end function

  real(dp) function marginal_quantile(u,model,j,control) result(x)
    real(dp), intent(in) :: u
    type(nvmix_model), intent(in) :: model
    integer, intent(in) :: j
    type(integration_control), intent(in), optional :: control
    type(nvmix_model) :: m
    real(dp) :: loc(1),scale(1,1)
    integer :: g
    loc(1)=model%loc(j); scale(1,1)=model%scale(j,j); g=model%groupings(j)
    m=make_nvmix_model(loc,scale,model%mix_family(g),model%mix_parameter(g))
    if(present(control))then; x=nvmix_quantile(u,m,control); else; x=nvmix_quantile(u,m); end if
  end function
  real(dp) function marginal_pdf(x,model,j,control) result(v)
    real(dp), intent(in) :: x
    type(nvmix_model), intent(in) :: model
    integer, intent(in) :: j
    type(integration_control), intent(in), optional :: control
    type(nvmix_model) :: m
    real(dp) :: loc(1),scale(1,1),xx(1)
    integer :: g
    loc(1)=model%loc(j); scale(1,1)=model%scale(j,j); xx(1)=x; g=model%groupings(j)
    m=make_nvmix_model(loc,scale,model%mix_family(g),model%mix_parameter(g))
    if(present(control))then; v=nvmix_pdf(xx,m,control); else; v=nvmix_pdf(xx,m); end if
  end function
  real(dp) function marginal_cdf(x,model,j,control) result(v)
    real(dp), intent(in) :: x
    type(nvmix_model), intent(in) :: model
    integer, intent(in) :: j
    type(integration_control), intent(in), optional :: control
    type(nvmix_model) :: m
    real(dp) :: loc(1),scale(1,1)
    integer :: g
    loc(1)=model%loc(j); scale(1,1)=model%scale(j,j); g=model%groupings(j)
    m=make_nvmix_model(loc,scale,model%mix_family(g),model%mix_parameter(g))
    if(present(control))then; v=nvmix_cdf_1d(x,m,control); else; v=nvmix_cdf_1d(x,m); end if
  end function

  real(dp) function dnvmix_copula(u,model,control,log_density) result(v)
    real(dp), intent(in) :: u(:)
    type(nvmix_model), intent(in) :: model
    type(integration_control), intent(in), optional :: control
    logical, intent(in), optional :: log_density
    real(dp), allocatable :: x(:)
    real(dp) :: ld
    logical :: lg
    integer :: j,d
    d=size(u); allocate(x(d)); lg=.false.; if(present(log_density))lg=log_density
    if(any(u<=0.0_dp) .or. any(u>=1.0_dp))then; if(lg)then; v=-huge(1.0_dp); else; v=0.0_dp; end if; return; end if
    do j=1,d
      if(present(control))then; x(j)=marginal_quantile(u(j),model,j,control); else; x(j)=marginal_quantile(u(j),model,j); end if
    end do
    if(present(control))then; ld=nvmix_logpdf(x,model,control); else; ld=nvmix_logpdf(x,model); end if
    do j=1,d
      if(present(control))then; ld=ld-log(marginal_pdf(x(j),model,j,control)); else; ld=ld-log(marginal_pdf(x(j),model,j)); end if
    end do
    if(lg)then; v=ld; else; v=exp(ld); end if
  end function
  function pnvmix_copula(lower,upper,model,control) result(v)
    real(dp), intent(in) :: lower(:),upper(:)
    type(nvmix_model), intent(in) :: model
    type(integration_control), intent(in), optional :: control
    type(probability_result) :: v
    real(dp), allocatable :: lo(:),hi(:)
    integer :: j,d
    d=size(lower); allocate(lo(d),hi(d))
    do j=1,d
      if(lower(j)<=0.0_dp)then; lo(j)=-huge(1.0_dp); else
        if(present(control))then
          lo(j)=marginal_quantile(lower(j),model,j,control)
        else
          lo(j)=marginal_quantile(lower(j),model,j)
        end if
      end if
      if(upper(j)>=1.0_dp)then; hi(j)=huge(1.0_dp); else
        if(present(control))then
          hi(j)=marginal_quantile(upper(j),model,j,control)
        else
          hi(j)=marginal_quantile(upper(j),model,j)
        end if
      end if
    end do
    if(present(control))then; v=nvmix_probability(lo,hi,model,control); else; v=nvmix_probability(lo,hi,model); end if
  end function
  function rnvmix_copula(n,model,seed,control) result(v)
    integer, intent(in) :: n
    type(nvmix_model), intent(in) :: model
    integer(i8), intent(in), optional :: seed
    type(integration_control), intent(in), optional :: control
    type(sample_result) :: v
    type(sample_result) :: x
    integer :: i,j
    if(present(seed))then; x=nvmix_random_sample(n,model,seed); else; x=nvmix_random_sample(n,model); end if
    allocate(v%x(n,model%dimension())); v%ok=x%ok; v%message=x%message; if(.not.x%ok)return
    do i=1,n; do j=1,model%dimension()
      if(present(control))then
        v%x(i,j)=marginal_cdf(x%x(i,j),model,j,control)
      else
        v%x(i,j)=marginal_cdf(x%x(i,j),model,j)
      end if
    end do; end do
  end function

  real(dp) function dstudent_copula(u,df,scale,log_density) result(v)
    real(dp), intent(in) :: u(:),df,scale(:,:)
    logical, intent(in), optional :: log_density
    type(nvmix_model) :: m
    real(dp), allocatable :: loc(:)
    allocate(loc(size(u))); loc=0.0_dp; m=make_nvmix_model(loc,scale,mix_inverse_gamma,df)
    v=dnvmix_copula(u,m,log_density=log_density)
  end function
  function pstudent_copula(lower,upper,df,scale,control) result(v)
    real(dp), intent(in) :: lower(:),upper(:),df,scale(:,:)
    type(integration_control), intent(in), optional :: control
    type(probability_result) :: v
    type(nvmix_model) :: m
    real(dp), allocatable :: loc(:)
    allocate(loc(size(lower))); loc=0.0_dp; m=make_nvmix_model(loc,scale,mix_inverse_gamma,df)
    if(present(control))then; v=pnvmix_copula(lower,upper,m,control); else; v=pnvmix_copula(lower,upper,m); end if
  end function
  function rstudent_copula(n,df,scale,seed) result(v)
    integer, intent(in) :: n
    real(dp), intent(in) :: df,scale(:,:)
    integer(i8), intent(in), optional :: seed
    type(sample_result) :: v
    type(nvmix_model) :: m
    real(dp), allocatable :: loc(:)
    allocate(loc(size(scale,1))); loc=0.0_dp; m=make_nvmix_model(loc,scale,mix_inverse_gamma,df)
    if(present(seed))then; v=rnvmix_copula(n,m,seed); else; v=rnvmix_copula(n,m); end if
  end function
  real(dp) function dgrouped_student_copula(u,groupings,df,scale,control,log_density) result(v)
    real(dp), intent(in) :: u(:),df(:),scale(:,:)
    integer, intent(in) :: groupings(:)
    type(integration_control), intent(in), optional :: control
    logical, intent(in), optional :: log_density
    type(nvmix_model) :: m
    real(dp), allocatable :: loc(:)
    integer, allocatable :: fam(:)
    allocate(loc(size(u)),fam(size(df))); loc=0.0_dp; fam=mix_inverse_gamma
    m=make_grouped_model(loc,scale,groupings,fam,df)
    if(present(control))then; v=dnvmix_copula(u,m,control,log_density); else; v=dnvmix_copula(u,m,log_density=log_density); end if
  end function
  function pgrouped_student_copula(lower,upper,groupings,df,scale,control) result(v)
    real(dp), intent(in) :: lower(:),upper(:),df(:),scale(:,:)
    integer, intent(in) :: groupings(:)
    type(integration_control), intent(in), optional :: control
    type(probability_result) :: v
    type(nvmix_model) :: m
    real(dp), allocatable :: loc(:)
    integer, allocatable :: fam(:)
    allocate(loc(size(lower)),fam(size(df))); loc=0.0_dp; fam=mix_inverse_gamma
    m=make_grouped_model(loc,scale,groupings,fam,df)
    if(present(control))then; v=pnvmix_copula(lower,upper,m,control); else; v=pnvmix_copula(lower,upper,m); end if
  end function
  function rgrouped_student_copula(n,groupings,df,scale,seed) result(v)
    integer, intent(in) :: n,groupings(:)
    real(dp), intent(in) :: df(:),scale(:,:)
    integer(i8), intent(in), optional :: seed
    type(sample_result) :: v
    type(nvmix_model) :: m
    real(dp), allocatable :: loc(:)
    integer, allocatable :: fam(:)
    allocate(loc(size(scale,1)),fam(size(df))); loc=0.0_dp; fam=mix_inverse_gamma
    m=make_grouped_model(loc,scale,groupings,fam,df)
    if(present(seed))then; v=rnvmix_copula(n,m,seed); else; v=rnvmix_copula(n,m); end if
  end function
end module nvmix_distributions
