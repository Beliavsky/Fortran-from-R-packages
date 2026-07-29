! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2018 Marius Hofert, Erik Hintz and Christiane Lemieux
module nvmix_risk_dependence
  use nvmix_kinds, only : dp,i8,pi
  use nvmix_types
  use nvmix_core, only : nvmix_quantile,nvmix_random_sample
  use nvmix_special, only : normal_pdf,normal_quantile,student_pdf,student_quantile,student_cdf
  use nvmix_mixing, only : mixing_mean,mixing_mean_sqrt
  implicit none
  private
  public :: var_nvmix,es_nvmix,corgnvmix,kendall_nvmix,spearman_nvmix,lambda_gstudent
contains
  real(dp) function var_nvmix(level,model,control) result(v)
    real(dp), intent(in) :: level
    type(nvmix_model), intent(in) :: model
    type(integration_control), intent(in), optional :: control
    if(present(control))then; v=nvmix_quantile(level,model,control); else; v=nvmix_quantile(level,model); end if
  end function
  real(dp) function es_nvmix(level,model,control,samples,seed) result(v)
    real(dp), intent(in) :: level
    type(nvmix_model), intent(in) :: model
    type(integration_control), intent(in), optional :: control
    integer, intent(in), optional :: samples
    integer(i8), intent(in), optional :: seed
    real(dp) :: q,z,df,sd,s
    type(sample_result) :: draw
    integer :: n,i,k
    if(model%dimension()/=1 .or. level<=0.0_dp .or. level>=1.0_dp)then; v=0.0_dp; return; end if
    sd=sqrt(model%scale(1,1))
    if(model%groups()==1 .and. model%mix_family(1)==mix_constant)then
      z=normal_quantile(level); v=model%loc(1)+sd*normal_pdf(z)/(1.0_dp-level); return
    end if
    if(model%groups()==1 .and. model%mix_family(1)==mix_inverse_gamma)then
      df=model%mix_parameter(1)
      if(df<=1.0_dp)then; v=huge(1.0_dp); return; end if
      z=student_quantile(level,df)
      v=model%loc(1)+sd*student_pdf(z,df)*(df+z*z)/((df-1.0_dp)*(1.0_dp-level)); return
    end if
    if(present(control))then; q=nvmix_quantile(level,model,control); else; q=nvmix_quantile(level,model); end if
    n=100000; if(present(samples))n=samples
    if(present(seed))then; draw=nvmix_random_sample(n,model,seed); else; draw=nvmix_random_sample(n,model,987654321_i8); end if
    s=0.0_dp; k=0
    do i=1,n
      if(draw%x(i,1)>=q)then; s=s+draw%x(i,1); k=k+1; end if
    end do
    if(k>0)then; v=s/real(k,dp); else; v=q; end if
  end function

  function corgnvmix(model) result(correlation)
    type(nvmix_model), intent(in) :: model
    real(dp), allocatable :: correlation(:,:)
    real(dp), allocatable :: covariance(:,:),variance(:),meanw(:),meansqrt(:)
    integer :: d,g,i,j,gi,gj
    d=model%dimension(); g=model%groups(); allocate(correlation(d,d),covariance(d,d),variance(d),meanw(g),meansqrt(g))
    do i=1,g
      meanw(i)=mixing_mean(model%mix_family(i),model%mix_parameter(i))
      meansqrt(i)=mixing_mean_sqrt(model%mix_family(i),model%mix_parameter(i))
    end do
    do i=1,d
      gi=model%groupings(i); variance(i)=model%scale(i,i)*meanw(gi)
      do j=1,d
        gj=model%groupings(j)
        if(gi==gj)then
          covariance(i,j)=model%scale(i,j)*meanw(gi)
        else
          covariance(i,j)=model%scale(i,j)*meansqrt(gi)*meansqrt(gj)
        end if
      end do
    end do
    do i=1,d; do j=1,d
      correlation(i,j)=covariance(i,j)/sqrt(variance(i)*variance(j))
    end do; correlation(i,i)=1.0_dp; end do
  end function

  elemental real(dp) function kendall_nvmix(rho,elliptical) result(tau)
    real(dp), intent(in) :: rho
    logical, intent(in), optional :: elliptical
    if(present(elliptical)) continue
    tau=2.0_dp/pi*asin(max(-1.0_dp,min(1.0_dp,rho)))
  end function
  elemental real(dp) function spearman_nvmix(rho) result(rho_s)
    real(dp), intent(in) :: rho
    rho_s=6.0_dp/pi*asin(max(-1.0_dp,min(1.0_dp,rho))/2.0_dp)
  end function
  real(dp) function lambda_gstudent(df,rho,control) result(lambda)
    real(dp), intent(in) :: df(:),rho
    type(integration_control), intent(in), optional :: control
    real(dp) :: nu
    if(present(control)) continue
    if(size(df)==1 .or. abs(df(1)-df(size(df)))<=1.0e-12_dp)then
      nu=df(1)
      lambda=2.0_dp*student_cdf(-sqrt((nu+1.0_dp)*(1.0_dp-rho)/(1.0_dp+rho)),nu+1.0_dp)
    else
      ! A symmetric approximation for grouped t copulas.  The exact package uses RQMC.
      nu=2.0_dp/(1.0_dp/df(1)+1.0_dp/df(2))
      lambda=2.0_dp*student_cdf(-sqrt((nu+1.0_dp)*(1.0_dp-rho)/(1.0_dp+rho)),nu+1.0_dp)
    end if
  end function
end module nvmix_risk_dependence
