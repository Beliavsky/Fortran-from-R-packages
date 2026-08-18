! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2012 Marius Hofert, Ivan Kojadinovic, Martin Maechler and Jun Yan
module copula_dependence
  use copula_kinds, only : dp, pi
  use copula_types
  use copula_special, only : debye1, student_cdf
  use copula_families, only : copula_cdf, copula_density
  implicit none
  private
  public :: kendall_tau_model, spearman_rho_model, tail_dependence
  public :: parameter_from_tau, parameter_from_rho
contains
  real(dp) function kendall_tau_model(model) result(value)
    type(copula_model), intent(in) :: model
    type(copula_model) :: base
    real(dp) :: sign_value
    base = model
    base%rotation = rotation_none
    select case (base%family)
    case (family_independence)
      value = 0.0_dp
    case (family_gaussian, family_student)
      value = 2.0_dp*asin(base%correlation(1,2))/pi
    case (family_clayton)
      value = base%theta/(base%theta+2.0_dp)
    case (family_gumbel)
      value = 1.0_dp-1.0_dp/base%theta
    case (family_frank)
      value = 1.0_dp-4.0_dp/base%theta+4.0_dp*debye1(base%theta)/base%theta
    case (family_fgm)
      value = 2.0_dp*base%theta/9.0_dp
    case (family_marshall_olkin)
      value = base%alpha1*base%alpha2/(base%alpha1+base%alpha2-base%alpha1*base%alpha2)
    case (family_lower_fh)
      value = -1.0_dp
    case (family_upper_fh)
      value = 1.0_dp
    case default
      value = numerical_tau(base)
    end select
    sign_value = 1.0_dp
    if (model%rotation == rotation_90 .or. model%rotation == rotation_270) sign_value = -1.0_dp
    value = sign_value*value
  end function kendall_tau_model

  real(dp) function spearman_rho_model(model) result(value)
    type(copula_model), intent(in) :: model
    type(copula_model) :: base
    real(dp) :: sign_value
    base = model
    base%rotation = rotation_none
    select case (base%family)
    case (family_independence)
      value = 0.0_dp
    case (family_gaussian)
      value = 6.0_dp*asin(base%correlation(1,2)/2.0_dp)/pi
    case (family_fgm)
      value = base%theta/3.0_dp
    case (family_lower_fh)
      value = -1.0_dp
    case (family_upper_fh)
      value = 1.0_dp
    case default
      value = numerical_rho(base)
    end select
    sign_value = 1.0_dp
    if (model%rotation == rotation_90 .or. model%rotation == rotation_270) sign_value = -1.0_dp
    value = sign_value*value
  end function spearman_rho_model

  function tail_dependence(model) result(lambda)
    type(copula_model), intent(in) :: model
    real(dp) :: lambda(2)
    type(copula_model) :: base
    real(dp) :: rho, df, lower, upper
    base = model
    base%rotation = rotation_none
    lower = 0.0_dp
    upper = 0.0_dp
    select case (base%family)
    case (family_student)
      rho = base%correlation(1,2)
      df = base%df
      lower = 2.0_dp*student_cdf(-sqrt((df+1.0_dp)*(1.0_dp-rho)/(1.0_dp+rho)),df+1.0_dp)
      upper = lower
    case (family_clayton)
      lower = 2.0_dp**(-1.0_dp/base%theta)
    case (family_gumbel, family_joe)
      upper = 2.0_dp-2.0_dp**(1.0_dp/base%theta)
    case (family_upper_fh)
      lower = 1.0_dp
      upper = 1.0_dp
    case (family_lower_fh)
      lower = 0.0_dp
      upper = 0.0_dp
    case (family_marshall_olkin)
      lower = 0.0_dp
      upper = min(base%alpha1,base%alpha2)
    case (family_galambos, family_husler_reiss, family_tawn)
      upper = 2.0_dp-2.0_dp*copula_cdf([0.5_dp,0.5_dp],base)
    end select
    select case (model%rotation)
    case (rotation_none)
      lambda = [lower,upper]
    case (rotation_180)
      lambda = [upper,lower]
    case default
      lambda = [0.0_dp,0.0_dp]
    end select
  end function tail_dependence

  real(dp) function parameter_from_tau(family, tau_value) result(theta)
    integer, intent(in) :: family
    real(dp), intent(in) :: tau_value
    type(copula_model) :: model
    real(dp) :: lo, hi, mid, target
    integer :: iteration
    target = tau_value
    select case (family)
    case (family_gaussian, family_student)
      theta = sin(0.5_dp*pi*target)
      return
    case (family_clayton)
      theta = 2.0_dp*target/max(1.0_dp-target,tiny(1.0_dp))
      return
    case (family_gumbel)
      theta = 1.0_dp/max(1.0_dp-target,tiny(1.0_dp))
      return
    case (family_fgm)
      theta = 4.5_dp*target
      return
    case (family_frank)
      lo = -50.0_dp
      hi = 50.0_dp
      if (target >= 0.0_dp) lo = 1.0e-7_dp
      if (target < 0.0_dp) hi = -1.0e-7_dp
    case (family_amh)
      lo = -0.999_dp
      hi = 0.999_dp
    case (family_joe)
      lo = 1.0_dp
      hi = 50.0_dp
    case (family_plackett)
      lo = 1.0e-4_dp
      hi = 1.0e4_dp
    case (family_galambos, family_husler_reiss)
      lo = 1.0e-3_dp
      hi = 50.0_dp
    case default
      theta = 0.0_dp
      return
    end select
    model%family = family
    model%dimension = 2
    do iteration = 1, 80
      mid = 0.5_dp*(lo+hi)
      model%theta = mid
      if (kendall_tau_model(model) < target) then
        lo = mid
      else
        hi = mid
      end if
    end do
    theta = 0.5_dp*(lo+hi)
  end function parameter_from_tau

  real(dp) function parameter_from_rho(family, rho_value) result(theta)
    integer, intent(in) :: family
    real(dp), intent(in) :: rho_value
    type(copula_model) :: model
    real(dp) :: lo, hi, mid
    integer :: iteration
    select case (family)
    case (family_gaussian)
      theta = 2.0_dp*sin(pi*rho_value/6.0_dp)
      return
    case (family_fgm)
      theta = 3.0_dp*rho_value
      return
    case (family_clayton)
      lo = 1.0e-5_dp
      hi = 100.0_dp
    case (family_gumbel, family_joe)
      lo = 1.0_dp
      hi = 50.0_dp
    case (family_frank)
      lo = -50.0_dp
      hi = 50.0_dp
    case (family_plackett)
      lo = 1.0e-4_dp
      hi = 1.0e4_dp
    case default
      theta = 0.0_dp
      return
    end select
    model%family = family
    model%dimension = 2
    do iteration = 1, 80
      mid = 0.5_dp*(lo+hi)
      if (abs(mid) < 1.0e-8_dp .and. family == family_frank) mid = 1.0e-8_dp
      model%theta = mid
      if (spearman_rho_model(model) < rho_value) then
        lo = mid
      else
        hi = mid
      end if
    end do
    theta = 0.5_dp*(lo+hi)
  end function parameter_from_rho

  real(dp) function numerical_rho(model) result(value)
    type(copula_model), intent(in) :: model
    integer, parameter :: n = 100
    integer :: i, j
    real(dp) :: u, v, total
    total = 0.0_dp
    do i = 1, n
      u = (real(i,dp)-0.5_dp)/real(n,dp)
      do j = 1, n
        v = (real(j,dp)-0.5_dp)/real(n,dp)
        total = total+copula_cdf([u,v],model)-u*v
      end do
    end do
    value = 12.0_dp*total/real(n*n,dp)
  end function numerical_rho

  real(dp) function numerical_tau(model) result(value)
    type(copula_model), intent(in) :: model
    integer, parameter :: n = 70
    integer :: i, j
    real(dp) :: u, v, total
    total = 0.0_dp
    do i = 1, n
      u = (real(i,dp)-0.5_dp)/real(n,dp)
      do j = 1, n
        v = (real(j,dp)-0.5_dp)/real(n,dp)
        total = total+copula_cdf([u,v],model)*copula_density([u,v],model)
      end do
    end do
    value = 4.0_dp*total/real(n*n,dp)-1.0_dp
    value = min(1.0_dp,max(-1.0_dp,value))
  end function numerical_tau
end module copula_dependence
