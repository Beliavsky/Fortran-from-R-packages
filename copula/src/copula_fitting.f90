! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2012 Marius Hofert, Ivan Kojadinovic, Martin Maechler and Jun Yan
module copula_fitting
  use copula_kinds, only : dp
  use copula_types
  use copula_families, only : copula_log_density, conditional_cdf
  use copula_dependence, only : parameter_from_tau
  use copula_empirical, only : sample_kendall_tau
  implicit none
  private
  public :: fit_copula, copula_log_likelihood, rosenblatt_transform, inverse_rosenblatt
contains
  function fit_copula(u, family, method, df) result(result)
    real(dp), intent(in) :: u(:,:)
    integer, intent(in) :: family
    character(len=*), intent(in), optional :: method
    real(dp), intent(in), optional :: df
    type(fit_result) :: result
    character(len=16) :: selected_method
    real(dp) :: tau_value, theta, lo, hi, x1, x2, f1, f2
    integer :: iteration
    result%model%family = family
    result%model%dimension = size(u,2)
    if (present(df)) result%model%df = df
    result%ok = size(u,1) >= 3 .and. size(u,2) == 2 .and. all(u > 0.0_dp) .and. all(u < 1.0_dp)
    if (.not. result%ok) then
      result%message = 'fitting requires at least three bivariate observations in (0,1)'
      return
    end if
    selected_method = 'mpl'
    if (present(method)) selected_method = lower(method)
    tau_value = sample_kendall_tau(u(:,1),u(:,2))
    theta = parameter_from_tau(family,tau_value)
    if (family == family_gaussian .or. family == family_student) then
      allocate(result%model%correlation(2,2))
      result%model%correlation = reshape([1.0_dp,theta,theta,1.0_dp],[2,2])
    else
      result%model%theta = theta
    end if
    if (trim(selected_method) == 'itau') then
      result%log_likelihood = copula_log_likelihood(u,result%model)
      result%converged = .true.
      result%iterations = 0
      result%standard_error = numerical_standard_error(u,result%model)
      return
    end if
    call parameter_bounds(family,lo,hi,result%ok)
    if (.not. result%ok) then
      result%message = 'family is not supported by the one-parameter fitter'
      return
    end if
    x1 = hi-0.6180339887498949_dp*(hi-lo)
    x2 = lo+0.6180339887498949_dp*(hi-lo)
    f1 = objective(x1)
    f2 = objective(x2)
    do iteration = 1, 200
      if (f1 > f2) then
        lo = x1
        x1 = x2
        f1 = f2
        x2 = lo+0.6180339887498949_dp*(hi-lo)
        f2 = objective(x2)
      else
        hi = x2
        x2 = x1
        f2 = f1
        x1 = hi-0.6180339887498949_dp*(hi-lo)
        f1 = objective(x1)
      end if
      if (abs(hi-lo) < 1.0e-8_dp*(1.0_dp+abs(0.5_dp*(hi+lo)))) exit
    end do
    theta = 0.5_dp*(lo+hi)
    call set_parameter(result%model,theta)
    result%log_likelihood = copula_log_likelihood(u,result%model)
    result%iterations = iteration
    result%converged = iteration < 200
    result%standard_error = numerical_standard_error(u,result%model)
    result%ok = result%log_likelihood > -0.5_dp*huge(1.0_dp)
  contains
    real(dp) function objective(parameter) result(value)
      real(dp), intent(in) :: parameter
      type(copula_model) :: trial
      trial = result%model
      call set_parameter(trial,parameter)
      if (.not. trial%valid()) then
        value = huge(1.0_dp)
      else
        value = -copula_log_likelihood(u,trial)
      end if
    end function objective
  end function fit_copula

  real(dp) function copula_log_likelihood(u, model) result(value)
    real(dp), intent(in) :: u(:,:)
    type(copula_model), intent(in) :: model
    real(dp) :: term
    integer :: i
    if (size(u,2) /= model%dimension .or. any(u <= 0.0_dp) .or. any(u >= 1.0_dp)) then
      value = -huge(1.0_dp)
      return
    end if
    value = 0.0_dp
    do i = 1, size(u,1)
      term = copula_log_density(u(i,:),model)
      if (term <= -0.5_dp*huge(1.0_dp)) then
        value = -huge(1.0_dp)
        return
      end if
      value = value+term
    end do
  end function copula_log_likelihood

  subroutine rosenblatt_transform(u, model, transformed, ok)
    real(dp), intent(in) :: u(:,:)
    type(copula_model), intent(in) :: model
    real(dp), allocatable, intent(out) :: transformed(:,:)
    logical, intent(out) :: ok
    integer :: i
    ok = size(u,2) == 2 .and. model%dimension == 2 .and. model%valid()
    allocate(transformed(size(u,1),2))
    transformed = 0.0_dp
    if (.not. ok) return
    transformed(:,1) = u(:,1)
    do i = 1, size(u,1)
      transformed(i,2) = conditional_cdf(u(i,2),u(i,1),model)
    end do
  end subroutine rosenblatt_transform

  subroutine inverse_rosenblatt(w, model, u, ok)
    real(dp), intent(in) :: w(:,:)
    type(copula_model), intent(in) :: model
    real(dp), allocatable, intent(out) :: u(:,:)
    logical, intent(out) :: ok
    real(dp) :: lo, hi, mid
    integer :: i, iteration
    ok = size(w,2) == 2 .and. model%dimension == 2 .and. model%valid()
    allocate(u(size(w,1),2))
    u = 0.0_dp
    if (.not. ok) return
    u(:,1) = w(:,1)
    do i = 1, size(w,1)
      lo = epsilon(1.0_dp)
      hi = 1.0_dp-epsilon(1.0_dp)
      do iteration = 1, 70
        mid = 0.5_dp*(lo+hi)
        if (conditional_cdf(mid,w(i,1),model) < w(i,2)) then
          lo = mid
        else
          hi = mid
        end if
      end do
      u(i,2) = 0.5_dp*(lo+hi)
    end do
  end subroutine inverse_rosenblatt

  real(dp) function numerical_standard_error(u, model) result(se)
    real(dp), intent(in) :: u(:,:)
    type(copula_model), intent(in) :: model
    type(copula_model) :: lower_model, upper_model
    real(dp) :: theta, h, center, lower_value, upper_value, curvature
    theta = get_parameter(model)
    h = 1.0e-4_dp*max(1.0_dp,abs(theta))
    lower_model = model
    upper_model = model
    call set_parameter(lower_model,theta-h)
    call set_parameter(upper_model,theta+h)
    center = copula_log_likelihood(u,model)
    lower_value = copula_log_likelihood(u,lower_model)
    upper_value = copula_log_likelihood(u,upper_model)
    curvature = -(upper_value-2.0_dp*center+lower_value)/(h*h)
    if (curvature > 0.0_dp) then
      se = 1.0_dp/sqrt(curvature)
    else
      se = huge(1.0_dp)
    end if
  end function numerical_standard_error

  subroutine parameter_bounds(family, lo, hi, ok)
    integer, intent(in) :: family
    real(dp), intent(out) :: lo, hi
    logical, intent(out) :: ok
    ok = .true.
    select case (family)
    case (family_gaussian, family_student)
      lo = -0.98_dp
      hi = 0.98_dp
    case (family_clayton)
      lo = 1.0e-4_dp
      hi = 30.0_dp
    case (family_gumbel, family_joe)
      lo = 1.0001_dp
      hi = 20.0_dp
    case (family_frank)
      lo = -30.0_dp
      hi = 30.0_dp
    case (family_amh, family_fgm)
      lo = -0.99_dp
      hi = 0.99_dp
    case (family_plackett)
      lo = 0.02_dp
      hi = 50.0_dp
    case (family_galambos, family_husler_reiss)
      lo = 0.02_dp
      hi = 20.0_dp
    case default
      lo = 0.0_dp
      hi = 0.0_dp
      ok = .false.
    end select
  end subroutine parameter_bounds

  subroutine set_parameter(model, parameter)
    type(copula_model), intent(inout) :: model
    real(dp), intent(in) :: parameter
    if (model%family == family_gaussian .or. model%family == family_student) then
      if (.not. allocated(model%correlation)) allocate(model%correlation(2,2))
      model%correlation = reshape([1.0_dp,parameter,parameter,1.0_dp],[2,2])
    else
      model%theta = parameter
    end if
  end subroutine set_parameter

  real(dp) function get_parameter(model) result(parameter)
    type(copula_model), intent(in) :: model
    if (model%family == family_gaussian .or. model%family == family_student) then
      parameter = model%correlation(1,2)
    else
      parameter = model%theta
    end if
  end function get_parameter

  pure function lower(text) result(output)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: output
    integer :: i, code
    output = text
    do i = 1, len(text)
      code = iachar(output(i:i))
      if (code >= iachar('A') .and. code <= iachar('Z')) output(i:i) = achar(code+32)
    end do
  end function lower
end module copula_fitting
