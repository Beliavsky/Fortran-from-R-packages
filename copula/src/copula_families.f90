! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2012 Marius Hofert, Ivan Kojadinovic, Martin Maechler and Jun Yan
module copula_families
  use copula_kinds, only : dp, pi
  use copula_types
  use copula_special, only : normal_pdf, normal_cdf, normal_quantile, student_pdf, student_cdf, &
    student_quantile, bivariate_normal_cdf, bivariate_student_cdf, gamma_quantile, &
    log_one_plus, exp_minus_one
  use copula_linalg, only : cholesky_lower, quadratic_form
  use copula_random, only : halton
  implicit none
  private
  public :: copula_cdf, copula_density, copula_log_density, conditional_cdf
  public :: pickands_function, pickands_derivative
contains
  real(dp) function copula_cdf(u, model) result(value)
    real(dp), intent(in) :: u(:)
    type(copula_model), intent(in) :: model
    real(dp) :: transformed(2), base
    if (size(u) /= model%dimension .or. .not. model%valid()) then
      value = 0.0_dp
      return
    end if
    if (any(u <= 0.0_dp)) then
      value = 0.0_dp
      return
    end if
    if (all(u >= 1.0_dp)) then
      value = 1.0_dp
      return
    end if
    if (model%rotation == rotation_none) then
      value = base_cdf(min(1.0_dp,max(0.0_dp,u)),model)
      return
    end if
    select case (model%rotation)
    case (rotation_90)
      transformed = [1.0_dp-u(1),u(2)]
      base = base_cdf(transformed,model)
      value = u(2)-base
    case (rotation_180)
      transformed = 1.0_dp-u
      base = base_cdf(transformed,model)
      value = u(1)+u(2)-1.0_dp+base
    case (rotation_270)
      transformed = [u(1),1.0_dp-u(2)]
      base = base_cdf(transformed,model)
      value = u(1)-base
    case default
      value = 0.0_dp
    end select
    value = min(minval(u),max(0.0_dp,value))
  end function copula_cdf

  real(dp) function base_cdf(u, model) result(value)
    real(dp), intent(in) :: u(:)
    type(copula_model), intent(in) :: model
    real(dp) :: s, a, b, disc, uv, t
    integer :: d
    d = size(u)
    select case (model%family)
    case (family_independence)
      value = product(u)
    case (family_gaussian)
      value = gaussian_cdf(u,model%correlation)
    case (family_student)
      value = student_copula_cdf(u,model%correlation,model%df)
    case (family_clayton)
      s = sum(u**(-model%theta))-real(d,dp)+1.0_dp
      value = max(0.0_dp,s)**(-1.0_dp/model%theta)
    case (family_gumbel)
      s = sum((-log(u))**model%theta)
      value = exp(-s**(1.0_dp/model%theta))
    case (family_frank)
      if (abs(model%theta) < 1.0e-8_dp) then
        value = product(u)
      else
        s = product(exp_minus_one(-model%theta*u))/exp_minus_one(-model%theta)**(d-1)
        value = -log_one_plus(s)/model%theta
      end if
    case (family_amh)
      value = u(1)*u(2)/(1.0_dp-model%theta*(1.0_dp-u(1))*(1.0_dp-u(2)))
    case (family_joe)
      a = (1.0_dp-u(1))**model%theta
      b = (1.0_dp-u(2))**model%theta
      value = 1.0_dp-(a+b-a*b)**(1.0_dp/model%theta)
    case (family_fgm)
      value = u(1)*u(2)*(1.0_dp+model%theta*(1.0_dp-u(1))*(1.0_dp-u(2)))
    case (family_plackett)
      if (abs(model%theta-1.0_dp) < 1.0e-8_dp) then
        value = u(1)*u(2)
      else
        a = 1.0_dp+(model%theta-1.0_dp)*(u(1)+u(2))
        disc = max(0.0_dp,a*a-4.0_dp*model%theta*(model%theta-1.0_dp)*u(1)*u(2))
        value = (a-sqrt(disc))/(2.0_dp*(model%theta-1.0_dp))
      end if
    case (family_marshall_olkin)
      value = min(u(1)**(1.0_dp-model%alpha1)*u(2), &
        u(1)*u(2)**(1.0_dp-model%alpha2))
    case (family_lower_fh)
      value = max(u(1)+u(2)-1.0_dp,0.0_dp)
    case (family_upper_fh)
      value = min(u(1),u(2))
    case (family_galambos, family_husler_reiss, family_tawn)
      uv = u(1)*u(2)
      if (uv <= 0.0_dp) then
        value = 0.0_dp
      else if (uv >= 1.0_dp) then
        value = 1.0_dp
      else
        t = log(u(2))/log(uv)
        value = exp(log(uv)*pickands_function(t,model))
      end if
    case default
      value = 0.0_dp
    end select
    value = min(minval(u),max(0.0_dp,value))
  end function base_cdf

  real(dp) function copula_density(u, model) result(value)
    real(dp), intent(in) :: u(:)
    type(copula_model), intent(in) :: model
    type(copula_model) :: base_model
    real(dp) :: transformed(2)
    if (size(u) /= model%dimension .or. .not. model%valid() .or. any(u <= 0.0_dp) .or. &
        any(u >= 1.0_dp)) then
      value = 0.0_dp
      return
    end if
    base_model = model
    base_model%rotation = rotation_none
    if (model%rotation == rotation_none) then
      value = base_density(u,base_model)
    else
      select case (model%rotation)
      case (rotation_90)
        transformed = [1.0_dp-u(1),u(2)]
      case (rotation_180)
        transformed = 1.0_dp-u
      case (rotation_270)
        transformed = [u(1),1.0_dp-u(2)]
      case default
        transformed = u
      end select
      value = base_density(transformed,base_model)
    end if
    value = max(0.0_dp,value)
  end function copula_density

  real(dp) function copula_log_density(u, model) result(value)
    real(dp), intent(in) :: u(:)
    type(copula_model), intent(in) :: model
    real(dp) :: density
    density = copula_density(u,model)
    if (density <= 0.0_dp) then
      value = -huge(1.0_dp)
    else
      value = log(density)
    end if
  end function copula_log_density

  real(dp) function base_density(u, model) result(value)
    real(dp), intent(in) :: u(:)
    type(copula_model), intent(in) :: model
    real(dp), allocatable :: z(:)
    real(dp) :: q, logdet, log_joint, log_marginal, h1, h2, l1, r1, l2, r2
    logical :: ok
    integer :: i, d
    d = size(u)
    select case (model%family)
    case (family_independence)
      value = 1.0_dp
    case (family_gaussian)
      allocate(z(d))
      do i = 1, d
        z(i) = normal_quantile(u(i))
      end do
      q = quadratic_form(model%correlation,z,logdet,ok)
      if (.not. ok) then
        value = 0.0_dp
      else
        value = exp(-0.5_dp*logdet-0.5_dp*(q-dot_product(z,z)))
      end if
    case (family_student)
      allocate(z(d))
      do i = 1, d
        z(i) = student_quantile(u(i),model%df)
      end do
      q = quadratic_form(model%correlation,z,logdet,ok)
      if (.not. ok) then
        value = 0.0_dp
      else
        log_joint = log_gamma(0.5_dp*(model%df+real(d,dp)))-log_gamma(0.5_dp*model%df) - &
          0.5_dp*logdet-0.5_dp*real(d,dp)*log(model%df*pi) - &
          0.5_dp*(model%df+real(d,dp))*log_one_plus(q/model%df)
        log_marginal = 0.0_dp
        do i = 1, d
          log_marginal = log_marginal+log(student_pdf(z(i),model%df))
        end do
        value = exp(log_joint-log_marginal)
      end if
    case (family_marshall_olkin, family_lower_fh, family_upper_fh)
      value = 0.0_dp
    case default
      if (d /= 2) then
        value = 0.0_dp
        return
      end if
      h1 = min(2.0e-5_dp,0.25_dp*min(u(1),1.0_dp-u(1)))
      h2 = min(2.0e-5_dp,0.25_dp*min(u(2),1.0_dp-u(2)))
      l1 = u(1)-h1
      r1 = u(1)+h1
      l2 = u(2)-h2
      r2 = u(2)+h2
      value = (base_cdf([r1,r2],model)-base_cdf([r1,l2],model)- &
        base_cdf([l1,r2],model)+base_cdf([l1,l2],model))/(4.0_dp*h1*h2)
    end select
  end function base_density

  real(dp) function conditional_cdf(v, given_u, model) result(value)
    real(dp), intent(in) :: v, given_u
    type(copula_model), intent(in) :: model
    type(copula_model) :: base_model
    real(dp) :: zu, zv, rho, h, left, right
    if (model%dimension /= 2 .or. .not. model%valid()) then
      value = 0.0_dp
      return
    end if
    if (v <= 0.0_dp) then
      value = 0.0_dp
      return
    else if (v >= 1.0_dp) then
      value = 1.0_dp
      return
    end if
    base_model = model
    base_model%rotation = rotation_none
    if (model%rotation /= rotation_none) then
      h = min(1.0e-5_dp,0.25_dp*min(given_u,1.0_dp-given_u))
      left = max(0.0_dp,given_u-h)
      right = min(1.0_dp,given_u+h)
      value = (copula_cdf([right,v],model)-copula_cdf([left,v],model))/(right-left)
      value = min(1.0_dp,max(0.0_dp,value))
      return
    end if
    select case (model%family)
    case (family_gaussian)
      rho = model%correlation(1,2)
      zu = normal_quantile(given_u)
      zv = normal_quantile(v)
      value = normal_cdf((zv-rho*zu)/sqrt(1.0_dp-rho*rho))
    case (family_student)
      rho = model%correlation(1,2)
      zu = student_quantile(given_u,model%df)
      zv = student_quantile(v,model%df)
      value = student_cdf((zv-rho*zu)*sqrt((model%df+1.0_dp)/ &
        ((model%df+zu*zu)*(1.0_dp-rho*rho))),model%df+1.0_dp)
    case default
      h = min(1.0e-5_dp,0.25_dp*min(given_u,1.0_dp-given_u))
      left = max(0.0_dp,given_u-h)
      right = min(1.0_dp,given_u+h)
      value = (base_cdf([right,v],base_model)-base_cdf([left,v],base_model))/(right-left)
    end select
    value = min(1.0_dp,max(0.0_dp,value))
  end function conditional_cdf

  real(dp) function gaussian_cdf(u, correlation) result(value)
    real(dp), intent(in) :: u(:), correlation(:,:)
    real(dp), allocatable :: z(:), l(:,:), y(:)
    real(dp) :: estimate, product_probability, conditional_probability, shift
    logical :: ok
    integer :: d, i, j, sample, n_samples
    d = size(u)
    if (d == 2) then
      value = bivariate_normal_cdf(normal_quantile(u(1)),normal_quantile(u(2)),correlation(1,2))
      return
    end if
    call cholesky_lower(correlation,l,ok)
    if (.not. ok) then
      value = 0.0_dp
      return
    end if
    allocate(z(d),y(d))
    do i = 1, d
      z(i) = normal_quantile(u(i))
    end do
    n_samples = 4096
    estimate = 0.0_dp
    do sample = 1, n_samples
      product_probability = 1.0_dp
      y = 0.0_dp
      do i = 1, d
        shift = 0.0_dp
        do j = 1, i-1
          shift = shift+l(i,j)*y(j)
        end do
        conditional_probability = normal_cdf((z(i)-shift)/l(i,i))
        conditional_probability = min(1.0_dp,max(0.0_dp,conditional_probability))
        product_probability = product_probability*conditional_probability
        if (i < d) y(i) = normal_quantile(halton(sample,prime_for(i))*conditional_probability)
      end do
      estimate = estimate+product_probability
    end do
    value = estimate/real(n_samples,dp)
  end function gaussian_cdf

  real(dp) function student_copula_cdf(u, correlation, df) result(value)
    real(dp), intent(in) :: u(:), correlation(:,:), df
    real(dp), allocatable :: z(:), scaled_u(:)
    real(dp) :: mixing, estimate
    integer :: d, sample, n_samples
    d = size(u)
    if (d == 2) then
      value = bivariate_student_cdf(student_quantile(u(1),df),student_quantile(u(2),df), &
        correlation(1,2),df)
      return
    end if
    allocate(z(d),scaled_u(d))
    do sample = 1, d
      z(sample) = student_quantile(u(sample),df)
    end do
    n_samples = 2048
    estimate = 0.0_dp
    do sample = 1, n_samples
      mixing = gamma_quantile(halton(sample,2),0.5_dp*df,2.0_dp)
      scaled_u = normal_cdf(z*sqrt(mixing/df))
      estimate = estimate+gaussian_cdf(scaled_u,correlation)
    end do
    value = estimate/real(n_samples,dp)
  end function student_copula_cdf

  integer pure function prime_for(index) result(prime)
    integer, intent(in) :: index
    integer, parameter :: primes(20) = [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71]
    prime = primes(min(index,size(primes)))
  end function prime_for

  real(dp) function pickands_function(t, model) result(value)
    real(dp), intent(in) :: t
    type(copula_model), intent(in) :: model
    real(dp) :: tt, left, right, lambda
    tt = min(1.0_dp,max(0.0_dp,t))
    if (tt <= 0.0_dp .or. tt >= 1.0_dp) then
      value = 1.0_dp
      return
    end if
    select case (model%family)
    case (family_gumbel)
      value = ((1.0_dp-tt)**model%theta+tt**model%theta)**(1.0_dp/model%theta)
    case (family_galambos)
      value = 1.0_dp-((1.0_dp-tt)**(-model%theta)+tt**(-model%theta))**(-1.0_dp/model%theta)
    case (family_husler_reiss)
      lambda = model%theta
      left = normal_cdf(lambda+0.5_dp*log((1.0_dp-tt)/tt)/lambda)
      right = normal_cdf(lambda+0.5_dp*log(tt/(1.0_dp-tt))/lambda)
      value = (1.0_dp-tt)*left+tt*right
    case (family_tawn)
      left = model%alpha1*(1.0_dp-tt)
      right = model%alpha2*tt
      value = (1.0_dp-model%alpha1)*(1.0_dp-tt)+(1.0_dp-model%alpha2)*tt + &
        (left**model%theta+right**model%theta)**(1.0_dp/model%theta)
    case default
      value = 1.0_dp
    end select
    value = max(max(tt,1.0_dp-tt),min(1.0_dp,value))
  end function pickands_function

  real(dp) function pickands_derivative(t, model) result(value)
    real(dp), intent(in) :: t
    type(copula_model), intent(in) :: model
    real(dp) :: h, left, right
    h = min(1.0e-5_dp,0.25_dp*min(t,1.0_dp-t))
    if (h <= 0.0_dp) then
      value = 0.0_dp
      return
    end if
    left = pickands_function(t-h,model)
    right = pickands_function(t+h,model)
    value = (right-left)/(2.0_dp*h)
  end function pickands_derivative
end module copula_families
