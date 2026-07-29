! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module dowd_copulas
  use dowd_kinds, only: dp, pi
  use dowd_math, only: normal_cdf, normal_pdf, normal_quantile
  implicit none
  private

  public :: product_copula, gaussian_copula, gumbel_copula
  public :: cdf_sum_product_copula, cdf_sum_gaussian_copula, cdf_sum_gumbel_copula
  public :: product_copula_var, gaussian_copula_var, gumbel_copula_var

contains

  pure real(dp) function product_copula(u, v) result(c)
    real(dp), intent(in) :: u, v
    c = min(max(u,0.0_dp),1.0_dp)*min(max(v,0.0_dp),1.0_dp)
  end function product_copula

  pure real(dp) function gumbel_copula(u, v, beta) result(c)
    real(dp), intent(in) :: u, v, beta
    real(dp) :: uc, vc, a
    if (beta < 1.0_dp) error stop "gumbel_copula: beta must be at least 1"
    uc = min(max(u,tiny(1.0_dp)),1.0_dp)
    vc = min(max(v,tiny(1.0_dp)),1.0_dp)
    if (uc >= 1.0_dp .and. vc >= 1.0_dp) then
      c = 1.0_dp
    else
      a = (-log(uc))**beta+(-log(vc))**beta
      c = exp(-a**(1.0_dp/beta))
    end if
  end function gumbel_copula

  real(dp) function gaussian_copula(u, v, rho, number_steps) result(c)
    real(dp), intent(in) :: u, v, rho
    integer, intent(in), optional :: number_steps
    real(dp) :: a, b, h, x, total, conditional_sd
    integer :: i, n
    if (abs(rho) >= 1.0_dp) error stop "gaussian_copula: abs(rho) must be less than 1"
    if (u <= 0.0_dp .or. v <= 0.0_dp) then
      c = 0.0_dp
      return
    else if (u >= 1.0_dp) then
      c = min(v,1.0_dp)
      return
    else if (v >= 1.0_dp) then
      c = min(u,1.0_dp)
      return
    end if
    a = normal_quantile(u)
    b = normal_quantile(v)
    n = 2000
    if (present(number_steps)) n = max(100,number_steps)
    if (mod(n,2) /= 0) n = n+1
    h = (a+9.0_dp)/real(n,dp)
    conditional_sd = sqrt(1.0_dp-rho*rho)
    total = 0.0_dp
    do i = 0, n
      x = -9.0_dp+real(i,dp)*h
      if (i == 0 .or. i == n) then
        total = total+normal_pdf(x)*normal_cdf((b-rho*x)/conditional_sd)
      else if (mod(i,2) == 0) then
        total = total+2.0_dp*normal_pdf(x)*normal_cdf((b-rho*x)/conditional_sd)
      else
        total = total+4.0_dp*normal_pdf(x)*normal_cdf((b-rho*x)/conditional_sd)
      end if
    end do
    c = min(max(total*h/3.0_dp,0.0_dp),1.0_dp)
  end function gaussian_copula

  pure real(dp) function cdf_sum_product_copula(quantile, mu1, mu2, sigma1, sigma2) result(p)
    real(dp), intent(in) :: quantile, mu1, mu2, sigma1, sigma2
    real(dp) :: variance
    variance = sigma1*sigma1+sigma2*sigma2
    if (variance <= 0.0_dp) then
      p = merge(1.0_dp,0.0_dp,quantile >= mu1+mu2)
    else
      p = normal_cdf((quantile-mu1-mu2)/sqrt(variance))
    end if
  end function cdf_sum_product_copula

  pure real(dp) function cdf_sum_gaussian_copula(quantile, mu1, mu2, sigma1, sigma2, rho) result(p)
    real(dp), intent(in) :: quantile, mu1, mu2, sigma1, sigma2, rho
    real(dp) :: variance
    if (abs(rho) > 1.0_dp) error stop "cdf_sum_gaussian_copula: invalid rho"
    variance = sigma1*sigma1+sigma2*sigma2+2.0_dp*rho*sigma1*sigma2
    if (variance <= 0.0_dp) then
      p = merge(1.0_dp,0.0_dp,quantile >= mu1+mu2)
    else
      p = normal_cdf((quantile-mu1-mu2)/sqrt(variance))
    end if
  end function cdf_sum_gaussian_copula

  pure real(dp) function gumbel_partial_u(u, v, beta) result(derivative)
    real(dp), intent(in) :: u, v, beta
    real(dp) :: uc, vc, a, c
    uc = min(max(u,1.0e-14_dp),1.0_dp-1.0e-14_dp)
    vc = min(max(v,1.0e-14_dp),1.0_dp-1.0e-14_dp)
    a = (-log(uc))**beta+(-log(vc))**beta
    c = exp(-a**(1.0_dp/beta))
    derivative = c*a**(1.0_dp/beta-1.0_dp)*(-log(uc))**(beta-1.0_dp)/uc
  end function gumbel_partial_u

  real(dp) function cdf_sum_gumbel_copula(quantile, mu1, mu2, sigma1, sigma2, beta, number_steps) result(p)
    real(dp), intent(in) :: quantile, mu1, mu2, sigma1, sigma2, beta
    integer, intent(in), optional :: number_steps
    real(dp) :: lo, hi, h, u, x, v, integrand, total
    integer :: i, n
    if (beta < 1.0_dp .or. sigma1 <= 0.0_dp .or. sigma2 <= 0.0_dp) &
      error stop "cdf_sum_gumbel_copula: invalid parameters"
    n = 4000
    if (present(number_steps)) n = max(200,number_steps)
    if (mod(n,2) /= 0) n = n+1
    lo = 1.0e-8_dp
    hi = 1.0_dp-1.0e-8_dp
    h = (hi-lo)/real(n,dp)
    total = 0.0_dp
    do i = 0, n
      u = lo+real(i,dp)*h
      x = mu1+sigma1*normal_quantile(u)
      v = normal_cdf((quantile-x-mu2)/sigma2)
      integrand = gumbel_partial_u(u,v,beta)
      if (i == 0 .or. i == n) then
        total = total+integrand
      else if (mod(i,2) == 0) then
        total = total+2.0_dp*integrand
      else
        total = total+4.0_dp*integrand
      end if
    end do
    p = min(max(total*h/3.0_dp,0.0_dp),1.0_dp)
  end function cdf_sum_gumbel_copula

  pure real(dp) function product_copula_var(mu1, mu2, sigma1, sigma2, cl) result(value)
    real(dp), intent(in) :: mu1, mu2, sigma1, sigma2, cl
    value = -(mu1+mu2)-sqrt(sigma1*sigma1+sigma2*sigma2)*normal_quantile(1.0_dp-cl)
  end function product_copula_var

  pure real(dp) function gaussian_copula_var(mu1, mu2, sigma1, sigma2, rho, cl) result(value)
    real(dp), intent(in) :: mu1, mu2, sigma1, sigma2, rho, cl
    value = -(mu1+mu2)-sqrt(max(0.0_dp,sigma1*sigma1+sigma2*sigma2+2.0_dp*rho*sigma1*sigma2))* &
            normal_quantile(1.0_dp-cl)
  end function gaussian_copula_var

  real(dp) function gumbel_copula_var(mu1, mu2, sigma1, sigma2, beta, cl, number_steps) result(value)
    real(dp), intent(in) :: mu1, mu2, sigma1, sigma2, beta, cl
    integer, intent(in), optional :: number_steps
    real(dp) :: p, lo, hi, mid, sd_scale, cdf_mid
    integer :: iter, n
    p = 1.0_dp-cl
    sd_scale = sqrt(sigma1*sigma1+sigma2*sigma2)
    lo = mu1+mu2-12.0_dp*sd_scale
    hi = mu1+mu2+12.0_dp*sd_scale
    n = 2000
    if (present(number_steps)) n = number_steps
    do iter = 1, 80
      mid = 0.5_dp*(lo+hi)
      cdf_mid = cdf_sum_gumbel_copula(mid,mu1,mu2,sigma1,sigma2,beta,n)
      if (cdf_mid < p) then
        lo = mid
      else
        hi = mid
      end if
    end do
    value = -0.5_dp*(lo+hi)
  end function gumbel_copula_var

end module dowd_copulas
