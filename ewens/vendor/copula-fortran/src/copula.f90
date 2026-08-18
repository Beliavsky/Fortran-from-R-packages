! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2012 Marius Hofert, Ivan Kojadinovic, Martin Maechler and Jun Yan
module copula
  use copula_kinds, only : dp, i8, pi
  use copula_types
  use copula_special, only : normal_pdf, normal_cdf, normal_quantile, student_pdf, student_cdf, &
    student_quantile, log1mexp, log1pexp, debye1
  use copula_families, only : copula_cdf, copula_density, copula_log_density, conditional_cdf, &
    pickands_function, pickands_derivative
  use copula_simulation, only : random_copula
  use copula_dependence, only : kendall_tau_model, spearman_rho_model, tail_dependence, &
    parameter_from_tau, parameter_from_rho
  use copula_empirical, only : pseudo_observations, empirical_copula, empirical_copula_values, &
    sample_kendall_tau, sample_spearman_rho, independence_test, exchangeability_test, radial_symmetry_test
  use copula_fitting, only : fit_copula, copula_log_likelihood, rosenblatt_transform, inverse_rosenblatt
  use copula_compositions, only : mixture_cdf, mixture_density, random_mixture, khoudraji_cdf, &
    khoudraji_density, nested_clayton_cdf
  use copula_special_discrete, only : stirling_first, stirling_second, eulerian_number, &
    sibuya_pmf, random_sibuya, logseries_pmf, random_logseries
  implicit none
  private
  public :: dp, i8, pi
  public :: copula_model, probability_control, probability_result, fit_result, test_result
  public :: family_independence, family_gaussian, family_student, family_clayton, family_gumbel
  public :: family_frank, family_amh, family_joe, family_fgm, family_plackett
  public :: family_marshall_olkin, family_lower_fh, family_upper_fh
  public :: family_galambos, family_husler_reiss, family_tawn
  public :: rotation_none, rotation_90, rotation_180, rotation_270
  public :: independence_copula, normal_copula, t_copula, clayton_copula, gumbel_copula
  public :: frank_copula, amh_copula, joe_copula, fgm_copula, plackett_copula
  public :: marshall_olkin_copula, lower_fh_copula, upper_fh_copula
  public :: galambos_copula, husler_reiss_copula, tawn_copula, rotated_copula
  public :: pCopula, dCopula, logdCopula, rCopula, cCopula
  public :: tau, rho, lambda, iTau, iRho
  public :: pobs, empirical_copula, empirical_copula_values
  public :: sample_kendall_tau, sample_spearman_rho
  public :: fit_copula, copula_log_likelihood, rosenblatt_transform, inverse_rosenblatt
  public :: independence_test, exchangeability_test, radial_symmetry_test
  public :: mixture_cdf, mixture_density, random_mixture, khoudraji_cdf, khoudraji_density
  public :: nested_clayton_cdf, pickands_function, pickands_derivative
  public :: stirling_first, stirling_second, eulerian_number
  public :: sibuya_pmf, random_sibuya, logseries_pmf, random_logseries
  public :: normal_pdf, normal_cdf, normal_quantile, student_pdf, student_cdf, student_quantile
  public :: log1mexp, log1pexp, debye1
contains
  function independence_copula(dimension) result(model)
    integer, intent(in), optional :: dimension
    type(copula_model) :: model
    model%family = family_independence
    model%dimension = 2
    if (present(dimension)) model%dimension = dimension
  end function independence_copula

  function normal_copula(correlation) result(model)
    real(dp), intent(in) :: correlation(:,:)
    type(copula_model) :: model
    model%family = family_gaussian
    model%dimension = size(correlation,1)
    allocate(model%correlation(size(correlation,1),size(correlation,2)))
    model%correlation = correlation
  end function normal_copula

  function t_copula(correlation, df) result(model)
    real(dp), intent(in) :: correlation(:,:), df
    type(copula_model) :: model
    model%family = family_student
    model%dimension = size(correlation,1)
    model%df = df
    allocate(model%correlation(size(correlation,1),size(correlation,2)))
    model%correlation = correlation
  end function t_copula

  function clayton_copula(theta, dimension) result(model)
    real(dp), intent(in) :: theta
    integer, intent(in), optional :: dimension
    type(copula_model) :: model
    model%family = family_clayton
    model%theta = theta
    model%dimension = 2
    if (present(dimension)) model%dimension = dimension
  end function clayton_copula

  function gumbel_copula(theta, dimension) result(model)
    real(dp), intent(in) :: theta
    integer, intent(in), optional :: dimension
    type(copula_model) :: model
    model%family = family_gumbel
    model%theta = theta
    model%dimension = 2
    if (present(dimension)) model%dimension = dimension
  end function gumbel_copula

  function frank_copula(theta, dimension) result(model)
    real(dp), intent(in) :: theta
    integer, intent(in), optional :: dimension
    type(copula_model) :: model
    model%family = family_frank
    model%theta = theta
    model%dimension = 2
    if (present(dimension)) model%dimension = dimension
  end function frank_copula

  function amh_copula(theta) result(model)
    real(dp), intent(in) :: theta
    type(copula_model) :: model
    model%family = family_amh
    model%theta = theta
    model%dimension = 2
  end function amh_copula

  function joe_copula(theta) result(model)
    real(dp), intent(in) :: theta
    type(copula_model) :: model
    model%family = family_joe
    model%theta = theta
    model%dimension = 2
  end function joe_copula

  function fgm_copula(theta) result(model)
    real(dp), intent(in) :: theta
    type(copula_model) :: model
    model%family = family_fgm
    model%theta = theta
    model%dimension = 2
  end function fgm_copula

  function plackett_copula(theta) result(model)
    real(dp), intent(in) :: theta
    type(copula_model) :: model
    model%family = family_plackett
    model%theta = theta
    model%dimension = 2
  end function plackett_copula

  function marshall_olkin_copula(alpha1, alpha2) result(model)
    real(dp), intent(in) :: alpha1, alpha2
    type(copula_model) :: model
    model%family = family_marshall_olkin
    model%alpha1 = alpha1
    model%alpha2 = alpha2
    model%dimension = 2
  end function marshall_olkin_copula

  function lower_fh_copula() result(model)
    type(copula_model) :: model
    model%family = family_lower_fh
    model%dimension = 2
  end function lower_fh_copula

  function upper_fh_copula() result(model)
    type(copula_model) :: model
    model%family = family_upper_fh
    model%dimension = 2
  end function upper_fh_copula

  function galambos_copula(theta) result(model)
    real(dp), intent(in) :: theta
    type(copula_model) :: model
    model%family = family_galambos
    model%theta = theta
    model%dimension = 2
  end function galambos_copula

  function husler_reiss_copula(theta) result(model)
    real(dp), intent(in) :: theta
    type(copula_model) :: model
    model%family = family_husler_reiss
    model%theta = theta
    model%dimension = 2
  end function husler_reiss_copula

  function tawn_copula(theta, alpha1, alpha2) result(model)
    real(dp), intent(in) :: theta, alpha1, alpha2
    type(copula_model) :: model
    model%family = family_tawn
    model%theta = theta
    model%alpha1 = alpha1
    model%alpha2 = alpha2
    model%dimension = 2
  end function tawn_copula

  function rotated_copula(base, rotation) result(model)
    type(copula_model), intent(in) :: base
    integer, intent(in) :: rotation
    type(copula_model) :: model
    model = base
    model%rotation = rotation
  end function rotated_copula

  real(dp) function pCopula(u, model) result(value)
    real(dp), intent(in) :: u(:)
    type(copula_model), intent(in) :: model
    value = copula_cdf(u,model)
  end function pCopula

  real(dp) function dCopula(u, model) result(value)
    real(dp), intent(in) :: u(:)
    type(copula_model), intent(in) :: model
    value = copula_density(u,model)
  end function dCopula

  real(dp) function logdCopula(u, model) result(value)
    real(dp), intent(in) :: u(:)
    type(copula_model), intent(in) :: model
    value = copula_log_density(u,model)
  end function logdCopula

  subroutine rCopula(n, model, samples, ok, seed)
    integer, intent(in) :: n
    type(copula_model), intent(in) :: model
    real(dp), allocatable, intent(out) :: samples(:,:)
    logical, intent(out) :: ok
    integer(i8), intent(in), optional :: seed
    call random_copula(n,model,samples,ok,seed)
  end subroutine rCopula

  real(dp) function cCopula(v, given_u, model) result(value)
    real(dp), intent(in) :: v, given_u
    type(copula_model), intent(in) :: model
    value = conditional_cdf(v,given_u,model)
  end function cCopula

  real(dp) function tau(model) result(value)
    type(copula_model), intent(in) :: model
    value = kendall_tau_model(model)
  end function tau

  real(dp) function rho(model) result(value)
    type(copula_model), intent(in) :: model
    value = spearman_rho_model(model)
  end function rho

  function lambda(model) result(value)
    type(copula_model), intent(in) :: model
    real(dp) :: value(2)
    value = tail_dependence(model)
  end function lambda

  real(dp) function iTau(family, tau_value) result(theta)
    integer, intent(in) :: family
    real(dp), intent(in) :: tau_value
    theta = parameter_from_tau(family,tau_value)
  end function iTau

  real(dp) function iRho(family, rho_value) result(theta)
    integer, intent(in) :: family
    real(dp), intent(in) :: rho_value
    theta = parameter_from_rho(family,rho_value)
  end function iRho

  subroutine pobs(x, u, ok)
    real(dp), intent(in) :: x(:,:)
    real(dp), allocatable, intent(out) :: u(:,:)
    logical, intent(out) :: ok
    call pseudo_observations(x,u,ok)
  end subroutine pobs
end module copula
