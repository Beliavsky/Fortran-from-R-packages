! SPDX-License-Identifier: GPL-3.0-or-later
module robstattm_compat
  use robstattm_kinds, only : dp
  use robstattm_types, only : robstattm_control, regression_result, logistic_result, &
    covariance_result, projection_result, location_scale_result, pca_result, linear_test_result
  use robstattm_psi, only : rho_value, rho_prime, rho_second, scale_m, &
    inverse_robust_r_squared, tuning_huber, tuning_bisquare, tuning_opt, tuning_mopt
  use robstattm_regression, only : mm_py_fit, sm_py_fit, lmrob_m_fit, lmrobdet_mm, &
    lmrobdet_dcml, lmrobdet_lin_test, lmrobdet_mm_rfpe
  use robstattm_logistic, only : logreg_by, logreg_wby, logreg_wml
  use robstattm_multivariate, only : loc_scale_m, fast_mve, init_pp, cov_classic, &
    cov_rob, cov_rob_mm, cov_rob_rocke
  use robstattm_pca, only : pca_rob_s, prcomp_rob
  implicit none
  private
  public :: rho, rhoprime, rhoprime2, scalem, invtr2
  public :: huber, bisquare, opt, mopt, optv0, moptv0
  public :: mmpy, smpy, lmrobm, lmrobdetmm, lmrobdetdcml
  public :: lmrobdetlintest, roblineartest, lmrobdetmm_rfpe
  public :: bylogreg, wbylogreg, wmllogreg
  public :: mlocdis, fastmve, initpp, covclassic, covrob, covrobmm, covrobrocke
  public :: multirobu, mmultishr, rockemulti, smpca, prcomprob
contains
  elemental function rho(u, family, cc, standardize) result(value)
    real(dp), intent(in) :: u, cc
    character(len=*), intent(in) :: family
    logical, intent(in), optional :: standardize
    real(dp) :: value
    value = rho_value(u, family, cc, standardize)
  end function rho

  elemental function rhoprime(u, family, cc, standardize) result(value)
    real(dp), intent(in) :: u, cc
    character(len=*), intent(in) :: family
    logical, intent(in), optional :: standardize
    real(dp) :: value
    value = rho_prime(u, family, cc, standardize)
  end function rhoprime

  elemental function rhoprime2(u, family, cc, standardize) result(value)
    real(dp), intent(in) :: u, cc
    character(len=*), intent(in) :: family
    logical, intent(in), optional :: standardize
    real(dp) :: value
    value = rho_second(u, family, cc, standardize)
  end function rhoprime2

  function scalem(u, delta, family, tuning_chi, max_iter, tol) result(value)
    real(dp), intent(in) :: u(:)
    real(dp), intent(in), optional :: delta, tuning_chi, tol
    character(len=*), intent(in), optional :: family
    integer, intent(in), optional :: max_iter
    real(dp) :: value
    value = scale_m(u, delta, family, tuning_chi, max_iter, tol)
  end function scalem

  function invtr2(robust_r_squared, family, cc) result(value)
    real(dp), intent(in) :: robust_r_squared, cc
    character(len=*), intent(in) :: family
    real(dp) :: value
    value = inverse_robust_r_squared(robust_r_squared, family, cc)
  end function invtr2

  function huber(efficiency) result(value)
    real(dp), intent(in) :: efficiency
    real(dp) :: value
    value = tuning_huber(efficiency)
  end function huber

  function bisquare(efficiency) result(value)
    real(dp), intent(in) :: efficiency
    real(dp) :: value
    value = tuning_bisquare(efficiency)
  end function bisquare

  function opt(efficiency) result(value)
    real(dp), intent(in) :: efficiency
    real(dp) :: value
    value = tuning_opt(efficiency)
  end function opt

  function mopt(efficiency) result(value)
    real(dp), intent(in) :: efficiency
    real(dp) :: value
    value = tuning_mopt(efficiency)
  end function mopt

  function optv0(efficiency) result(value)
    real(dp), intent(in) :: efficiency
    real(dp) :: value
    value = tuning_opt(efficiency)
  end function optv0

  function moptv0(efficiency) result(value)
    real(dp), intent(in) :: efficiency
    real(dp) :: value
    value = tuning_mopt(efficiency)
  end function moptv0

  subroutine mmpy(x, y, result, control)
    real(dp), intent(in) :: x(:, :), y(:)
    type(regression_result), intent(out) :: result
    type(robstattm_control), intent(in), optional :: control
    call mm_py_fit(x, y, result, control)
  end subroutine mmpy

  subroutine smpy(x, y, result, control)
    real(dp), intent(in) :: x(:, :), y(:)
    type(regression_result), intent(out) :: result
    type(robstattm_control), intent(in), optional :: control
    call sm_py_fit(x, y, result, control)
  end subroutine smpy

  subroutine lmrobm(x, y, result, control)
    real(dp), intent(in) :: x(:, :), y(:)
    type(regression_result), intent(out) :: result
    type(robstattm_control), intent(in), optional :: control
    call lmrob_m_fit(x, y, result, control)
  end subroutine lmrobm

  subroutine lmrobdetmm(x, y, result, control)
    real(dp), intent(in) :: x(:, :), y(:)
    type(regression_result), intent(out) :: result
    type(robstattm_control), intent(in), optional :: control
    call lmrobdet_mm(x, y, result, control)
  end subroutine lmrobdetmm

  subroutine lmrobdetdcml(x, y, result, control)
    real(dp), intent(in) :: x(:, :), y(:)
    type(regression_result), intent(out) :: result
    type(robstattm_control), intent(in), optional :: control
    call lmrobdet_dcml(x, y, result, control)
  end subroutine lmrobdetdcml

  subroutine lmrobdetlintest(full_fit, restricted_fit, result, control)
    type(regression_result), intent(in) :: full_fit, restricted_fit
    type(linear_test_result), intent(out) :: result
    type(robstattm_control), intent(in), optional :: control
    call lmrobdet_lin_test(full_fit, restricted_fit, result, control)
  end subroutine lmrobdetlintest

  subroutine roblineartest(full_fit, restricted_fit, result, control)
    type(regression_result), intent(in) :: full_fit, restricted_fit
    type(linear_test_result), intent(out) :: result
    type(robstattm_control), intent(in), optional :: control
    call lmrobdet_lin_test(full_fit, restricted_fit, result, control)
  end subroutine roblineartest

  function lmrobdetmm_rfpe(fit, control, scale, minimum_rho, penalty) result(value)
    type(regression_result), intent(in) :: fit
    type(robstattm_control), intent(in), optional :: control
    real(dp), intent(in), optional :: scale
    real(dp), intent(out), optional :: minimum_rho, penalty
    real(dp) :: value
    value = lmrobdet_mm_rfpe(fit, control, scale, minimum_rho, penalty)
  end function lmrobdetmm_rfpe

  subroutine bylogreg(x, y, result, intercept, const, max_iter, tol)
    real(dp), intent(in) :: x(:, :), y(:)
    type(logistic_result), intent(out) :: result
    logical, intent(in), optional :: intercept
    real(dp), intent(in), optional :: const, tol
    integer, intent(in), optional :: max_iter
    call logreg_by(x, y, result, intercept, const, max_iter, tol)
  end subroutine bylogreg

  subroutine wbylogreg(x, y, result, intercept, const, max_iter, tol)
    real(dp), intent(in) :: x(:, :), y(:)
    type(logistic_result), intent(out) :: result
    logical, intent(in), optional :: intercept
    real(dp), intent(in), optional :: const, tol
    integer, intent(in), optional :: max_iter
    call logreg_wby(x, y, result, intercept, const, max_iter, tol)
  end subroutine wbylogreg

  subroutine wmllogreg(x, y, result, intercept, max_iter, tol)
    real(dp), intent(in) :: x(:, :), y(:)
    type(logistic_result), intent(out) :: result
    logical, intent(in), optional :: intercept
    integer, intent(in), optional :: max_iter
    real(dp), intent(in), optional :: tol
    call logreg_wml(x, y, result, intercept, max_iter, tol)
  end subroutine wmllogreg

  subroutine mlocdis(x, result, family, efficiency, max_iter, tol)
    real(dp), intent(in) :: x(:)
    type(location_scale_result), intent(out) :: result
    character(len=*), intent(in), optional :: family
    real(dp), intent(in), optional :: efficiency, tol
    integer, intent(in), optional :: max_iter
    call loc_scale_m(x, result, family, efficiency, max_iter, tol)
  end subroutine mlocdis

  subroutine fastmve(x, result, nsamp, seed)
    real(dp), intent(in) :: x(:, :)
    type(covariance_result), intent(out) :: result
    integer, intent(in), optional :: nsamp, seed
    call fast_mve(x, result, nsamp, seed)
  end subroutine fastmve

  subroutine initpp(x, result, random_multiplier, fixed_multiplier, minimum_directions, seed)
    real(dp), intent(in) :: x(:, :)
    type(projection_result), intent(out) :: result
    integer, intent(in), optional :: random_multiplier, fixed_multiplier, minimum_directions, seed
    call init_pp(x, result, random_multiplier, fixed_multiplier, minimum_directions, seed)
  end subroutine initpp

  subroutine covclassic(x, result, unbiased, correlation)
    real(dp), intent(in) :: x(:, :)
    type(covariance_result), intent(out) :: result
    logical, intent(in), optional :: unbiased, correlation
    call cov_classic(x, result, unbiased, correlation)
  end subroutine covclassic

  subroutine covrob(x, result, method, max_iter, tol, correlation)
    real(dp), intent(in) :: x(:, :)
    type(covariance_result), intent(out) :: result
    character(len=*), intent(in), optional :: method
    integer, intent(in), optional :: max_iter
    real(dp), intent(in), optional :: tol
    logical, intent(in), optional :: correlation
    call cov_rob(x, result, method, max_iter, tol, correlation)
  end subroutine covrob

  subroutine multirobu(x, result, method, max_iter, tol, correlation)
    real(dp), intent(in) :: x(:, :)
    type(covariance_result), intent(out) :: result
    character(len=*), intent(in), optional :: method
    integer, intent(in), optional :: max_iter
    real(dp), intent(in), optional :: tol
    logical, intent(in), optional :: correlation
    call cov_rob(x, result, method, max_iter, tol, correlation)
  end subroutine multirobu

  subroutine covrobmm(x, result, max_iter, tol, correlation, seed)
    real(dp), intent(in) :: x(:, :)
    type(covariance_result), intent(out) :: result
    integer, intent(in), optional :: max_iter, seed
    real(dp), intent(in), optional :: tol
    logical, intent(in), optional :: correlation
    call cov_rob_mm(x, result, max_iter, tol, correlation, seed)
  end subroutine covrobmm

  subroutine mmultishr(x, result, max_iter, tol, correlation, seed)
    real(dp), intent(in) :: x(:, :)
    type(covariance_result), intent(out) :: result
    integer, intent(in), optional :: max_iter, seed
    real(dp), intent(in), optional :: tol
    logical, intent(in), optional :: correlation
    call cov_rob_mm(x, result, max_iter, tol, correlation, seed)
  end subroutine mmultishr

  subroutine covrobrocke(x, result, initial_method, max_steps, proportion_minimum, q, &
      max_iter, tol, correlation, seed)
    real(dp), intent(in) :: x(:, :)
    type(covariance_result), intent(out) :: result
    character(len=*), intent(in), optional :: initial_method
    integer, intent(in), optional :: max_steps, proportion_minimum, q, max_iter, seed
    real(dp), intent(in), optional :: tol
    logical, intent(in), optional :: correlation
    call cov_rob_rocke(x, result, initial_method, max_steps, proportion_minimum, q, &
      max_iter, tol, correlation, seed)
  end subroutine covrobrocke

  subroutine rockemulti(x, result, initial_method, max_steps, proportion_minimum, q, &
      max_iter, tol, correlation, seed)
    real(dp), intent(in) :: x(:, :)
    type(covariance_result), intent(out) :: result
    character(len=*), intent(in), optional :: initial_method
    integer, intent(in), optional :: max_steps, proportion_minimum, q, max_iter, seed
    real(dp), intent(in), optional :: tol
    logical, intent(in), optional :: correlation
    call cov_rob_rocke(x, result, initial_method, max_steps, proportion_minimum, q, &
      max_iter, tol, correlation, seed)
  end subroutine rockemulti

  subroutine smpca(x, result, ncomp, desired_proportion, delta_scale, max_iter, tolerance)
    real(dp), intent(in) :: x(:, :)
    type(pca_result), intent(out) :: result
    integer, intent(in), optional :: ncomp, max_iter
    real(dp), intent(in), optional :: desired_proportion, delta_scale, tolerance
    call pca_rob_s(x, result, ncomp, desired_proportion, delta_scale, max_iter, tolerance)
  end subroutine smpca

  subroutine prcomprob(x, result, rank, delta_scale, max_iter, tolerance)
    real(dp), intent(in) :: x(:, :)
    type(pca_result), intent(out) :: result
    integer, intent(in), optional :: rank, max_iter
    real(dp), intent(in), optional :: delta_scale, tolerance
    call prcomp_rob(x, result, rank, delta_scale, max_iter, tolerance)
  end subroutine prcomprob
end module robstattm_compat
