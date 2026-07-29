! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2016 Marius Hofert, Kurt Hornik and Alexander J. McNeil
module qrmtools_garch
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
  use qrmtools_kinds, only : dp, pi
  use qrmtools_types, only : garch_result
  use qrmtools_stats, only : student_t_pdf
  use qrmtools_optimization, only : nelder_mead
  implicit none
  private

  public :: loglik_garch_11, fit_garch_11, tail_index_garch_11

  real(dp), allocatable, save :: context_x(:)
  real(dp), save :: context_sig2 = 1.0_dp
  real(dp), save :: context_delta = 1.0_dp
  logical, save :: context_student = .false.

contains

  real(dp) function loglik_garch_11(parameters, x, sig2, delta, student) result(value)
    real(dp), intent(in) :: parameters(:)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in) :: sig2
    real(dp), intent(in) :: delta
    logical, intent(in), optional :: student
    real(dp) :: mu_cor
    real(dp) :: mu_ema
    real(dp) :: variance
    real(dp) :: df
    real(dp) :: xlag2
    integer :: i
    logical :: use_student

    use_student = .false.
    if (present(student)) use_student = student

    if (size(parameters) < 2 .or. size(x) < 1 .or. sig2 <= 0.0_dp .or. delta <= 0.0_dp) then
      value = -huge(1.0_dp)
      return
    end if

    if (use_student) then
      if (size(parameters) < 3 .or. parameters(3) <= 2.0_dp) then
        value = -huge(1.0_dp)
        return
      end if
      df = parameters(3)
    else
      df = 0.0_dp
    end if

    mu_cor = exp(-delta / exp(parameters(1)))
    mu_ema = exp(-delta / exp(parameters(2)))
    variance = sig2
    value = 0.0_dp

    xlag2 = sum(x*x) / real(size(x), dp)
    variance = sig2 * (1.0_dp - mu_cor) + &
      mu_cor * (1.0_dp - mu_ema) * xlag2 + mu_cor * mu_ema * variance
    call add_loglik_contribution(x(1), variance, use_student, df, value)
    if (.not. ieee_is_finite(value)) return

    do i = 2, size(x)
      xlag2 = x(i-1)**2
      variance = sig2 * (1.0_dp - mu_cor) + &
        mu_cor * (1.0_dp - mu_ema) * xlag2 + mu_cor * mu_ema * variance
      call add_loglik_contribution(x(i), variance, use_student, df, value)
      if (.not. ieee_is_finite(value)) return
    end do
  end function loglik_garch_11


  subroutine add_loglik_contribution(observation, variance, use_student, df, value)
    real(dp), intent(in) :: observation
    real(dp), intent(in) :: variance
    logical, intent(in) :: use_student
    real(dp), intent(in) :: df
    real(dp), intent(inout) :: value
    real(dp) :: scaled_sigma
    real(dp) :: z

    if (variance <= 0.0_dp .or. .not. ieee_is_finite(variance)) then
      value = -huge(1.0_dp)
      return
    end if

    if (use_student) then
      scaled_sigma = sqrt(variance) * sqrt((df - 2.0_dp) / df)
      z = observation / scaled_sigma
      value = value + log(max(student_t_pdf(z, df), tiny(1.0_dp))) - log(scaled_sigma)
    else
      z = observation / sqrt(variance)
      value = value - 0.5_dp*log(2.0_dp*pi) - 0.5_dp*z*z - 0.5_dp*log(variance)
    end if
  end subroutine add_loglik_contribution

  function fit_garch_11(x, initial, sig2, delta, student, max_iterations) result(output)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in), optional :: initial(:)
    real(dp), intent(in), optional :: sig2
    real(dp), intent(in), optional :: delta
    logical, intent(in), optional :: student
    integer, intent(in), optional :: max_iterations
    type(garch_result) :: output
    real(dp), allocatable :: start(:)
    real(dp), allocatable :: optimum(:)
    real(dp), allocatable :: variance(:)
    real(dp) :: variance0
    real(dp) :: dt
    real(dp) :: best
    real(dp) :: mu_cor
    real(dp) :: mu_ema
    real(dp) :: mx2
    real(dp) :: xlag2
    integer :: iterations
    integer :: evaluations
    integer :: maxit
    integer :: i
    logical :: use_student
    logical :: converged

    if (size(x) < 2) then
      output%message = 'At least two observations are required.'
      return
    end if

    use_student = .false.
    if (present(student)) use_student = student

    variance0 = sum(x*x) / real(size(x), dp)
    if (present(sig2)) variance0 = sig2
    dt = 1.0_dp
    if (present(delta)) dt = delta

    if (variance0 <= 0.0_dp .or. dt <= 0.0_dp) then
      output%message = 'Variance and time interval must be positive.'
      return
    end if

    if (present(initial)) then
      allocate(start(size(initial)))
      start = initial
    else if (use_student) then
      allocate(start(3))
      start = [1.0_dp, 1.0_dp, 4.0_dp]
    else
      allocate(start(2))
      start = [1.0_dp, 1.0_dp]
    end if

    if ((.not. use_student .and. size(start) /= 2) .or. &
        (use_student .and. size(start) /= 3)) then
      output%message = 'Invalid initial parameter vector.'
      return
    end if

    context_x = x
    context_sig2 = variance0
    context_delta = dt
    context_student = use_student

    maxit = 2000
    if (present(max_iterations)) maxit = max_iterations
    call nelder_mead(garch_objective, start, optimum, best, iterations, evaluations, &
      converged, maxit, 1.0e-9_dp)

    if (.not. ieee_is_finite(best) .or. best >= 1.0e149_dp) then
      output%message = 'GARCH optimization failed.'
      return
    end if

    mu_cor = exp(-dt / exp(optimum(1)))
    mu_ema = exp(-dt / exp(optimum(2)))
    mx2 = sum(x*x) / real(size(x), dp)

    if (use_student) then
      output%coefficients = [mx2*(1.0_dp-mu_cor), &
        mu_cor*(1.0_dp-mu_ema), mu_cor*mu_ema, optimum(3)]
    else
      output%coefficients = [mx2*(1.0_dp-mu_cor), &
        mu_cor*(1.0_dp-mu_ema), mu_cor*mu_ema]
    end if

    allocate(variance(size(x)))
    allocate(output%sigma(size(x)))
    allocate(output%residuals(size(x)))

    xlag2 = mx2
    variance(1) = variance0*(1.0_dp-mu_cor) + &
      mu_cor*(1.0_dp-mu_ema)*xlag2 + mu_cor*mu_ema*variance0
    do i = 2, size(x)
      xlag2 = x(i-1)**2
      variance(i) = variance0*(1.0_dp-mu_cor) + &
        mu_cor*(1.0_dp-mu_ema)*xlag2 + mu_cor*mu_ema*variance(i-1)
    end do

    if (any(variance <= 0.0_dp)) then
      output%message = 'Fitted conditional variance is not positive.'
      return
    end if

    output%sigma = sqrt(variance)
    output%residuals = x / output%sigma
    output%log_likelihood = -best
    output%iterations = iterations
    output%evaluations = evaluations
    output%converged = converged
    output%ok = .true.
  end function fit_garch_11

  real(dp) function garch_objective(parameters) result(value)
    real(dp), intent(in) :: parameters(:)
    real(dp) :: log_likelihood

    log_likelihood = loglik_garch_11(parameters, context_x, context_sig2, &
      context_delta, context_student)
    if (ieee_is_finite(log_likelihood)) then
      value = -log_likelihood
    else
      value = 1.0e150_dp
    end if
  end function garch_objective

  real(dp) function tail_index_garch_11(innovations, alpha1, beta1, lower, upper) result(value)
    real(dp), intent(in) :: innovations(:)
    real(dp), intent(in) :: alpha1
    real(dp), intent(in) :: beta1
    real(dp), intent(in), optional :: lower
    real(dp), intent(in), optional :: upper
    real(dp) :: lo
    real(dp) :: hi
    real(dp) :: mid
    real(dp) :: flo
    real(dp) :: fmid
    integer :: i

    if (size(innovations) < 1 .or. alpha1 < 0.0_dp .or. beta1 < 0.0_dp .or. &
        alpha1 + beta1 >= 1.0_dp) then
      value = ieee_value(1.0_dp, ieee_quiet_nan)
      return
    end if

    lo = 1.0e-6_dp
    hi = 100.0_dp
    if (present(lower)) lo = lower
    if (present(upper)) hi = upper

    flo = tail_equation(lo, innovations, alpha1, beta1)
    if (flo * tail_equation(hi, innovations, alpha1, beta1) >= 0.0_dp) then
      value = ieee_value(1.0_dp, ieee_quiet_nan)
      return
    end if

    do i = 1, 160
      mid = 0.5_dp * (lo + hi)
      fmid = tail_equation(mid, innovations, alpha1, beta1)
      if (flo*fmid <= 0.0_dp) then
        hi = mid
      else
        lo = mid
        flo = fmid
      end if
    end do
    value = 0.5_dp * (lo + hi)
  end function tail_index_garch_11

  pure real(dp) function tail_equation(a, innovations, alpha1, beta1) result(value)
    real(dp), intent(in) :: a
    real(dp), intent(in) :: innovations(:)
    real(dp), intent(in) :: alpha1
    real(dp), intent(in) :: beta1

    value = sum((alpha1*innovations**2 + beta1)**(0.5_dp*a)) / &
      real(size(innovations), dp) - 1.0_dp
  end function tail_equation

end module qrmtools_garch
