module rcppnumerical_logistic
  use rcppnumerical_kinds, only : dp
  use rcppnumerical_optimization, only : optimization_result_t, optim_lbfgs
  implicit none
  private

  type, public :: logistic_fit_t
    real(dp), allocatable :: coefficients(:)
    real(dp), allocatable :: fitted_values(:)
    real(dp), allocatable :: linear_predictors(:)
    real(dp) :: log_likelihood = -huge(1.0_dp)
    integer :: status = -1
    integer :: iterations = 0
    integer :: evaluations = 0
    logical :: converged = .false.
  end type logistic_fit_t

  public :: fast_lr

contains

  subroutine fast_lr(x, y, fit, start, eps_f, eps_g, maxit)
    real(dp), intent(in) :: x(:,:), y(:)
    type(logistic_fit_t), intent(out) :: fit
    real(dp), intent(in), optional :: start(:)
    real(dp), intent(in), optional :: eps_f, eps_g
    integer, intent(in), optional :: maxit

    type(optimization_result_t) :: opt
    real(dp), allocatable :: beta(:)
    integer :: n, p, iterations
    real(dp) :: f_tol, g_tol

    n = size(x, 1)
    p = size(x, 2)
    allocate(fit%coefficients(p), fit%fitted_values(n), fit%linear_predictors(n))
    fit%coefficients = 0.0_dp
    fit%fitted_values = 0.0_dp
    fit%linear_predictors = 0.0_dp
    if (size(y) /= n .or. any(y < 0.0_dp) .or. any(y > 1.0_dp)) return
    allocate(beta(p))
    beta = 0.0_dp
    if (present(start)) then
      if (size(start) /= p) return
      beta = start
    end if
    f_tol = 1.0e-8_dp
    if (present(eps_f)) f_tol = eps_f
    g_tol = 1.0e-5_dp
    if (present(eps_g)) g_tol = eps_g
    iterations = 300
    if (present(maxit)) iterations = maxit

    call optim_lbfgs(logistic_objective, beta, opt, iterations, f_tol, g_tol)
    fit%coefficients = beta
    fit%linear_predictors = matmul(x, beta)
    fit%fitted_values = logistic(fit%linear_predictors)
    fit%log_likelihood = -negative_log_likelihood(fit%linear_predictors, y)
    fit%status = opt%status
    fit%iterations = opt%iterations
    fit%evaluations = opt%evaluations
    fit%converged = opt%converged

  contains

    subroutine logistic_objective(beta_current, value, gradient, user_data)
      real(dp), intent(in) :: beta_current(:)
      real(dp), intent(out) :: value
      real(dp), intent(out) :: gradient(:)
      class(*), intent(inout), optional :: user_data
      real(dp) :: eta(n), probability(n)
      eta = matmul(x, beta_current)
      probability = logistic(eta)
      value = negative_log_likelihood(eta, y)
      gradient = matmul(transpose(x), probability - y)
    end subroutine logistic_objective

  end subroutine fast_lr

  pure elemental real(dp) function logistic(eta)
    real(dp), intent(in) :: eta
    if (eta >= 0.0_dp) then
      logistic = 1.0_dp/(1.0_dp + exp(-eta))
    else
      logistic = exp(eta)/(1.0_dp + exp(eta))
    end if
  end function logistic

  pure real(dp) function negative_log_likelihood(eta, y)
    real(dp), intent(in) :: eta(:), y(:)
    integer :: i
    negative_log_likelihood = 0.0_dp
    do i = 1, size(eta)
      if (eta(i) > 0.0_dp) then
        negative_log_likelihood = negative_log_likelihood + &
          eta(i) + log(1.0_dp + exp(-eta(i))) - y(i)*eta(i)
      else
        negative_log_likelihood = negative_log_likelihood + &
          log(1.0_dp + exp(eta(i))) - y(i)*eta(i)
      end if
    end do
  end function negative_log_likelihood

end module rcppnumerical_logistic
