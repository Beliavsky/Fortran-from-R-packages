! SPDX-License-Identifier: GPL-2.0-or-later
module quantreg_nonlinear
  use quantreg_kinds, only : dp
  use quantreg_types, only : rq_result, nlrq_result
  use quantreg_dense, only : rq_fit_fnb, check_loss_sum
  implicit none
  private
  public :: nlrq_fit, nlrq_model_fn

  abstract interface
    subroutine nlrq_model_fn(theta, fitted, jacobian)
      import dp
      real(dp), intent(in) :: theta(:)
      real(dp), intent(out) :: fitted(:)
      real(dp), intent(out) :: jacobian(:,:)
    end subroutine nlrq_model_fn
  end interface
contains

  subroutine nlrq_fit(y, theta0, tau, model, result, maxiter, tol)
    real(dp), intent(in) :: y(:), theta0(:), tau
    procedure(nlrq_model_fn) :: model
    type(nlrq_result), intent(out) :: result
    integer, intent(in), optional :: maxiter
    real(dp), intent(in), optional :: tol
    real(dp), allocatable :: theta(:), trial(:), fitted(:), fitted_trial(:)
    real(dp), allocatable :: jac(:,:), jac_trial(:,:), pseudo_y(:)
    type(rq_result) :: stepfit
    real(dp) :: ftol, obj, trial_obj, alpha
    integer :: it, maxit, n, p

    n = size(y)
    p = size(theta0)
    maxit = 100
    if (present(maxiter)) maxit = maxiter
    ftol = 1.0e-7_dp
    if (present(tol)) ftol = tol
    allocate(theta(p), trial(p), fitted(n), fitted_trial(n), jac(n,p), jac_trial(n,p), pseudo_y(n))
    theta = theta0
    call model(theta, fitted, jac)
    obj = check_loss_sum(y-fitted, tau)
    result%info = 1

    do it = 1, maxit
      pseudo_y = y - fitted + matmul(jac,theta)
      call rq_fit_fnb(jac, pseudo_y, tau, stepfit)
      if (stepfit%info /= 0) then
        result%info = stepfit%info
        exit
      end if
      alpha = 1.0_dp
      do
        trial = theta + alpha*(stepfit%coefficients-theta)
        call model(trial, fitted_trial, jac_trial)
        trial_obj = check_loss_sum(y-fitted_trial, tau)
        if (trial_obj <= obj .or. alpha <= 1.0e-6_dp) exit
        alpha = 0.5_dp*alpha
      end do
      if (maxval(abs(trial-theta)) <= ftol*(1.0_dp+maxval(abs(theta)))) then
        theta = trial
        fitted = fitted_trial
        obj = trial_obj
        result%info = 0
        exit
      end if
      theta = trial
      fitted = fitted_trial
      jac = jac_trial
      obj = trial_obj
    end do
    allocate(result%coefficients(p), result%residuals(n))
    result%coefficients = theta
    result%residuals = y-fitted
    result%objective = obj
    result%iterations = min(it,maxit)
  end subroutine nlrq_fit

end module quantreg_nonlinear
