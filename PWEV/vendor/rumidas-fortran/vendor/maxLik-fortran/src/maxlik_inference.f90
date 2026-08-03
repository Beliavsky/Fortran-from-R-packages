! SPDX-License-Identifier: GPL-2.0-or-later
module maxlik_inference
  use maxlik_kinds, only: dp
  use maxlik_types, only: maxlik_problem, maxlik_control, maxlik_result, derivative_comparison
  use maxlik_status, only: MAXLIK_EVALUATION_ERROR, MAXLIK_SINGULAR_HESSIAN
  use maxlik_linalg, only: invert_matrix, symmetric_condition_number
  use maxlik_evaluation, only: evaluate_value, evaluate_gradient, evaluate_hessian, evaluate_scores, &
    numeric_gradient, numeric_hessian
  implicit none
  private

  public :: finalize_result, covariance_matrix, robust_covariance_matrix
  public :: standard_errors, normal_confidence_intervals, maxlik_aic
  public :: compare_derivatives, condition_number

contains

  subroutine finalize_result(problem, control, result)
    type(maxlik_problem), intent(in) :: problem
    type(maxlik_control), intent(in) :: control
    type(maxlik_result), intent(inout) :: result

    integer :: status
    real(dp) :: value

    if (.not. allocated(result%estimate)) return
    call evaluate_value(problem, result%estimate, value, result%function_count, status, include_penalty=.false.)
    if (status == 0) result%maximum = value
    call evaluate_gradient(problem, result%estimate, result%gradient, result%function_count, &
      result%gradient_count, status, control%use_central_differences, include_penalty=.false.)
    if (status /= 0) then
      result%code = MAXLIK_EVALUATION_ERROR
      result%converged = .false.
      return
    end if

    if (associated(problem%scores) .and. problem%nobs > 0) then
      allocate(result%gradient_obs(problem%nobs, problem%npar))
      call evaluate_scores(problem, result%estimate, result%gradient_obs, result%gradient_count, status)
      if (status /= 0) deallocate(result%gradient_obs)
    end if

    if (.not. control%final_hessian) return
    call evaluate_hessian(problem, result%estimate, result%hessian, result%function_count, &
      result%gradient_count, result%hessian_count, status, include_penalty=.false.)
    if (status /= 0) then
      result%code = MAXLIK_EVALUATION_ERROR
      result%converged = .false.
      return
    end if
    call covariance_matrix(result%hessian, result%active, result%covariance, status)
    if (status == 0) then
      allocate(result%std_error(problem%npar))
      call standard_errors(result%covariance, result%std_error)
      result%condition_number = condition_number(result%hessian, result%active, status)
    else
      result%condition_number = huge(1.0_dp)
    end if
  end subroutine finalize_result

  subroutine covariance_matrix(hessian, active, covariance, status)
    real(dp), intent(in) :: hessian(:, :)
    logical, intent(in) :: active(:)
    real(dp), allocatable, intent(out) :: covariance(:, :)
    integer, intent(out) :: status

    real(dp), allocatable :: information(:, :), inverse(:, :)
    integer, allocatable :: index(:)
    integer :: i, j, m, n

    n = size(active)
    status = MAXLIK_SINGULAR_HESSIAN
    allocate(covariance(n, n))
    covariance = 0.0_dp
    if (size(hessian, 1) /= n .or. size(hessian, 2) /= n) return
    m = count(active)
    if (m == 0) then
      status = 0
      return
    end if
    allocate(index(m))
    index = pack([(i, i=1,n)], active)
    allocate(information(m, m), inverse(m, m))
    do j = 1, m
      do i = 1, m
        information(i, j) = -hessian(index(i), index(j))
      end do
    end do
    call invert_matrix(information, inverse, status)
    if (status /= 0) then
      status = MAXLIK_SINGULAR_HESSIAN
      return
    end if
    do j = 1, m
      do i = 1, m
        covariance(index(i), index(j)) = inverse(i, j)
      end do
    end do
    covariance = 0.5_dp * (covariance + transpose(covariance))
    status = 0
  end subroutine covariance_matrix

  subroutine robust_covariance_matrix(hessian, scores, active, covariance, status)
    real(dp), intent(in) :: hessian(:, :), scores(:, :)
    logical, intent(in) :: active(:)
    real(dp), allocatable, intent(out) :: covariance(:, :)
    integer, intent(out) :: status

    real(dp), allocatable :: bread(:, :), meat(:, :)

    call covariance_matrix(hessian, active, bread, status)
    if (status /= 0) return
    if (size(scores, 2) /= size(active)) then
      status = MAXLIK_SINGULAR_HESSIAN
      return
    end if
    allocate(meat(size(active), size(active)))
    meat = matmul(transpose(scores), scores)
    allocate(covariance(size(active), size(active)))
    covariance = matmul(bread, matmul(meat, bread))
    covariance = 0.5_dp * (covariance + transpose(covariance))
    status = 0
  end subroutine robust_covariance_matrix

  pure subroutine standard_errors(covariance, std_error)
    real(dp), intent(in) :: covariance(:, :)
    real(dp), intent(out) :: std_error(:)
    integer :: i

    do i = 1, size(std_error)
      std_error(i) = sqrt(max(0.0_dp, covariance(i, i)))
    end do
  end subroutine standard_errors

  pure subroutine normal_confidence_intervals(estimate, std_error, z_value, lower, upper)
    real(dp), intent(in) :: estimate(:), std_error(:), z_value
    real(dp), intent(out) :: lower(:), upper(:)
    lower = estimate - z_value * std_error
    upper = estimate + z_value * std_error
  end subroutine normal_confidence_intervals

  pure real(dp) function maxlik_aic(maximum, number_parameters) result(value)
    real(dp), intent(in) :: maximum
    integer, intent(in) :: number_parameters
    value = 2.0_dp * real(number_parameters, dp) - 2.0_dp * maximum
  end function maxlik_aic

  real(dp) function condition_number(hessian, active, status) result(value)
    real(dp), intent(in) :: hessian(:, :)
    logical, intent(in) :: active(:)
    integer, intent(out) :: status

    real(dp), allocatable :: information(:, :)
    integer, allocatable :: index(:)
    integer :: i, j, m, n

    n = size(active)
    m = count(active)
    if (m == 0) then
      value = 1.0_dp
      status = 0
      return
    end if
    allocate(index(m), information(m, m))
    index = pack([(i, i=1,n)], active)
    do j = 1, m
      do i = 1, m
        information(i, j) = -hessian(index(i), index(j))
      end do
    end do
    value = symmetric_condition_number(information, status)
  end function condition_number

  subroutine compare_derivatives(problem, x, comparison, tolerance)
    type(maxlik_problem), intent(in) :: problem
    real(dp), intent(in) :: x(:)
    type(derivative_comparison), intent(out) :: comparison
    real(dp), intent(in), optional :: tolerance

    type(maxlik_problem) :: raw_problem
    real(dp) :: tol
    integer :: function_count, gradient_count, status

    tol = 1.0e-5_dp
    if (present(tolerance)) tol = tolerance
    function_count = 0
    gradient_count = 0
    raw_problem = problem
    raw_problem%penalty_rho = 0.0_dp

    allocate(comparison%numeric_gradient(problem%npar), comparison%gradient_error(problem%npar))
    call numeric_gradient(raw_problem, x, comparison%numeric_gradient, function_count, status, central=.true.)
    if (status /= 0) then
      comparison%status = status
      return
    end if
    if (associated(problem%gradient)) then
      allocate(comparison%analytic_gradient(problem%npar))
      call problem%gradient(x, comparison%analytic_gradient, status)
      if (status /= 0) then
        comparison%status = status
        return
      end if
      comparison%gradient_error = comparison%analytic_gradient - comparison%numeric_gradient
      comparison%max_gradient_error = maxval(abs(comparison%gradient_error))
      comparison%gradient_ok = comparison%max_gradient_error <= tol
    else
      comparison%gradient_error = 0.0_dp
      comparison%max_gradient_error = 0.0_dp
      comparison%gradient_ok = .true.
    end if

    allocate(comparison%numeric_hessian(problem%npar, problem%npar), &
      comparison%hessian_error(problem%npar, problem%npar))
    call numeric_hessian(raw_problem, x, comparison%numeric_hessian, function_count, gradient_count, status)
    if (status /= 0) then
      comparison%status = status
      return
    end if
    if (associated(problem%hessian)) then
      allocate(comparison%analytic_hessian(problem%npar, problem%npar))
      call problem%hessian(x, comparison%analytic_hessian, status)
      if (status /= 0) then
        comparison%status = status
        return
      end if
      comparison%hessian_error = comparison%analytic_hessian - comparison%numeric_hessian
      comparison%max_hessian_error = maxval(abs(comparison%hessian_error))
      comparison%hessian_ok = comparison%max_hessian_error <= 10.0_dp * tol
    else
      comparison%hessian_error = 0.0_dp
      comparison%max_hessian_error = 0.0_dp
      comparison%hessian_ok = .true.
    end if
    comparison%status = 0
  end subroutine compare_derivatives

end module maxlik_inference
