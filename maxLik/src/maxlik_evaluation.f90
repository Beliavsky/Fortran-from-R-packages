! SPDX-License-Identifier: GPL-2.0-or-later
module maxlik_evaluation
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use maxlik_kinds, only: dp
  use maxlik_types, only: maxlik_problem, maxlik_objective
  use maxlik_status, only: MAXLIK_INVALID_INPUT, MAXLIK_EVALUATION_ERROR
  implicit none
  private

  public :: initialize_problem, valid_problem, project_parameters
  public :: set_fixed, clear_fixed, set_bounds
  public :: set_equality_constraints, set_inequality_constraints, clear_constraints
  public :: evaluate_value, evaluate_gradient, evaluate_hessian, evaluate_scores
  public :: numeric_gradient, numeric_hessian, constraint_violation

contains

  subroutine initialize_problem(problem, npar, objective, nobs)
    type(maxlik_problem), intent(out) :: problem
    integer, intent(in) :: npar
    procedure(maxlik_objective) :: objective
    integer, intent(in), optional :: nobs

    problem%npar = npar
    problem%nobs = 0
    if (present(nobs)) problem%nobs = nobs
    problem%objective => objective
    allocate(problem%active(npar), problem%lower(npar), problem%upper(npar))
    problem%active = .true.
    problem%lower = -huge(1.0_dp)
    problem%upper = huge(1.0_dp)
  end subroutine initialize_problem

  logical function valid_problem(problem, x) result(ok)
    type(maxlik_problem), intent(in) :: problem
    real(dp), intent(in) :: x(:)

    ok = problem%npar > 0 .and. size(x) == problem%npar
    ok = ok .and. associated(problem%objective)
    ok = ok .and. allocated(problem%active) .and. allocated(problem%lower) .and. allocated(problem%upper)
    if (.not. ok) return
    ok = size(problem%active) == problem%npar
    ok = ok .and. size(problem%lower) == problem%npar .and. size(problem%upper) == problem%npar
    ok = ok .and. all(problem%lower <= problem%upper)
    ok = ok .and. all(ieee_is_finite(x))
    if (.not. ok) return
    ok = all(x >= problem%lower) .and. all(x <= problem%upper)
    if (allocated(problem%eq_a)) then
      ok = ok .and. allocated(problem%eq_b)
      ok = ok .and. size(problem%eq_a, 2) == problem%npar
      ok = ok .and. size(problem%eq_a, 1) == size(problem%eq_b)
    end if
    if (allocated(problem%ineq_a)) then
      ok = ok .and. allocated(problem%ineq_b)
      ok = ok .and. size(problem%ineq_a, 2) == problem%npar
      ok = ok .and. size(problem%ineq_a, 1) == size(problem%ineq_b)
    end if
  end function valid_problem

  pure subroutine project_parameters(problem, x)
    type(maxlik_problem), intent(in) :: problem
    real(dp), intent(inout) :: x(:)
    x = max(problem%lower, min(problem%upper, x))
  end subroutine project_parameters

  subroutine set_fixed(problem, indices)
    type(maxlik_problem), intent(inout) :: problem
    integer, intent(in) :: indices(:)
    integer :: i
    if (.not. allocated(problem%active)) return
    do i = 1, size(indices)
      if (indices(i) >= 1 .and. indices(i) <= problem%npar) problem%active(indices(i)) = .false.
    end do
  end subroutine set_fixed

  subroutine clear_fixed(problem)
    type(maxlik_problem), intent(inout) :: problem
    if (allocated(problem%active)) problem%active = .true.
  end subroutine clear_fixed

  subroutine set_bounds(problem, lower, upper, status)
    type(maxlik_problem), intent(inout) :: problem
    real(dp), intent(in) :: lower(:), upper(:)
    integer, intent(out) :: status
    status = MAXLIK_INVALID_INPUT
    if (size(lower) /= problem%npar .or. size(upper) /= problem%npar) return
    if (any(lower > upper)) return
    problem%lower = lower
    problem%upper = upper
    status = 0
  end subroutine set_bounds

  subroutine set_equality_constraints(problem, a, b, status)
    type(maxlik_problem), intent(inout) :: problem
    real(dp), intent(in) :: a(:, :), b(:)
    integer, intent(out) :: status
    status = MAXLIK_INVALID_INPUT
    if (size(a, 2) /= problem%npar .or. size(a, 1) /= size(b)) return
    if (allocated(problem%eq_a)) deallocate(problem%eq_a)
    if (allocated(problem%eq_b)) deallocate(problem%eq_b)
    allocate(problem%eq_a(size(a, 1), size(a, 2)), problem%eq_b(size(b)))
    problem%eq_a = a
    problem%eq_b = b
    status = 0
  end subroutine set_equality_constraints

  subroutine set_inequality_constraints(problem, a, b, status)
    type(maxlik_problem), intent(inout) :: problem
    real(dp), intent(in) :: a(:, :), b(:)
    integer, intent(out) :: status
    status = MAXLIK_INVALID_INPUT
    if (size(a, 2) /= problem%npar .or. size(a, 1) /= size(b)) return
    if (allocated(problem%ineq_a)) deallocate(problem%ineq_a)
    if (allocated(problem%ineq_b)) deallocate(problem%ineq_b)
    allocate(problem%ineq_a(size(a, 1), size(a, 2)), problem%ineq_b(size(b)))
    problem%ineq_a = a
    problem%ineq_b = b
    status = 0
  end subroutine set_inequality_constraints

  subroutine clear_constraints(problem)
    type(maxlik_problem), intent(inout) :: problem
    if (allocated(problem%eq_a)) deallocate(problem%eq_a)
    if (allocated(problem%eq_b)) deallocate(problem%eq_b)
    if (allocated(problem%ineq_a)) deallocate(problem%ineq_a)
    if (allocated(problem%ineq_b)) deallocate(problem%ineq_b)
    problem%penalty_rho = 0.0_dp
  end subroutine clear_constraints

  subroutine evaluate_value(problem, x, value, count, status, include_penalty)
    type(maxlik_problem), intent(in) :: problem
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: value
    integer, intent(inout) :: count
    integer, intent(out) :: status
    logical, intent(in), optional :: include_penalty

    logical :: use_penalty
    real(dp) :: penalty

    value = -huge(1.0_dp)
    status = MAXLIK_INVALID_INPUT
    if (.not. associated(problem%objective) .or. size(x) /= problem%npar) return
    call problem%objective(x, value, status)
    count = count + 1
    if (status /= 0 .or. .not. ieee_is_finite(value)) then
      status = MAXLIK_EVALUATION_ERROR
      return
    end if
    use_penalty = .true.
    if (present(include_penalty)) use_penalty = include_penalty
    if (use_penalty .and. problem%penalty_rho > 0.0_dp) then
      call penalty_value(problem, x, penalty)
      value = value - problem%penalty_rho * penalty
    end if
    status = 0
  end subroutine evaluate_value

  subroutine evaluate_gradient(problem, x, gradient, function_count, gradient_count, status, central, include_penalty)
    type(maxlik_problem), intent(in) :: problem
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: gradient(:)
    integer, intent(inout) :: function_count, gradient_count
    integer, intent(out) :: status
    logical, intent(in), optional :: central, include_penalty

    real(dp), allocatable :: scores(:, :), pg(:)
    logical :: use_penalty

    status = MAXLIK_INVALID_INPUT
    gradient = 0.0_dp
    if (size(gradient) /= problem%npar) return
    if (associated(problem%gradient)) then
      call problem%gradient(x, gradient, status)
      gradient_count = gradient_count + 1
    else if (associated(problem%scores) .and. problem%nobs > 0) then
      allocate(scores(problem%nobs, problem%npar))
      call problem%scores(x, scores, status)
      gradient_count = gradient_count + 1
      if (status == 0) gradient = sum(scores, dim=1)
    else
      call numeric_gradient(problem, x, gradient, function_count, status, central)
      gradient_count = gradient_count + 1
      if (status /= 0) return
      gradient = merge(gradient, 0.0_dp, problem%active)
      return
    end if
    if (status /= 0 .or. .not. all(ieee_is_finite(gradient))) then
      status = MAXLIK_EVALUATION_ERROR
      return
    end if
    use_penalty = .true.
    if (present(include_penalty)) use_penalty = include_penalty
    if (use_penalty .and. problem%penalty_rho > 0.0_dp) then
      allocate(pg(problem%npar))
      call penalty_gradient(problem, x, pg)
      gradient = gradient - problem%penalty_rho * pg
    end if
    where (.not. problem%active) gradient = 0.0_dp
    status = 0
  end subroutine evaluate_gradient

  subroutine evaluate_hessian(problem, x, hessian, function_count, gradient_count, hessian_count, status, include_penalty)
    type(maxlik_problem), intent(in) :: problem
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: hessian(:, :)
    integer, intent(inout) :: function_count, gradient_count, hessian_count
    integer, intent(out) :: status
    logical, intent(in), optional :: include_penalty

    real(dp), allocatable :: ph(:, :)
    logical :: use_penalty
    integer :: i

    hessian = 0.0_dp
    status = MAXLIK_INVALID_INPUT
    if (size(hessian, 1) /= problem%npar .or. size(hessian, 2) /= problem%npar) return
    if (associated(problem%hessian)) then
      call problem%hessian(x, hessian, status)
      hessian_count = hessian_count + 1
    else
      call numeric_hessian(problem, x, hessian, function_count, gradient_count, status)
      hessian_count = hessian_count + 1
      if (status /= 0) return
      return
    end if
    if (status /= 0 .or. .not. all(ieee_is_finite(hessian))) then
      status = MAXLIK_EVALUATION_ERROR
      return
    end if
    use_penalty = .true.
    if (present(include_penalty)) use_penalty = include_penalty
    if (use_penalty .and. problem%penalty_rho > 0.0_dp) then
      allocate(ph(problem%npar, problem%npar))
      call penalty_hessian(problem, x, ph)
      hessian = hessian - problem%penalty_rho * ph
    end if
    hessian = 0.5_dp * (hessian + transpose(hessian))
    do i = 1, problem%npar
      if (.not. problem%active(i)) then
        hessian(i, :) = 0.0_dp
        hessian(:, i) = 0.0_dp
        hessian(i, i) = -1.0_dp
      end if
    end do
    status = 0
  end subroutine evaluate_hessian

  subroutine evaluate_scores(problem, x, scores, gradient_count, status)
    type(maxlik_problem), intent(in) :: problem
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: scores(:, :)
    integer, intent(inout) :: gradient_count
    integer, intent(out) :: status

    status = MAXLIK_INVALID_INPUT
    scores = 0.0_dp
    if (.not. associated(problem%scores)) return
    if (size(scores, 1) /= problem%nobs .or. size(scores, 2) /= problem%npar) return
    call problem%scores(x, scores, status)
    gradient_count = gradient_count + 1
    if (status /= 0 .or. .not. all(ieee_is_finite(scores))) then
      status = MAXLIK_EVALUATION_ERROR
      return
    end if
    where (spread(.not. problem%active, 1, problem%nobs)) scores = 0.0_dp
    status = 0
  end subroutine evaluate_scores

  subroutine numeric_gradient(problem, x, gradient, function_count, status, central)
    type(maxlik_problem), intent(in) :: problem
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: gradient(:)
    integer, intent(inout) :: function_count
    integer, intent(out) :: status
    logical, intent(in), optional :: central

    real(dp) :: xp(size(x)), xm(size(x)), fp, fm, f0, step
    logical :: use_central
    integer :: i

    gradient = 0.0_dp
    status = 0
    use_central = .true.
    if (present(central)) use_central = central
    if (.not. use_central) then
      call evaluate_value(problem, x, f0, function_count, status)
      if (status /= 0) return
    end if
    do i = 1, size(x)
      if (.not. problem%active(i)) cycle
      if (use_central) then
        step = epsilon(1.0_dp)**(1.0_dp / 3.0_dp) * max(1.0_dp, abs(x(i)))
        xp = x
        xm = x
        xp(i) = min(problem%upper(i), x(i) + step)
        xm(i) = max(problem%lower(i), x(i) - step)
        if (abs(xp(i) - xm(i)) <= tiny(1.0_dp)) cycle
        call evaluate_value(problem, xp, fp, function_count, status)
        if (status /= 0) return
        call evaluate_value(problem, xm, fm, function_count, status)
        if (status /= 0) return
        gradient(i) = (fp - fm) / (xp(i) - xm(i))
      else
        step = sqrt(epsilon(1.0_dp)) * max(1.0_dp, abs(x(i)))
        xp = x
        xp(i) = min(problem%upper(i), x(i) + step)
        if (abs(xp(i) - x(i)) <= tiny(1.0_dp)) xp(i) = max(problem%lower(i), x(i) - step)
        if (abs(xp(i) - x(i)) <= tiny(1.0_dp)) cycle
        call evaluate_value(problem, xp, fp, function_count, status)
        if (status /= 0) return
        gradient(i) = (fp - f0) / (xp(i) - x(i))
      end if
    end do
  end subroutine numeric_gradient

  subroutine numeric_hessian(problem, x, hessian, function_count, gradient_count, status)
    type(maxlik_problem), intent(in) :: problem
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: hessian(:, :)
    integer, intent(inout) :: function_count, gradient_count
    integer, intent(out) :: status

    real(dp) :: xp(size(x)), xm(size(x)), gp(size(x)), gm(size(x)), step
    integer :: i, j

    hessian = 0.0_dp
    status = 0
    do j = 1, size(x)
      if (.not. problem%active(j)) then
        hessian(j, j) = -1.0_dp
        cycle
      end if
      step = epsilon(1.0_dp)**0.25_dp * max(1.0_dp, abs(x(j)))
      xp = x
      xm = x
      xp(j) = min(problem%upper(j), x(j) + step)
      xm(j) = max(problem%lower(j), x(j) - step)
      if (abs(xp(j) - xm(j)) <= tiny(1.0_dp)) cycle
      call evaluate_gradient(problem, xp, gp, function_count, gradient_count, status)
      if (status /= 0) return
      call evaluate_gradient(problem, xm, gm, function_count, gradient_count, status)
      if (status /= 0) return
      hessian(:, j) = (gp - gm) / (xp(j) - xm(j))
    end do
    hessian = 0.5_dp * (hessian + transpose(hessian))
    do i = 1, size(x)
      if (.not. problem%active(i)) then
        hessian(i, :) = 0.0_dp
        hessian(:, i) = 0.0_dp
        hessian(i, i) = -1.0_dp
      end if
    end do
  end subroutine numeric_hessian

  real(dp) function constraint_violation(problem, x) result(value)
    type(maxlik_problem), intent(in) :: problem
    real(dp), intent(in) :: x(:)
    real(dp), allocatable :: residual(:)

    value = 0.0_dp
    if (allocated(problem%eq_a)) then
      residual = matmul(problem%eq_a, x) + problem%eq_b
      value = max(value, maxval(abs(residual)))
    end if
    if (allocated(problem%ineq_a)) then
      residual = matmul(problem%ineq_a, x) + problem%ineq_b
      value = max(value, maxval(max(0.0_dp, -residual)))
    end if
  end function constraint_violation

  subroutine penalty_value(problem, x, value)
    type(maxlik_problem), intent(in) :: problem
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: value
    real(dp), allocatable :: residual(:)

    value = 0.0_dp
    if (allocated(problem%eq_a)) then
      residual = matmul(problem%eq_a, x) + problem%eq_b
      value = value + dot_product(residual, residual)
    end if
    if (allocated(problem%ineq_a)) then
      residual = min(matmul(problem%ineq_a, x) + problem%ineq_b, 0.0_dp)
      value = value + dot_product(residual, residual)
    end if
  end subroutine penalty_value

  subroutine penalty_gradient(problem, x, gradient)
    type(maxlik_problem), intent(in) :: problem
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: gradient(:)
    real(dp), allocatable :: residual(:)

    gradient = 0.0_dp
    if (allocated(problem%eq_a)) then
      residual = matmul(problem%eq_a, x) + problem%eq_b
      gradient = gradient + 2.0_dp * matmul(transpose(problem%eq_a), residual)
    end if
    if (allocated(problem%ineq_a)) then
      residual = min(matmul(problem%ineq_a, x) + problem%ineq_b, 0.0_dp)
      gradient = gradient + 2.0_dp * matmul(transpose(problem%ineq_a), residual)
    end if
  end subroutine penalty_gradient

  subroutine penalty_hessian(problem, x, hessian)
    type(maxlik_problem), intent(in) :: problem
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: hessian(:, :)
    real(dp), allocatable :: residual(:), active_a(:, :)
    integer :: i

    hessian = 0.0_dp
    if (allocated(problem%eq_a)) then
      hessian = hessian + 2.0_dp * matmul(transpose(problem%eq_a), problem%eq_a)
    end if
    if (allocated(problem%ineq_a)) then
      residual = matmul(problem%ineq_a, x) + problem%ineq_b
      allocate(active_a(size(problem%ineq_a, 1), problem%npar))
      active_a = 0.0_dp
      do i = 1, size(residual)
        if (residual(i) < 0.0_dp) active_a(i, :) = problem%ineq_a(i, :)
      end do
      hessian = hessian + 2.0_dp * matmul(transpose(active_a), active_a)
    end if
  end subroutine penalty_hessian

end module maxlik_evaluation
