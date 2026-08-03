! SPDX-License-Identifier: GPL-2.0-only
module optimx_eval
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use optimx_kinds, only: dp
  use optimx_types
  implicit none
  private
  public :: initialize_problem, valid_problem, project_bounds, projected_gradient
  public :: evaluate_value, evaluate_gradient, evaluate_hessian
  public :: grfwd, grback, grcentral, grnd, grpracma

contains
  subroutine initialize_problem(problem, n)
    type(optimx_problem), intent(inout) :: problem
    integer, intent(in) :: n
    problem%n = n
    if (allocated(problem%lower)) deallocate(problem%lower)
    if (allocated(problem%upper)) deallocate(problem%upper)
    if (allocated(problem%mask)) deallocate(problem%mask)
    allocate(problem%lower(n), problem%upper(n), problem%mask(n))
    problem%lower = -huge(1.0_dp)
    problem%upper = huge(1.0_dp)
    problem%mask = 1
  end subroutine initialize_problem

  logical function valid_problem(problem, x) result(ok)
    type(optimx_problem), intent(in) :: problem
    real(dp), intent(in) :: x(:)
    ok = problem%n == size(x) .and. associated(problem%objective)
    if (.not. ok) return
    ok = allocated(problem%lower) .and. allocated(problem%upper) .and. allocated(problem%mask)
    if (.not. ok) return
    ok = size(problem%lower) == problem%n .and. size(problem%upper) == problem%n
    ok = ok .and. size(problem%mask) == problem%n
    if (.not. ok) return
    ok = all(problem%lower <= problem%upper) .and. all(ieee_is_finite(x))
  end function valid_problem

  pure subroutine project_bounds(problem, x)
    type(optimx_problem), intent(in) :: problem
    real(dp), intent(inout) :: x(:)
    integer :: i
    do i = 1, size(x)
      if (problem%mask(i) == 0) cycle
      x(i) = max(problem%lower(i), min(problem%upper(i), x(i)))
    end do
  end subroutine project_bounds

  pure subroutine projected_gradient(problem, x, g, pg)
    type(optimx_problem), intent(in) :: problem
    real(dp), intent(in) :: x(:), g(:)
    real(dp), intent(out) :: pg(:)
    integer :: i
    pg = g
    do i = 1, size(x)
      if (problem%mask(i) == 0) then
        pg(i) = 0.0_dp
      else if (x(i) <= problem%lower(i) + 10.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(x(i)))) then
        if (g(i) > 0.0_dp) pg(i) = 0.0_dp
      else if (x(i) >= problem%upper(i) - 10.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(x(i)))) then
        if (g(i) < 0.0_dp) pg(i) = 0.0_dp
      end if
    end do
  end subroutine projected_gradient

  subroutine evaluate_value(problem, x, value, count, status)
    type(optimx_problem), intent(in) :: problem
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: value
    integer, intent(inout) :: count
    integer, intent(out) :: status
    real(dp) :: g(size(x)), h(size(x),size(x))
    g = 0.0_dp; h = 0.0_dp
    call problem%objective(x, value, g, h, .false., .false., status)
    count = count + 1
    if (status == 0 .and. .not. ieee_is_finite(value)) status = OPTIMX_BAD_EVALUATION
  end subroutine evaluate_value

  subroutine evaluate_gradient(problem, x, g, fcount, gcount, status, central)
    type(optimx_problem), intent(in) :: problem
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: g(:)
    integer, intent(inout) :: fcount, gcount
    integer, intent(out) :: status
    logical, intent(in), optional :: central
    real(dp) :: f, hess(size(x),size(x))
    logical :: use_central
    use_central = .true.; if (present(central)) use_central = central
    if (problem%has_gradient) then
      hess = 0.0_dp; g = 0.0_dp
      call problem%objective(x, f, g, hess, .true., .false., status)
      gcount = gcount + 1
      if (status == 0 .and. .not. all(ieee_is_finite(g))) status = OPTIMX_BAD_EVALUATION
    else
      if (use_central) then
        call grcentral(problem, x, g, fcount, status)
      else
        call grfwd(problem, x, g, fcount, status)
      end if
      gcount = gcount + 1
    end if
  end subroutine evaluate_gradient

  subroutine evaluate_hessian(problem, x, h, fcount, gcount, hcount, status)
    type(optimx_problem), intent(in) :: problem
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: h(:, :)
    integer, intent(inout) :: fcount, gcount, hcount
    integer, intent(out) :: status
    real(dp) :: f, g(size(x)), gp(size(x)), gm(size(x)), xp(size(x)), xm(size(x)), step
    integer :: j
    if (problem%has_hessian) then
      g = 0.0_dp; h = 0.0_dp
      call problem%objective(x, f, g, h, .true., .true., status)
      gcount = gcount + 1; hcount = hcount + 1
      if (status == 0 .and. .not. all(ieee_is_finite(h))) status = OPTIMX_BAD_EVALUATION
      return
    end if
    status = 0
    do j = 1, size(x)
      step = epsilon(1.0_dp)**(1.0_dp/3.0_dp) * max(1.0_dp, abs(x(j)))
      xp = x; xm = x
      xp(j) = min(problem%upper(j), x(j) + step)
      xm(j) = max(problem%lower(j), x(j) - step)
      if (abs(xp(j) - xm(j)) <= tiny(1.0_dp)) then
        h(:,j) = 0.0_dp
      else
        call evaluate_gradient(problem, xp, gp, fcount, gcount, status)
        if (status /= 0) return
        call evaluate_gradient(problem, xm, gm, fcount, gcount, status)
        if (status /= 0) return
        h(:,j) = (gp - gm) / (xp(j) - xm(j))
      end if
    end do
    h = 0.5_dp * (h + transpose(h))
    hcount = hcount + 1
  end subroutine evaluate_hessian

  subroutine grfwd(problem, x, g, fcount, status)
    type(optimx_problem), intent(in) :: problem
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: g(:)
    integer, intent(inout) :: fcount
    integer, intent(out) :: status
    real(dp) :: xp(size(x)), f0, fp, step
    integer :: i
    call evaluate_value(problem, x, f0, fcount, status); if (status /= 0) return
    do i = 1, size(x)
      if (problem%mask(i) == 0) then; g(i)=0.0_dp; cycle; end if
      step = sqrt(epsilon(1.0_dp))*max(1.0_dp,abs(x(i)))
      xp = x; xp(i) = min(problem%upper(i), x(i)+step)
      if (abs(xp(i) - x(i)) <= tiny(1.0_dp)) then
        xp(i) = max(problem%lower(i), x(i)-step)
      end if
      call evaluate_value(problem, xp, fp, fcount, status); if (status /= 0) return
      g(i) = (fp-f0)/(xp(i)-x(i))
    end do
  end subroutine grfwd

  subroutine grback(problem, x, g, fcount, status)
    type(optimx_problem), intent(in) :: problem
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: g(:)
    integer, intent(inout) :: fcount
    integer, intent(out) :: status
    real(dp) :: xm(size(x)), f0, fm, step
    integer :: i
    call evaluate_value(problem, x, f0, fcount, status); if (status /= 0) return
    do i = 1, size(x)
      if (problem%mask(i) == 0) then; g(i)=0.0_dp; cycle; end if
      step = sqrt(epsilon(1.0_dp))*max(1.0_dp,abs(x(i)))
      xm = x; xm(i) = max(problem%lower(i), x(i)-step)
      if (abs(xm(i) - x(i)) <= tiny(1.0_dp)) xm(i) = min(problem%upper(i), x(i)+step)
      call evaluate_value(problem, xm, fm, fcount, status); if (status /= 0) return
      g(i) = (f0-fm)/(x(i)-xm(i))
    end do
  end subroutine grback

  subroutine grcentral(problem, x, g, fcount, status)
    type(optimx_problem), intent(in) :: problem
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: g(:)
    integer, intent(inout) :: fcount
    integer, intent(out) :: status
    real(dp) :: xp(size(x)), xm(size(x)), fp, fm, step
    integer :: i
    status = 0
    do i = 1, size(x)
      if (problem%mask(i) == 0) then; g(i)=0.0_dp; cycle; end if
      step = epsilon(1.0_dp)**(1.0_dp/3.0_dp)*max(1.0_dp,abs(x(i)))
      xp = x; xm = x
      xp(i) = min(problem%upper(i), x(i)+step)
      xm(i) = max(problem%lower(i), x(i)-step)
      if (abs(xp(i) - xm(i)) <= tiny(1.0_dp)) then; g(i)=0.0_dp; cycle; end if
      call evaluate_value(problem, xp, fp, fcount, status); if (status /= 0) return
      call evaluate_value(problem, xm, fm, fcount, status); if (status /= 0) return
      g(i) = (fp-fm)/(xp(i)-xm(i))
    end do
  end subroutine grcentral

  subroutine grnd(problem, x, g, fcount, status)
    type(optimx_problem), intent(in) :: problem
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: g(:)
    integer, intent(inout) :: fcount
    integer, intent(out) :: status
    call grcentral(problem, x, g, fcount, status)
  end subroutine grnd

  subroutine grpracma(problem, x, g, fcount, status)
    type(optimx_problem), intent(in) :: problem
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: g(:)
    integer, intent(inout) :: fcount
    integer, intent(out) :: status
    call grcentral(problem, x, g, fcount, status)
  end subroutine grpracma
end module optimx_eval
