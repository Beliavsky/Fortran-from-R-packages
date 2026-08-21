! SPDX-License-Identifier: GPL-2.0-or-later
module rootsolve_derivatives
  use rootsolve_kinds, only : dp
  use rootsolve_types, only : root_func, scalar_objective, steady_rhs
  implicit none
  private
  public :: perturb_value, gradient, hessian, jacobian_full, jacobian_band
contains

  pure real(dp) function perturb_value(value, pert) result(delta)
    real(dp), intent(in) :: value
    real(dp), intent(in), optional :: pert
    real(dp) :: p
    p = 1.0e-8_dp
    if (present(pert)) p = pert
    delta = max(abs(value) * p, p)
  end function perturb_value

  subroutine gradient(f, x, jac, centered, pert)
    procedure(root_func) :: f
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: jac(:,:)
    logical, intent(in), optional :: centered
    real(dp), intent(in), optional :: pert
    real(dp), allocatable :: xx(:), f0(:), fp(:), fm(:)
    real(dp) :: d
    logical :: ctr
    integer :: j, m, n

    n = size(x)
    m = size(jac,1)
    if (size(jac,2) /= n) error stop 'gradient: shape mismatch'
    allocate(xx(n), f0(m), fp(m), fm(m))
    xx = x
    call f(x, f0)
    ctr = .false.
    if (present(centered)) ctr = centered
    do j = 1, n
      d = perturb_value(x(j), pert)
      xx = x
      xx(j) = x(j) + d
      call f(xx, fp)
      if (ctr) then
        xx(j) = x(j) - d
        call f(xx, fm)
        jac(:,j) = (fp - fm) / (2.0_dp * d)
      else
        jac(:,j) = (fp - f0) / d
      end if
    end do
  end subroutine gradient

  subroutine hessian(f, x, hess, centered, pert)
    procedure(scalar_objective) :: f
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: hess(:,:)
    logical, intent(in), optional :: centered
    real(dp), intent(in), optional :: pert
    real(dp), allocatable :: xx(:), gp(:), gm(:), g0(:)
    real(dp) :: d
    logical :: ctr
    integer :: i, j, n

    n = size(x)
    if (size(hess,1) /= n .or. size(hess,2) /= n) error stop 'hessian: shape mismatch'
    allocate(xx(n), gp(n), gm(n), g0(n))
    ctr = .false.
    if (present(centered)) ctr = centered
    call scalar_gradient(f, x, g0, ctr, pert)
    do j = 1, n
      d = perturb_value(x(j), pert)
      xx = x
      xx(j) = x(j) + d
      call scalar_gradient(f, xx, gp, ctr, pert)
      if (ctr) then
        xx(j) = x(j) - d
        call scalar_gradient(f, xx, gm, ctr, pert)
        hess(:,j) = (gp - gm) / (2.0_dp * d)
      else
        hess(:,j) = (gp - g0) / d
      end if
    end do
    do j = 1, n
      do i = j+1, n
        d = 0.5_dp * (hess(i,j) + hess(j,i))
        hess(i,j) = d
        hess(j,i) = d
      end do
    end do
  end subroutine hessian

  subroutine scalar_gradient(f, x, g, centered, pert)
    procedure(scalar_objective) :: f
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: g(:)
    logical, intent(in) :: centered
    real(dp), intent(in), optional :: pert
    real(dp), allocatable :: xx(:)
    real(dp) :: f0, fp, fm, d
    integer :: j
    allocate(xx(size(x)))
    f0 = f(x)
    do j = 1, size(x)
      d = perturb_value(x(j), pert)
      xx = x
      xx(j) = x(j) + d
      fp = f(xx)
      if (centered) then
        xx(j) = x(j) - d
        fm = f(xx)
        g(j) = (fp - fm) / (2.0_dp * d)
      else
        g(j) = (fp - f0) / d
      end if
    end do
  end subroutine scalar_gradient

  subroutine jacobian_full(y, func, jac, time, dy, pert, centered)
    real(dp), intent(in) :: y(:)
    procedure(steady_rhs) :: func
    real(dp), intent(out) :: jac(:,:)
    real(dp), intent(in), optional :: time
    real(dp), intent(in), optional :: dy(:), pert
    logical, intent(in), optional :: centered
    real(dp), allocatable :: yy(:), ref(:), fp(:), fm(:)
    real(dp) :: t, d
    logical :: ctr
    integer :: j, n

    n = size(y)
    if (size(jac,1) /= n .or. size(jac,2) /= n) error stop 'jacobian_full: shape mismatch'
    allocate(yy(n), ref(n), fp(n), fm(n))
    t = 0.0_dp
    if (present(time)) t = time
    if (present(dy)) then
      if (size(dy) /= n) error stop 'jacobian_full: dy shape mismatch'
      ref = dy
    else
      call func(t, y, ref)
    end if
    ctr = .false.
    if (present(centered)) ctr = centered
    do j = 1, n
      d = perturb_value(y(j), pert)
      yy = y
      yy(j) = y(j) + d
      call func(t, yy, fp)
      if (ctr) then
        yy(j) = y(j) - d
        call func(t, yy, fm)
        jac(:,j) = (fp - fm) / (2.0_dp * d)
      else
        jac(:,j) = (fp - ref) / d
      end if
    end do
  end subroutine jacobian_full

  subroutine jacobian_band(y, func, bandup, banddown, jac, time, dy, pert)
    real(dp), intent(in) :: y(:)
    procedure(steady_rhs) :: func
    integer, intent(in) :: bandup, banddown
    real(dp), intent(out) :: jac(:,:)
    real(dp), intent(in), optional :: time
    real(dp), intent(in), optional :: dy(:), pert
    real(dp), allocatable :: yy(:), ref(:), fp(:), delta(:)
    real(dp) :: t
    integer :: j, k, i, first, last, n, nband

    n = size(y)
    nband = bandup + banddown + 1
    if (size(jac,1) /= nband .or. size(jac,2) /= n) error stop 'jacobian_band: shape mismatch'
    allocate(yy(n), ref(n), fp(n), delta(n))
    t = 0.0_dp
    if (present(time)) t = time
    if (present(dy)) then
      if (size(dy) /= n) error stop 'jacobian_band: dy shape mismatch'
      ref = dy
    else
      call func(t, y, ref)
    end if
    do i = 1, n
      delta(i) = perturb_value(y(i), pert)
    end do
    jac = 0.0_dp
    do j = 1, min(n, nband)
      yy = y
      do k = j, n, nband
        yy(k) = yy(k) + delta(k)
      end do
      call func(t, yy, fp)
      do k = j, n, nband
        first = max(1, k-bandup)
        last = min(n, k+banddown)
        do i = first, last
          jac(i-k+bandup+1,k) = (fp(i) - ref(i)) / delta(k)
        end do
      end do
    end do
  end subroutine jacobian_band
end module rootsolve_derivatives
