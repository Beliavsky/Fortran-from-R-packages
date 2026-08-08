! SPDX-License-Identifier: GPL-2.0-only
module calibrar_gradient
  use calibrar_kinds, only : dp
  use calibrar_interfaces, only : scalar_objective
  implicit none
  private
  public :: gradient_options, numerical_gradient, gradient_forward, gradient_backward
  public :: gradient_central, gradient_richardson

  type :: gradient_options
    real(dp) :: eps = 1.0e-8_dp
    real(dp) :: d = 1.0e-4_dp
    real(dp) :: zero_tol = sqrt(epsilon(1.0_dp)/7.0e-7_dp)
    integer :: r = 4
    real(dp) :: v = 2.0_dp
  end type gradient_options

contains

  subroutine numerical_gradient(fn, x, g, method, options)
    procedure(scalar_objective) :: fn
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: g(:)
    character(len=*), intent(in), optional :: method
    type(gradient_options), intent(in), optional :: options
    character(len=16) :: m
    m = "richardson"
    if (present(method)) m = adjustl(method)
    select case (trim(m))
    case ("forward")
      call gradient_forward(fn, x, g, options)
    case ("backward")
      call gradient_backward(fn, x, g, options)
    case ("central")
      call gradient_central(fn, x, g, options)
    case ("richardson")
      call gradient_richardson(fn, x, g, options)
    case default
      error stop "numerical_gradient: undefined method"
    end select
  end subroutine numerical_gradient

  subroutine gradient_forward(fn, x, g, options)
    procedure(scalar_objective) :: fn
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: g(:)
    type(gradient_options), intent(in), optional :: options
    type(gradient_options) :: op
    real(dp), allocatable :: xx(:), h(:)
    real(dp) :: fx
    integer :: i
    op = gradient_options()
    if (present(options)) op = options
    allocate(xx(size(x)), h(size(x)))
    call simple_h(x, op, h)
    fx = eval_scalar(fn, x)
    do i = 1, size(x)
      xx = x; xx(i) = xx(i) + h(i)
      g(i) = (eval_scalar(fn, xx)-fx)/h(i)
    end do
  end subroutine gradient_forward

  subroutine gradient_backward(fn, x, g, options)
    procedure(scalar_objective) :: fn
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: g(:)
    type(gradient_options), intent(in), optional :: options
    type(gradient_options) :: op
    real(dp), allocatable :: xx(:), h(:)
    real(dp) :: fx
    integer :: i
    op = gradient_options()
    if (present(options)) op = options
    allocate(xx(size(x)), h(size(x)))
    call simple_h(x, op, h)
    fx = eval_scalar(fn, x)
    do i = 1, size(x)
      xx = x; xx(i) = xx(i) - h(i)
      g(i) = (fx-eval_scalar(fn, xx))/h(i)
    end do
  end subroutine gradient_backward

  subroutine gradient_central(fn, x, g, options)
    procedure(scalar_objective) :: fn
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: g(:)
    type(gradient_options), intent(in), optional :: options
    type(gradient_options) :: op
    real(dp), allocatable :: xp(:), xm(:), h(:)
    integer :: i
    op = gradient_options()
    if (present(options)) op = options
    allocate(xp(size(x)), xm(size(x)), h(size(x)))
    call simple_h(x, op, h)
    do i = 1, size(x)
      xp = x; xm = x
      xp(i) = xp(i) + h(i); xm(i) = xm(i) - h(i)
      g(i) = (eval_scalar(fn, xp)-eval_scalar(fn, xm))/(2.0_dp*h(i))
    end do
  end subroutine gradient_central

  subroutine gradient_richardson(fn, x, g, options)
    procedure(scalar_objective) :: fn
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: g(:)
    type(gradient_options), intent(in), optional :: options
    type(gradient_options) :: op
    real(dp), allocatable :: xp(:), xm(:), h(:), df(:,:)
    real(dp) :: fac
    integer :: i, k, m, rr
    op = gradient_options()
    if (present(options)) op = options
    rr = max(1, op%r)
    allocate(xp(size(x)), xm(size(x)), h(size(x)), df(rr,size(x)))
    call start_h(x, op, h)
    do k = 1, rr
      do i = 1, size(x)
        if (k > 1) then
          if (abs(df(k-1,i)) < 1.0e-20_dp) then
            df(k,i) = 0.0_dp
            cycle
          end if
        end if
        xp = x
        xm = x
        xp(i) = xp(i) + h(i)
        xm(i) = xm(i) - h(i)
        df(k,i) = (eval_scalar(fn, xp)-eval_scalar(fn, xm))/(2.0_dp*h(i))
      end do
      h = h/op%v
    end do
    do m = 1, rr-1
      fac = 4.0_dp**m
      do k = 1, rr-m
        df(k,:) = (fac*df(k+1,:)-df(k,:))/(fac-1.0_dp)
      end do
    end do
    g = df(1,:)
  end subroutine gradient_richardson

  subroutine simple_h(x, op, h)
    real(dp), intent(in) :: x(:)
    type(gradient_options), intent(in) :: op
    real(dp), intent(out) :: h(:)
    integer :: i
    do i = 1, size(x)
      if (abs(x(i)) < op%zero_tol) then
        h(i) = op%eps
      else
        h(i) = op%eps*abs(x(i))
      end if
      if (h(i) <= tiny(1.0_dp)) h(i) = op%eps
    end do
  end subroutine simple_h

  subroutine start_h(x, op, h)
    real(dp), intent(in) :: x(:)
    type(gradient_options), intent(in) :: op
    real(dp), intent(out) :: h(:)
    integer :: i
    do i = 1, size(x)
      if (abs(x(i)) < op%zero_tol) then
        h(i) = max(op%eps, 1.0e-4_dp)
      else
        h(i) = op%d*abs(x(i))
      end if
    end do
  end subroutine start_h

  function eval_scalar(fn, x) result(f)
    procedure(scalar_objective) :: fn
    real(dp), intent(in) :: x(:)
    real(dp) :: f
    f = fn(x)
  end function eval_scalar
end module calibrar_gradient
