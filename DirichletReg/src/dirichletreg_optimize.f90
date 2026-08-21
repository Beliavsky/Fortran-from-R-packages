! SPDX-License-Identifier: GPL-2.0-or-later
module dirichletreg_optimize
  use dirichletreg_kinds, only : dp
  use dirichletreg_linalg, only : solve_linear
  implicit none
  private
  public :: maximize_bfgs, maximize_newton

  abstract interface
    subroutine objective_grad(theta, f, g)
      import dp
      real(dp), intent(in) :: theta(:)
      real(dp), intent(out) :: f
      real(dp), intent(out) :: g(:)
    end subroutine objective_grad

    subroutine objective_hess(theta, f, g, h)
      import dp
      real(dp), intent(in) :: theta(:)
      real(dp), intent(out) :: f
      real(dp), intent(out) :: g(:)
      real(dp), intent(out) :: h(:,:)
    end subroutine objective_hess
  end interface

contains

  subroutine maximize_bfgs(fn, theta, f, iterations, convergence, iterlim, tol)
    procedure(objective_grad) :: fn
    real(dp), intent(inout) :: theta(:)
    real(dp), intent(out) :: f
    integer, intent(out) :: iterations, convergence
    integer, intent(in), optional :: iterlim
    real(dp), intent(in), optional :: tol

    integer :: n, it, maxit
    real(dp) :: eps, fnew, step, gd, rho, relchg
    real(dp), allocatable :: g(:), gnew(:), p(:), s(:), y(:), h(:,:), v(:,:), ident(:,:), trial(:)

    n = size(theta)
    maxit = 10000
    if (present(iterlim)) maxit = iterlim
    eps = sqrt(epsilon(1.0_dp))
    if (present(tol)) eps = tol

    allocate(g(n), gnew(n), p(n), s(n), y(n), h(n,n), v(n,n), ident(n,n), trial(n))
    ident = 0.0_dp
    do it = 1, n
      ident(it,it) = 1.0_dp
    end do
    h = ident
    call fn(theta, f, g)
    convergence = 1

    do it = 1, maxit
      if (maxval(abs(g)) <= eps) then
        convergence = 0
        exit
      end if
      p = matmul(h,g)
      gd = dot_product(g,p)
      if (gd <= 0.0_dp) then
        p = g
        gd = dot_product(g,p)
        h = ident
      end if

      step = 1.0_dp
      do
        trial = theta + step*p
        call fn(trial, fnew, gnew)
        if (fnew >= f + 1.0e-4_dp*step*gd) exit
        step = 0.5_dp*step
        if (step < 1.0e-12_dp) exit
      end do
      if (step < 1.0e-12_dp) then
        convergence = 2
        exit
      end if

      s = trial - theta
      y = g - gnew ! gradient change for minimizing -f
      relchg = abs(fnew-f)/max(1.0_dp,abs(f))
      if (dot_product(y,s) > 1.0e-12_dp*sqrt(dot_product(y,y)*dot_product(s,s))) then
        rho = 1.0_dp/dot_product(y,s)
        v = ident - rho*outer(s,y)
        h = matmul(v,matmul(h,transpose(v))) + rho*outer(s,s)
      else
        h = ident
      end if
      theta = trial
      f = fnew
      g = gnew
      if (relchg <= eps .and. maxval(abs(g)) <= sqrt(eps)) then
        convergence = 0
        exit
      end if
    end do
    iterations = min(it,maxit)

  contains
    pure function outer(a,b) result(c)
      real(dp), intent(in) :: a(:), b(:)
      real(dp) :: c(size(a),size(b))
      integer :: i
      do i = 1, size(a)
        c(i,:) = a(i)*b
      end do
    end function outer
  end subroutine maximize_bfgs


  subroutine maximize_newton(fn, theta, f, hfinal, iterations, convergence, iterlim, tol)
    procedure(objective_hess) :: fn
    real(dp), intent(inout) :: theta(:)
    real(dp), intent(out) :: f
    real(dp), intent(out) :: hfinal(:,:)
    integer, intent(out) :: iterations, convergence
    integer, intent(in), optional :: iterlim
    real(dp), intent(in), optional :: tol

    integer :: n, it, maxit, ierr
    real(dp) :: eps, fnew, step, gd, relchg
    real(dp), allocatable :: g(:), gnew(:), h(:,:), hnew(:,:), p(:), trial(:)

    n = size(theta)
    maxit = 10000
    if (present(iterlim)) maxit = iterlim
    eps = epsilon(1.0_dp)**0.75_dp
    if (present(tol)) eps = tol
    allocate(g(n),gnew(n),h(n,n),hnew(n,n),p(n),trial(n))
    call fn(theta,f,g,h)
    convergence = 1

    do it = 1, maxit
      if (maxval(abs(g)) <= eps) then
        convergence = 0
        exit
      end if
      call solve_linear(h,-g,p,ierr)
      if (ierr /= 0 .or. dot_product(g,p) <= 0.0_dp) p = g/max(1.0_dp,maxval(abs(g)))
      gd = dot_product(g,p)
      step = 1.0_dp
      do
        trial = theta + step*p
        call fn(trial,fnew,gnew,hnew)
        if (fnew >= f + 1.0e-4_dp*step*gd) exit
        step = 0.5_dp*step
        if (step < 1.0e-12_dp) exit
      end do
      if (step < 1.0e-12_dp) then
        convergence = 2
        exit
      end if
      relchg = abs(fnew-f)/max(1.0_dp,abs(f))
      theta = trial
      f = fnew
      g = gnew
      h = hnew
      if (relchg <= eps .and. maxval(abs(g)) <= sqrt(eps)) then
        convergence = 0
        exit
      end if
    end do
    iterations = min(it,maxit)
    hfinal = h
  end subroutine maximize_newton

end module dirichletreg_optimize
